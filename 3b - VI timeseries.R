# 6 - create cleaned timeseries
library(terra)
library(future)
library(future.apply)
library(progressr)

# create output directory
out_dir <- "D:/PilbaraLaiProject/FliteredTimeseries"
dir.create(out_dir,showWarnings = F)

# One stack per site-year with dates and geometry in sidecars; R/stack_io.R does
# the indexing.  DEA NBART (script 2a)
# alt: point sat_root to D:/PilbaraLaiProject/S2_STACK_ROOT to run from the raw Sen2Cor reflectance of script 2b instead.
source("R/stack_io.R")
sat_root  <- "D:/PilbaraLaiProject/DEA_NBART_stacks"
stack_idx <- stack_index(sat_root)

# minimum fraction of the site AOI carrying valid, unflagged data on a date
min_valid_frac <- 0.80

# NDVI screen: a pixel-date is flagged when its NDVI lies more than mad_k MADs below (or above)
# the median of that pixel's dates within screen_days either side, itself included. The flag is
# ignored when fewer than screen_side dates lie on either side.
mad_k <- c(below = 2.5, above = 3)
screen_days <- 30
screen_side <- 2

# One observation per date. Sentinel-2 delivers each overpass once per MGRS
# tile, so a site near a boundary carries the same pass twice. Most pixels
# covered wins, lowest view angle breaks ties.
one_pass_a_day <- TRUE

# DEA masks to apply, from the bit flags script 2a writes to qc/.
#   cloud_dea   DEA's s2cloudless-mask, their own cutoff
#   cloud_prob  s2cloudless-prob above 3b's cloud_prob_thresh
#   fmask       fmask cloud or cloud shadow
#   contig      at least one band missing
#   terrain     terrain shadow (flags no pixel at these sites)
# DEA ships its masks undilated, so mask_dilate_m buffers cloud edges.
mask_apply <- c(cloud_dea = TRUE, cloud_prob = TRUE, fmask = TRUE, contig = TRUE, terrain = TRUE)
mask_dilate_m <- 50
if (is.null(stack_idx)) stop("no annual stacks in ", S2_STACK_ROOT)

site_list <- sort(unique(stack_idx$Site))

S2_metadata <- do.call(rbind, lapply(sort(unique(stack_idx$Site)), function(s) {
  m <- site_meta(s, stack_idx)
  if (is.null(m)) return(NULL)
  data.frame(Site = s, Date = m$Date, DateTime = format(m$DateTime, usetz = TRUE),
             vza = m$MEAN_INCIDENCE_ZENITH_ANGLE_B8A,   # used for the one_pass_a_day tie-break
             stringsAsFactors = FALSE)
}))

# range each index can physically take; values outside it (ratio indices with a
# near-zero denominator) are set to NA and counted in the log
vi_valid <- list(NDVI   = c(-1, 1),    EVI  = c(-1, 2.5), EVI2 = c(-1, 2.5),
                 kNDVI  = c(0, 1),     NDGI = c(-1, 1),   MSAVI2 = c(-1, 1),
                 NIRv   = c(-1, 1),    SR   = c(0, 30),   MSR  = c(-1, 10),
                 GCVI   = c(-1, 30),   DVI  = c(-1, 1))

# Departure of each date from the median of the pixel's dates within screen_days either side,
# in MADs, for one pixel's series x on dates d (days). NA where fewer than screen_side dates lie
# on either side.
zdev <- function(x, d) {
  v <- which(!is.na(x)); out <- rep(NA_real_, length(x)); m <- length(v)
  if (m < 2 * screen_side + 1) return(out)
  xv <- x[v]; dv <- d[v]
  lo <- findInterval(dv - screen_days - 1e-9, dv) + 1L; hi <- findInterval(dv + screen_days, dv)
  i <- seq_len(m); ok <- which(i - lo >= screen_side & hi - i >= screen_side)
  out[v[ok]] <- vapply(ok, function(j) {
    w <- xv[lo[j]:hi[j]]; s <- mad(w)
    if (s > 0) (xv[j] - median(w)) / s else NA_real_
  }, numeric(1))
  out
}

# create a log file and a function for writing to it - tried the progressr progress bar but wouldn't work for site-level mulitsessions
log_file <- file.path(out_dir, "qc_log.txt")
writeLines("", log_file)  # initialise empty

log_msg <- function(...) {
  msg <- paste0(paste(...), ': ', format(Sys.time(), "%d/%m/%y %X"), '\n')
  cat(msg, file = log_file, append = TRUE)
}

site = 'YU' # (for debugging)
PairedTsQc <- function(site){
  terra::terraOptions(memfrac = 0.15, progress = 0)
  source("R/stack_io.R")  # workers start without the parent session's functions
  log_msg('running ',site)

  s2md_sub <- S2_metadata[S2_metadata$Site==site,]

  log_msg('creating ',site,' directory: ',format(Sys.time(), "%d/%m/%y %X"))
  dir.create(file.path(out_dir,site),showWarnings = F)

  #########################################
  # pixelwise implementation
  #########################################
  # B4 across every year for this site, in date order (layers picked by name suffix)
  B4_stack  <- site_band_series(site, 'B4', stack_idx)
  stopifnot(!is.null(B4_stack))

  # the dates come from the meta sidecars and the pixels from the stacks, so
  # check they describe the same images before anything relies on the pairing
  if (nlyr(B4_stack) != nrow(s2md_sub))
    stop("site ", site, ": ", nlyr(B4_stack), " stack layers but ",
         nrow(s2md_sub), " metadata rows - a stack or sidecar is missing")
  # DEA masks, per mask_apply
  bad_mask <- NULL
  qc_rast  <- site_qc_series(site, stack_idx)
  if (is.null(qc_rast)) {
    log_msg('no QC stacks for ',site,' - re-run script 2a to apply the DEA masks')
  } else {
    stopifnot(nlyr(qc_rast) == nlyr(B4_stack))
    # bit values are set in script 2a; keep the two in step
    bits <- c(cloud_dea = 1L, cloud_prob = 2L, fmask = 4L,
              contig = 8L, terrain = 16L)[names(mask_apply)[mask_apply]]
    bits <- c(bits, nodata = 32L)
    # only cloud and shadow edges are dilated; the rest are hard edges
    cloudy <- c("cloud_dea","cloud_prob","fmask")
    soft <- sum(bits[names(bits) %in% cloudy])
    hard <- sum(bits[!names(bits) %in% cloudy])
    bad_mask <- if (soft > 0) bits_set(qc_rast, soft) else NULL
    if (!is.null(bad_mask) && mask_dilate_m > 0) {
      w <- 2 * ceiling(mask_dilate_m / res(B4_stack)[1]) + 1
      bad_mask <- focal(bad_mask, w = w, fun = "max", na.rm = TRUE) > 0
    }
    if (hard > 0) {
      h <- bits_set(qc_rast, hard)
      bad_mask <- if (is.null(bad_mask)) h else (bad_mask | h)
    }
    log_msg('DEA masks applied (',paste(names(bits), collapse = '+'),
            ', dilate ',mask_dilate_m,'m) ',site)
  }
  if (!is.null(bad_mask)) B4_stack <- mask(B4_stack, bad_mask, maskvalues = 1)

  # sel also subsets the VI stacks below, so every output shares a layer set
  sel <- NULL
  if (one_pass_a_day) {
    cov <- 1 - colMeans(is.na(values(B4_stack)))
    sel <- sort(unname(vapply(split(seq_along(cov), s2md_sub$Date), function(j)
      j[order(-cov[j], s2md_sub$vza[j])][1], integer(1))))
    if (length(sel) < nlyr(B4_stack)) {
      log_msg('collapsed ',nlyr(B4_stack),' layers to ',length(sel),
              ' one-per-date ',site)
      B4_stack <- B4_stack[[sel]]
      s2md_sub <- s2md_sub[sel, ]
    } else sel <- NULL
  }
  # the VI stacks are read fresh below, so the mask needs the same layer subset
  if (!is.null(bad_mask) && !is.null(sel)) bad_mask <- bad_mask[[sel]]

  names(B4_stack) <- as.character(s2md_sub$Date)
  log_msg('B4 stack loaded ',site)

  # NDVI after the masks and the range check as a matrix (pixels x dates), and the pixel-dates the screen flags
  ndvi_stack <- site_vi_series(site, 'NDVI', stack_idx, lyrs = sel)
  stopifnot(!is.null(ndvi_stack))
  if (!is.null(bad_mask)) ndvi_stack <- mask(ndvi_stack, bad_mask, maskvalues = 1)
  ndvi_mat <- values(ndvi_stack)
  ndvi_mat[!is.na(ndvi_mat) & (ndvi_mat < vi_valid$NDVI[1] | ndvi_mat > vi_valid$NDVI[2])] <- NA
  z <- t(apply(ndvi_mat, 1, zdev, d = as.numeric(s2md_sub$Date)))
  nd_qc <- !is.na(z) & (z < -mad_k[["below"]] | z > mad_k[["above"]])
  log_msg('NDVI screen flagged ',sum(nd_qc),' of ',sum(!is.na(ndvi_mat)),' pixel-dates ',site)

  # one layer filter shared by every band, so all 11 outputs carry the same
  # dates: at most 1 - min_valid_frac of the pixels masked, missing or flagged
  keep <- which(colMeans(is.na(ndvi_mat) | nd_qc) <= 1 - min_valid_frac)
  stopifnot(length(keep) > 0)
  log_msg('keeping ',length(keep),' of ',nlyr(B4_stack),' dates ',site)

  # back to raster
  nd_qc_rast <- B4_stack
  values(nd_qc_rast) <- nd_qc

  # tidy
  rm(ndvi_stack,ndvi_mat,z,nd_qc,B4_stack)
  gc()

  # ok - do i save this or create a new TS that isn't so slow. this method rocks way harder than the lists i'm using in timeseries analysis....

  # indices derived from the reflectance by the helper in R/stack_io.R
  band_names <- VI_NAMES

  for(band in band_names){

    log_msg('running mask for ',band,' ',site)

    vi_stack <- site_vi_series(site, band, stack_idx, lyrs = sel)
    stopifnot(!is.null(vi_stack))
    if (!is.null(bad_mask)) vi_stack <- mask(vi_stack, bad_mask, maskvalues = 1)
    names(vi_stack) <- as.character(s2md_sub$DateTime)

    # triple check the names align - before anything relies on them
    stopifnot(nlyr(vi_stack) == nlyr(nd_qc_rast))
    stopifnot(as.Date(names(vi_stack)) == as.Date(names(nd_qc_rast)))

    # mask the dodgy pixels, then apply the shared layer filter
    vi_stack <- mask(vi_stack,nd_qc_rast,maskvalues=1)
    vi_stack <- vi_stack[[keep]]

    # drop values the index cannot physically take
    rng <- vi_valid[[band]]
    if (!is.null(rng)) {
      n_bad <- sum(values(vi_stack) < rng[1] | values(vi_stack) > rng[2], na.rm = TRUE)
      if (n_bad > 0) {
        vi_stack <- clamp(vi_stack, rng[1], rng[2], values = FALSE)
        log_msg('flagged ', n_bad, ' ', band, ' values outside [', rng[1], ', ',
                rng[2], '] ', site)
      }
    }

    # write the raster
    writeRaster(vi_stack,
                filename = file.path(out_dir,site, paste0(band, "-ts.tif")),
                overwrite = TRUE)
    rm(vi_stack)
    gc()
  }

  return(paste0('Completed ',site))
  log_msg('Completed ',site)
}

# for(i in 1:length(site_list)){
#   PairedTsQc(site_list[i])
# }

log_file_norm <- normalizePath(file.path(out_dir, "qc_log.txt"), mustWork = FALSE)
# log window only in an interactive session (it hangs under Rscript)
if (interactive()) {
  message('Automatically opening power shell to get log updates:')
  try(shell(paste0('start "" powershell -NoExit -Command "Get-Content \\"',
                   log_file_norm, '\\" -Wait"'), wait = FALSE), silent = TRUE)
}

# capped at 8: a session per site exhausted memory, and the work is I/O bound
n_workers <- min(8L, parallelly::availableCores() - 2L, length(site_list))
log_msg('starting multisession ncores =', n_workers)
message("copy and paste this into powershell for logging updates if that didn't work: Get-Content ",log_file_norm," -Wait")

# ---------------------- start multisession -----------------------
plan(multisession, workers = n_workers)

future_lapply(site_list,
              PairedTsQc,
              future.globals = TRUE,
              future.packages = c("terra"))

plan(sequential)
# ---------------------- end multisession -----------------------
log_msg('timeseries cleaning completed')

