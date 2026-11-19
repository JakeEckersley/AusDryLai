################################################################################
#### 3c - assemble and filter the SL2P retrievals                           ####
################################################################################
# Every post-download QC step for the LEAF products, so no later script repeats
# it. Script 2f applies only what it can apply in Earth Engine.
#
# In:  <LeafDownloads_*>/stacks/<SITE>_<start>_<end>.tif  one per site-year
#      <LeafDownloads_*>/dates/<SITE>_<start>_<end>.csv   acquisition times, view
#                                                         zenith angles (2f)
# Out: <out_dir>/<SITE>/leaf_<tag>-ts.tif       estimateLAI, QC applied
#      <out_dir>/<SITE>/leaf_<tag>_err-ts.tif   errorLAI
#      <out_dir>/<SITE>/leaf_<tag>_qc-ts.tif    QC band, unapplied
#
# A site is skipped once all three outputs exist; delete them to rebuild.

library(terra)
source("R/stack_io.R")

## ---- config -----------------------------------------------------------------
leaf_dirs <- c(leaf_10m = "D:/PilbaraLaiProject/LeafDownloads_10m",
               leaf_20m = "D:/PilbaraLaiProject/LeafDownloads_20m")
out_dir   <- "D:/PilbaraLaiProject/FilteredLeaf"

# 1 = reflectance outside the domain SL2P was trained on. 0 applies nothing.
qc_bits <- 1L

# Keep one retrieval per date. Sentinel-2 delivers each overpass once per MGRS
# tile and SL2P runs on each separately, so a site near a boundary carries the
# same pass twice - about half of all layers at MS, YU, CW, HC and the MCKs.
# Most AOI pixels retrieved after the domain mask wins, the smaller view zenith
# angle breaks ties, as in script 3b.
one_pass_a_day <- TRUE

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- assemble ---------------------------------------------------------------
# GEE strips band names, so the sidecar carries the times. toBands() emits the
# bands per image in order - estimateLAI, errorLAI, QC, then the s2cloudless
# probability, which is not used here - so sidecar row i is layers nb*(i-1) + 1:3
# with nb = 3 or 4 bands per date.
read_site <- function(stacks) {
  lay <- list(); qly <- list(); ely <- list(); vza <- numeric(0)
  utc <- as.POSIXct(character(0), tz = "UTC")
  for (p in stacks) {
    d <- file.path(dirname(dirname(p)), "dates",
                   sub(".tif", ".csv", basename(p), fixed = TRUE))
    if (!file.exists(d)) {
      message("  no dates sidecar for ", basename(p), " - re-run script 2f"); next }
    dt <- read.csv(d, stringsAsFactors = FALSE)
    if (!"vza" %in% names(dt)) stop("no view zenith angles in ", basename(d), " - re-download with script 2f")
    x  <- try(rast(p), silent = TRUE)
    if (inherits(x, "try-error")) { message("  unreadable: ", basename(p)); next }
    nb <- nlyr(x) / nrow(dt)
    if (!nb %in% c(3, 4)) { message("  layer/date mismatch in ", basename(p)); next }
    i <- seq(1L, by = nb, length.out = nrow(dt))
    lay <- c(lay, list(x[[i]]))
    ely <- c(ely, list(x[[i + 1L]]))
    qly <- c(qly, list(x[[i + 2L]]))
    utc <- c(utc, as.POSIXct(dt$stamp, format = "%Y_%m_%d_%H_%M_%S", tz = "UTC"))
    vza <- c(vza, dt$vza)
  }
  if (!length(lay)) return(NULL)
  o <- order(utc)
  list(lai = do.call(c, lay)[[o]],
       err = if (length(ely) == length(lay)) do.call(c, ely)[[o]] else NULL,
       qc  = if (length(qly) == length(lay)) do.call(c, qly)[[o]] else NULL,
       vza = vza[o],
       utc = utc[o])
}

# Coverage is counted after the domain mask. It decides first, which settles the
# partial-clip case. Where both deliveries cover the site equally the retrievals
# still differ, because SL2P takes per-tile view geometry as an input, so the
# smaller view zenith angle breaks the tie. A missing angle loses to any value.
collapse_to_one_per_date <- function(x, label) {
  dte <- as.Date(format(x$utc, tz = "Australia/Perth"))
  if (!anyDuplicated(dte)) return(x)
  cov <- 1 - colMeans(is.na(values(x$lai)))
  vz  <- x$vza; vz[!is.finite(vz)] <- Inf
  sel <- sort(unname(vapply(split(seq_along(dte), dte), function(j)
    j[order(-cov[j], vz[j])][1], integer(1))))
  message("  collapsed ", nlyr(x$lai), " layers to ", length(sel), " for ", label)
  list(lai = x$lai[[sel]],
       err = if (is.null(x$err)) NULL else x$err[[sel]],
       qc  = if (is.null(x$qc))  NULL else x$qc[[sel]],
       vza = x$vza[sel],
       utc = x$utc[sel])
}

## ---- run --------------------------------------------------------------------
for (nm in names(leaf_dirs)) {
  stacks <- Sys.glob(file.path(leaf_dirs[[nm]], "stacks", "*.tif"))
  stacks <- stacks[!grepl(".part", stacks, fixed = TRUE)]
  if (!length(stacks)) { message(nm, ": no annual stacks, skipping"); next }
  stack_site <- sub("_.*", "", basename(stacks))

  for (s in sort(unique(stack_site))) {
    f_out <- file.path(out_dir, s, paste0(nm, "-ts.tif"))
    f_err <- sub("-ts.tif", "_err-ts.tif", f_out, fixed = TRUE)
    f_qc  <- sub("-ts.tif", "_qc-ts.tif",  f_out, fixed = TRUE)
    if (all(file.exists(f_out, f_err, f_qc))) next
    dir.create(dirname(f_out), showWarnings = FALSE, recursive = TRUE)

    x <- read_site(stacks[stack_site == s])
    if (is.null(x)) { message("  nothing assembled for ", s); next }
    # mask first, so the one-per-date choice sees the pixels that survive; the
    # QC band itself is still written unapplied below
    x$lai <- apply_qc_mask(x$lai, bits = qc_bits, qc = x$qc)
    if (one_pass_a_day) x <- collapse_to_one_per_date(x, paste(s, nm))

    nms <- format(x$utc, tz = "Australia/Perth", usetz = TRUE)
    if (!is.null(x$qc)) {
      names(x$qc) <- nms
      writeRaster(x$qc, f_qc, overwrite = TRUE)
    }
    lai <- x$lai
    names(lai) <- nms
    writeRaster(lai, f_out, overwrite = TRUE)
    if (!is.null(x$err)) {
      names(x$err) <- nms
      writeRaster(x$err, f_err, overwrite = TRUE)
    }
    message("filtered ", nm, " for ", s, ": ", nlyr(lai), " dates")
  }
}
