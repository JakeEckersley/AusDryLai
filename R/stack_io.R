################################################################################
#### Reading the annual stacks written by scripts 2a, 2b and 2g              ####
################################################################################
# Each exports one GeoTIFF per site-year. toBands() lays the collection out
# image-major, so the file holds n_dates x n_bands layers in the order
#   img1_B1, img1_B2, ... img1_timestamp, img2_B1, ...
# and terra names each layer <granuleID>_<band>.
# used endsWith() so "_B8" does not match "_B8A", nor "_B1" match "_B11"

# Indexes are held in two sidecars:
#   dates/<key>.csv  layer index -> UTC stamp
#   meta/<key>.csv   per-image solar and view geometry, plus rel_azimuth // used for local RTMs
# Sourced by scripts 3b, 3e and 4.

S2_STACK_ROOT <- "D:/PilbaraLaiProject/S2_annual_stacks"

# Every site-year stack on disk, with its sidecars
# data.frame(Site, Year, key, path, dates_csv, meta_csv)
stack_index <- function(root = S2_STACK_ROOT) {
  
  f <- Sys.glob(file.path(root, "stacks", "*.tif"))
  if (!length(f)) return(NULL)
  
  key <- tools::file_path_sans_ext(basename(f))
  
  # Keys are <SITE>_<YEAR>, with a trailing a/b chain if script 2a had to halve
  # the window to keep exports under the GEE band limit (YU_2019a, YU_2019ba)
  # Sites can contain digits - strip the part suffix, then the year
  
  base <- sub("[ab]+$", "", key)
  yr   <- as.integer(substr(base, nchar(base) - 3L, nchar(base)))
  st   <- substr(base, 1L, nchar(base) - 5L)
  
  d <- data.frame(Site = st, Year = yr, key = key, path = f,
                  dates_csv = file.path(root, "dates", paste0(key, ".csv")),
                  meta_csv  = file.path(root, "meta",  paste0(key, ".csv")),
                  stringsAsFactors = FALSE)
  
  d[order(d$Site, d$Year), ]
}

# Capture dates for one stack, in layer order
# @param tz timezone to report in; stamps are written UTC
stack_dates <- function(dates_csv, tz = "Etc/GMT-8") {
  if (!file.exists(dates_csv)) return(NULL)
  d <- read.csv(dates_csv, stringsAsFactors = FALSE)
  as.POSIXct(d$stamp, format = "%Y_%m_%d_%H_%M_%S", tz = "UTC") |>
    (\(x) as.POSIXct(format(x, tz = tz), tz = tz))()
}

# Per-image geometry for one stack
stack_meta <- function(meta_csv) {
  if (!file.exists(meta_csv)) return(NULL)
  read.csv(meta_csv, stringsAsFactors = FALSE)
}

# Single band timeseries (ordered by date)
# path= stack GeoTIFF
# band = band name as exported, e.g. "B4" or "NDVI"
# dates = optional POSIXct to name the layers with
# returns SpatRaster with one layer per date, in date order
stack_band <- function(path, band, dates = NULL) {
  r <- terra::rast(path)
  i <- which(endsWith(names(r), paste0("_", band)))
  if (!length(i))
    stop("band ", band, " not found in ", basename(path))
  out <- r[[i]]
  if (!is.null(dates)) {
    if (length(dates) != terra::nlyr(out))
      stop("date/layer mismatch in ", basename(path), ": ", terra::nlyr(out),
           " layers vs ", length(dates), " dates")
    names(out) <- format(dates, usetz = TRUE)
  }
  out
}



# Open a stack once keeping its layer names
# returns: list(r, nm)
# The cache lives inside the closure, so it travels with the function into a
# parallel worker instead of having to exist in that worker's global environment.
stack_open <- local({
  cache <- new.env(parent = emptyenv())
  function(path) {
    z <- get0(path, envir = cache, inherits = FALSE)

    if (!is.null(z)) {
      ok <- tryCatch({ terra::nlyr(z$r); TRUE }, error = function(e) FALSE)
      if (ok) return(z)
    }

    r <- terra::rast(path)
    z <- list(r = r, nm = names(r))

    # each worker only loads one site's stacks, but safer to clear if it grows anyway
    if (length(ls(cache)) > 24L) rm(list = ls(cache), envir = cache)
    assign(path, z, envir = cache)
    z
  }
})

# Several bands for one image, as a single multi-layer raster
#
# Used by the RTM inversion, which needs all the reflectance bands of one
# image together.
# i = image index within the stack
stack_image <- function(path, bands, i) {
  
  z <- stack_open(path)
  
  sel <- vapply(bands, function(b) {
    k <- which(endsWith(z$nm, paste0("_", b)))
    if (length(k) < i) stop("image ", i, " missing band ", b)
    k[i]
  }, integer(1))
  
  out <- z$r[[sel]]
  names(out) <- bands
  out
  
}

# All dates of one band for one site, across every year, in date order
# The site-year files are separate rasters on the same grid, so they can be
# concatenated. Returns NULL when the site has no stacks.
site_band_series <- function(site, band, idx = stack_index(), tz = "Etc/GMT-8") {
  
  if (is.null(idx)) return(NULL)
  
  s <- idx[idx$Site == site, ]
  if (!nrow(s)) return(NULL)
  
  lay <- list(); tm <- as.POSIXct(character(0), tz = tz)
  for (k in seq_len(nrow(s))) {
    
    d <- stack_dates(s$dates_csv[k], tz = tz)
    
    if (is.null(d)) next
    
    r <- try(stack_band(s$path[k], band), silent = TRUE)
    
    if (inherits(r, "try-error")) next
    
    if (terra::nlyr(r) != length(d)) {
      warning("date/layer mismatch in ", s$key[k]); next
    }
    
    lay <- c(lay, list(r))
    tm <- c(tm, d)
  }
  
  if (!length(lay)) return(NULL)
  
  out <- do.call(c, lay)
  
  o <- order(tm)
  out <- out[[o]]
  tm <- tm[o]
  
  names(out) <- format(tm, usetz = TRUE)
  out
}

# Geometry for one site across every year, in date order
site_meta <- function(site, idx = stack_index(), tz = "Etc/GMT-8") {
  if (is.null(idx)) return(NULL)
  s <- idx[idx$Site == site, ]
  if (!nrow(s)) return(NULL)
  m <- do.call(rbind, lapply(seq_len(nrow(s)), function(k) {
    a <- stack_meta(s$meta_csv[k]); if (is.null(a)) return(NULL)
    a$key <- s$key[k]
    a
  }))
  if (is.null(m)) return(NULL)
  m$DateTime <- as.POSIXct(format(as.POSIXct(m$stamp, format = "%Y_%m_%d_%H_%M_%S",
                                             tz = "UTC"), tz = tz), tz = tz)
  m$Date <- as.Date(m$DateTime)
  m[order(m$DateTime), ]
}

# Vegetation indices from reflectance, scale factor 10000.
# Formulas as in script 2b's add_indices(), which computed them in Earth Engine
# and exported them as extra bands; DEA ships reflectance only. Called by both
# scripts 3b and 4 so the raw and filtered arms cannot drift apart.
VI_NAMES <- c("NDVI","EVI","MSAVI2","NIRv","SR","DVI","EVI2","MSR","kNDVI",
              "GCVI","NDGI")

# site = site code, band = one of VI_NAMES, idx = stack_index()
# lyrs = optional layer subset applied to the reflectance before computing
site_vi_series <- function(site, band, idx = stack_index(), lyrs = NULL) {
  g <- function(b) {
    r <- site_band_series(site, b, idx)
    if (is.null(r)) return(NULL)
    if (!is.null(lyrs)) r <- r[[lyrs]]
    r / 10000
  }
  NIR <- g("B8"); RED <- g("B4"); BLUE <- g("B2"); GREEN <- g("B3")
  if (any(vapply(list(NIR, RED, BLUE, GREEN), is.null, logical(1)))) return(NULL)
  ndvi <- function() (NIR - RED) / (NIR + RED)
  sr   <- function() NIR / RED
  switch(band,
    NDVI   = ndvi(),
    EVI    = 2.5 * (NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1),
    MSAVI2 = (2 * NIR + 1 - sqrt((2 * NIR + 1)^2 - 8 * (NIR - RED))) / 2,
    NIRv   = ndvi() * NIR,
    SR     = sr(),
    DVI    = NIR - RED,
    EVI2   = 2.5 * (NIR - RED) / (NIR + 2.4 * RED + 1),
    MSR    = { s <- sr(); (s - 1) / sqrt(s + 1) },
    kNDVI  = tanh(ndvi()^2),
    GCVI   = NIR / GREEN - 1,
    NDGI   = (GREEN - RED) / (GREEN + RED),
    stop("unknown index: ", band))
}

# TRUE where any bit in `mask` is set in the integer raster `qc`.
# terra has no bitwAnd method for SpatRaster, so each bit is tested with
# integer division and modulo.
bits_set <- function(qc, mask) {
  out <- NULL
  for (b in 2L^(0:7)) {
    if (bitwAnd(as.integer(mask), as.integer(b)) == 0L) next
    t <- (qc %/% b) %% 2 >= 1
    out <- if (is.null(out)) t else (out | t)
  }
  out
}

# Apply a model's QC sibling (_qc-ts.tif) to its LAI stack. Bit 1 is
# reflectance outside the model's training range, bit 2 PROSAIL's RMSE flag.
# No sibling means nothing to apply.
# qc   = in-memory QC stack, used instead of reading the sibling
# lyrs = layer subset already applied to r, to match the QC stack to it
apply_qc_mask <- function(r, path = NULL, bits = 1L, lyrs = NULL, qc = NULL) {
  if (is.null(bits) || !length(bits) || bits <= 0) return(r)
  lab <- if (is.null(path)) "QC stack" else basename(path)
  if (is.null(qc)) {
    if (is.null(path)) return(r)
    qc_path <- sub("-ts.tif", "_qc-ts.tif", path, fixed = TRUE)
    if (!file.exists(qc_path)) return(r)
    lab <- basename(qc_path)
    qc <- try(terra::rast(qc_path), silent = TRUE)
    if (inherits(qc, "try-error")) return(r)
  }
  if (!is.null(lyrs)) {
    if (terra::nlyr(qc) < max(lyrs)) {
      warning("QC layer count does not match ", lab)
      return(r)
    }
    qc <- qc[[lyrs]]
  }
  if (terra::nlyr(qc) != terra::nlyr(r)) {
    warning("QC layer count does not match ", lab)
    return(r)
  }
  bad <- bits_set(qc, bits)
  if (is.null(bad)) return(r)
  terra::mask(r, bad, maskvalues = 1)
}

# DEA QC bit flags for one site across every year, in date order.
# Script 2a writes qc/<SITE>_<year>.tif beside stacks/, one uint8 band per date:
#   1 cloud, DEA s2cloudless-mask    2 cloud, s2cloudless-prob threshold
#   4 fmask cloud or shadow          8 non-contiguous
#  16 terrain shadow                 32 no data
# NULL means the site has no QC stacks, i.e. a download predating them.
site_qc_series <- function(site, idx = stack_index(), tz = "Etc/GMT-8") {
  if (is.null(idx)) return(NULL)
  s <- idx[idx$Site == site, ]
  if (!nrow(s)) return(NULL)
  lay <- list(); tm <- as.POSIXct(character(0), tz = tz)
  for (k in seq_len(nrow(s))) {
    f <- file.path(dirname(dirname(s$path[k])), "qc", basename(s$path[k]))
    if (!file.exists(f)) next
    d <- stack_dates(s$dates_csv[k], tz = tz)
    if (is.null(d)) next
    r <- try(terra::rast(f), silent = TRUE)
    if (inherits(r, "try-error")) next
    if (terra::nlyr(r) != length(d)) {
      warning("date/layer mismatch in QC stack ", basename(f)); next
    }
    lay <- c(lay, list(r)); tm <- c(tm, d)
  }
  if (!length(lay)) return(NULL)
  out <- do.call(c, lay)
  o <- order(tm); out <- out[[o]]; tm <- tm[o]
  names(out) <- format(tm, usetz = TRUE)
  out
}
