################################################################################
#### 3b - Sentinel-2 NBART annual stacks from Digital Earth Australia       ####
################################################################################
# source: https://knowledge.dea.ga.gov.au/data/product/dea-surface-reflectance-nbart-sentinel-2c-msi/?tab=access
# Same output layout as script 2b, so scripts 3b and 4 can read either:
#   <out_dir>/stacks/<SITE>_<year>.tif   annual stack, n_dates x n_bands
#   <out_dir>/dates/<SITE>_<year>.csv    acquisition times, one row per date
#   <out_dir>/meta/<SITE>_<year>.csv     per-granule geometry and cloud cover
#   <out_dir>/done/<SITE>_<year>.done    completion marker, for safe restarts
#
# What this gives that 3a does not:
#   NBART  = BRDF-adjusted to nadir AND terrain-illumination corrected. Script
#            3a downloads COPERNICUS/S2_SR_HARMONIZED, which is Sen2Cor surface
#            reflectance with neither correction applied.
#   angles = per-pixel solar/view rasters rather than one scene mean, though
#            over a 500 m footprint they barely vary (0.01 deg across YU)
#   masks  = fmask, s2cloudless and a combined terrain shadow band, supplied
#            as data rather than computed in Earth Engine
#
# Three collections, one per instrument, searched together:
#   ga_s2am_ard_3  Sentinel-2A
#   ga_s2bm_ard_3  Sentinel-2B
#   ga_s2cm_ard_3  Sentinel-2C
#
# Public S3 COGs, read through /vsicurl with no account, key or quota. Only the
# site footprint is fetched, so a date costs about 5 s rather than a scene.

library(terra)
library(jsonlite)
library(future.apply)
library(parallelly)
source("R/stack_io.R")   # bits_set

## ---- config -----------------------------------------------------------------
out_dir    <- "D:/PilbaraLaiProject/DEA_NBART_stacks"
sites_path <- "data/SiteBoundsForSatDownloads"
start_date <- "2015-03-01"
end_date   <- "2024-09-01"
max_cloud  <- 50          # scene cloud cover %, as in 3a
out_crs    <- "EPSG:32750"
out_scale  <- 10          # all bands resampled to the 10 m grid, as GEE did
n_workers  <- 6

# band01 and band09 are 60 m atmospheric bands that nothing downstream uses.
# DEA names the files bandNN; layers are written with the ESA names so that
# R/stack_io.R, which matches on "_B4", reads these stacks unchanged.
bands <- c(band02 = "B2",  band03 = "B3",  band04 = "B4", band05 = "B5",
           band06 = "B6",  band07 = "B7",  band08 = "B8", band08a = "B8A",
           band11 = "B11", band12 = "B12")
oa_angle_bands <- c("solar-zenith","satellite-view","relative-azimuth")

# Categorical masks, packed into a per-date bit flag in qc/ and left unapplied;
# script 3b decides which to use.
oa_mask_bands <- c("s2cloudless-mask","s2cloudless-prob","fmask",
                   "nbart-contiguity","combined-terrain-shadow")

# Cloud probability cutoff for the s2cloudless-prob band (0-1). DEA's shipped
# s2cloudless-mask uses an undocumented cutoff, so both are written as separate
# bits. 0.10 is script 2b's value, chosen visually to catch Pilbara smoke.
# NA uses DEA's mask only. DEA reports s2cloudless over-classifies cloud over
# bare ground from 2022-01-25 to 2024-06-07 and recommends fmask there.
cloud_prob_thresh <- 0.10

# Bit flags in qc/<SITE>_<year>.tif, one uint8 band per date
QC_CLOUD_DEA  <- 1L    # s2cloudless-mask == 2
QC_CLOUD_PROB <- 2L    # s2cloudless-prob > cloud_prob_thresh
QC_FMASK      <- 4L    # fmask cloud (2) or cloud shadow (3)
QC_NOCONTIG   <- 8L    # nbart-contiguity == 0
QC_TERRAIN    <- 16L   # combined-terrain-shadow == 0
QC_NODATA     <- 32L   # s2cloudless-mask == 0 or fmask == 0

# Bits counted against the per-date frac_valid in the meta sidecar. One cloud
# bit only, or cloud is charged twice where the two masks agree. Terrain shadow
# is left out: it measures 1.000 lit across the whole record.
qc_valid_bits <- QC_CLOUD_PROB + QC_FMASK + QC_NOCONTIG + QC_NODATA
collections <- c("ga_s2am_ard_3","ga_s2bm_ard_3","ga_s2cm_ard_3")
stac <- "https://explorer.dea.ga.gov.au/stac/search"

Sys.setenv(GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
           CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".tif",
           GDAL_HTTP_MULTIPLEX = "YES", VSI_CACHE = "TRUE")
options(timeout = 300)   # default 60 s is too short for a busy STAC page
for (d in c("stacks","dates","meta","done","qc"))
  dir.create(file.path(out_dir, d), showWarnings = FALSE, recursive = TRUE)

## ---- sites ------------------------------------------------------------------
# one shapefile per site, named <SITE>.shp, as script 2b expects
site_list <- sort(tools::file_path_sans_ext(
  basename(Sys.glob(file.path(sites_path, "*.shp")))))
site_aoi <- function(site) terra::aggregate(
  vect(file.path(sites_path, paste0(site, ".shp"))))
years <- seq(as.integer(format(as.Date(start_date), "%Y")),
             as.integer(format(as.Date(end_date),   "%Y")))

## ---- STAC -------------------------------------------------------------------
# DEA paginates by offset through a rel="next" link (?_o=100), not a page
# parameter - passing page= is silently ignored and every request returns the
# same first features, so follow the link the server gives.
stac_items <- function(bb, d0, d1, max_pages = 50, tries = 4) {
  q <- sprintf("%s?collections=%s&bbox=%s&datetime=%s/%s&limit=100",
               stac, paste(collections, collapse = ","),
               paste(bb, collapse = ","), d0, d1)
  out <- list(); page <- 0
  # A page that times out must abort the whole site-year, not end the loop
  # quietly: silent truncation here cost 157 of 252 dates on a test run and
  # still wrote a .done marker, so the gap would never have been noticed.
  while (!is.null(q) && page < max_pages) {
    j <- NULL
    for (a in seq_len(tries)) {
      j <- try(fromJSON(q, simplifyVector = FALSE), silent = TRUE)
      if (!inherits(j, "try-error")) break
      Sys.sleep(2 * a)
    }
    if (inherits(j, "try-error")) {
      warning("STAC page failed after ", tries, " tries: ", q, call. = FALSE)
      return(structure(list(), incomplete = TRUE))
    }
    if (!length(j$features)) break
    out <- c(out, j$features)
    nx <- Filter(function(l) identical(l$rel, "next"), j$links)
    q <- if (length(nx)) nx[[1]]$href else NULL
    page <- page + 1
  }
  if (page >= max_pages) return(structure(list(), incomplete = TRUE))
  if (!length(out)) return(NULL)
  out <- out[!duplicated(vapply(out, function(f) f$id, character(1)))]
  do.call(rbind, lapply(out, function(f) {
    a <- f$assets[["nbart_red"]]$href
    if (is.null(a)) return(NULL)
    data.frame(
      id       = f$id,
      datetime = f$properties$datetime,
      cloud    = f$properties[["eo:cloud_cover"]] %||% NA_real_,
      platform = f$properties$platform %||% NA_character_,
      tile     = f$properties[["odc:region_code"]] %||% NA_character_,
      # strip the band suffix to leave the per-scene prefix
      prefix   = sub("_band04\\.tif$", "", sub("^s3://dea-public-data/",
                     "/vsicurl/https://data.dea.ga.gov.au/", a)),
      stringsAsFactors = FALSE)
  }))
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- one scene over one site -------------------------------------------------
# Reads only the site window. 20 m and 60 m bands are resampled onto the 10 m
# grid of band04, which is what GEE did implicitly in script 2b.
# g is a fixed per-site template. Several sites straddle two MGRS tiles, and a
# scene from one tile covers only part of the AOI - cropping then returns a
# clipped extent, so dates from different tiles disagree and rast() refuses to
# stack them. Resampling every scene onto one template makes partial coverage
# come back as NA instead, at a consistent size.
read_scene <- function(prefix, aoi, g, keep_bands = names(bands)) {
  nb <- sub("_oa_", "_nbart_", prefix)
  r0 <- try(rast(paste0(nb, "_band04.tif")), silent = TRUE)
  if (inherits(r0, "try-error")) return(NULL)
  if (is.null(intersect(ext(r0), ext(g)))) return(NULL)   # no overlap at all
  # categorical bands need nearest neighbour; bilinear invents class codes
  grab <- function(f, method = "bilinear") {
    r <- try(rast(f), silent = TRUE); if (inherits(r, "try-error")) return(NULL)
    if (is.null(intersect(ext(r), ext(g)))) return(NULL)
    x <- try(resample(crop(r, ext(g), snap = "out"), g, method = method),
             silent = TRUE)
    if (inherits(x, "try-error") || !ncell(x)) return(NULL)
    x
  }
  refl <- lapply(keep_bands, function(b) grab(paste0(nb, "_", b, ".tif")))
  if (any(vapply(refl, is.null, logical(1)))) return(NULL)
  refl <- rast(refl); names(refl) <- keep_bands

  oa  <- sub("_nbart_", "_oa_", nb)
  ang <- lapply(oa_angle_bands, function(b) grab(paste0(oa, "_", b, ".tif")))
  m   <- lapply(oa_mask_bands, function(b)
           grab(paste0(oa, "_", b, ".tif"),
                method = if (b == "s2cloudless-prob") "bilinear" else "near"))
  names(m) <- oa_mask_bands

  # Pack the masks into one uint8 bit flag. NA contributes no bit, so a band DEA
  # did not ship cannot silently mask the whole scene.
  z    <- refl[[1]] * 0
  flag <- function(x, bit) if (is.null(x)) z else bit * ifel(is.na(x), 0, x)
  scl  <- m[["s2cloudless-mask"]]; prb <- m[["s2cloudless-prob"]]
  fm   <- m[["fmask"]]
  cld  <- if (is.null(prb) || is.na(cloud_prob_thresh)) NULL
          else prb > cloud_prob_thresh
  qc <- flag(if (is.null(scl)) NULL else scl == 2, QC_CLOUD_DEA) +
        flag(cld, QC_CLOUD_PROB) +
        flag(if (is.null(fm)) NULL else fm == 2 | fm == 3, QC_FMASK) +
        flag(if (is.null(m[["nbart-contiguity"]])) NULL else
               m[["nbart-contiguity"]] == 0, QC_NOCONTIG) +
        flag(if (is.null(m[["combined-terrain-shadow"]])) NULL else
               m[["combined-terrain-shadow"]] == 0, QC_TERRAIN) +
        flag(if (is.null(scl) && is.null(fm)) NULL
             else if (is.null(scl)) fm == 0
             else if (is.null(fm)) scl == 0 else scl == 0 | fm == 0, QC_NODATA)

  fr <- function(x, v) if (is.null(x)) NA else mean(values(x) == v, na.rm = TRUE)
  # AOI fraction carrying data and unflagged by qc_valid_bits. The scene-level
  # max_cloud prefilter cannot give this: a scene under 50% cloud tile-wide can
  # still be wholly clouded over a 100 m site. Recorded, not used to drop dates.
  valid <- ifel(is.na(refl[[1]]), 0, 1) * ifel(bits_set(qc, qc_valid_bits), 0, 1)
  list(refl = refl, qc = qc,
       frac_valid = mean(values(valid), na.rm = TRUE),
       stats = c(sza = if (is.null(ang[[1]])) NA else mean(values(ang[[1]]), na.rm = TRUE),
                 vza = if (is.null(ang[[2]])) NA else mean(values(ang[[2]]), na.rm = TRUE),
                 raa = if (is.null(ang[[3]])) NA else mean(values(ang[[3]]), na.rm = TRUE),
                 frac_lit   = fr(m[["combined-terrain-shadow"]], 1),
                 frac_clear = fr(fm, 1),
                 frac_cloud_dea  = fr(scl, 2),
                 frac_cloud_prob = if (is.null(cld)) NA else
                                     mean(values(cld), na.rm = TRUE),
                 frac_contig = fr(m[["nbart-contiguity"]], 1)))
}

## ---- one site-year -----------------------------------------------------------
do_site_year <- function(site, year, aoi) {
  done <- file.path(out_dir, "done", sprintf("%s_%d.done", site, year))
  if (file.exists(done)) return(sprintf("%s %d skip", site, year))
  d0 <- max(as.Date(start_date), as.Date(sprintf("%d-01-01", year)))
  d1 <- min(as.Date(end_date),   as.Date(sprintf("%d-12-31", year)))
  if (d1 <= d0) { file.create(done); return(sprintf("%s %d empty", site, year)) }

  bb <- as.vector(ext(project(aoi, "EPSG:4326")))[c(1, 3, 2, 4)]
  it <- stac_items(round(bb, 5), format(d0), format(d1))
  # no .done marker on an incomplete search, so the site-year is retried
  if (isTRUE(attr(it, "incomplete")))
    return(sprintf("%s %d STAC INCOMPLETE - not marked done, rerun", site, year))
  if (is.null(it)) { file.create(done); return(sprintf("%s %d no scenes", site, year)) }
  it <- it[is.na(it$cloud) | it$cloud < max_cloud, ]
  it <- it[order(it$datetime), ]
  if (!nrow(it)) { file.create(done); return(sprintf("%s %d all cloudy", site, year)) }

  # Per-site template. Reuse the grid of a year already downloaded so new years
  # align with it; site_band_series() stacks a site's years together and would
  # fail otherwise. Only derive a fresh grid when the site has nothing yet.
  prev <- Sys.glob(file.path(out_dir, "stacks", paste0(site, "_*.tif")))
  tmpl <- if (length(prev)) {
    rast(prev[1])[[1]]
  } else {
    a <- ext(project(aoi, out_crs))
    rast(ext(floor(a[1] / out_scale) * out_scale, ceiling(a[2] / out_scale) * out_scale,
             floor(a[3] / out_scale) * out_scale, ceiling(a[4] / out_scale) * out_scale),
         resolution = out_scale, crs = out_crs)
  }

  got <- lapply(seq_len(nrow(it)), function(i) read_scene(it$prefix[i], aoi, tmpl))
  ok  <- !vapply(got, is.null, logical(1))
  if (!any(ok)) { file.create(done); return(sprintf("%s %d no reads", site, year)) }
  it <- it[ok, ]; got <- got[ok]

  st <- rast(lapply(got, `[[`, "refl"))
  # ESA band names, so stack_band()'s endsWith("_B4") match works
  names(st) <- paste0(rep(it$id, each = length(bands)), "_", bands)
  writeRaster(st, file.path(out_dir, "stacks", sprintf("%s_%d.tif", site, year)),
              overwrite = TRUE, datatype = "INT2S")

  # QC bit flags, one band per date, in the same layer order as the stack.
  # Script 3b applies whichever bits it is configured for.
  qc <- rast(lapply(got, `[[`, "qc"))
  names(qc) <- it$id
  writeRaster(qc, file.path(out_dir, "qc", sprintf("%s_%d.tif", site, year)),
              overwrite = TRUE, datatype = "INT1U")

  # sidecars in script 2b's shape: stamp as %Y_%m_%d_%H_%M_%S, and the same
  # metadata column names, so scripts 3b, 3e and 4 need no changes
  stamp <- format(as.POSIXct(sub("Z$", "", it$datetime), format = "%Y-%m-%dT%H:%M:%OS",
                             tz = "UTC"), "%Y_%m_%d_%H_%M_%S", tz = "UTC")
  write.csv(data.frame(layer = seq_along(stamp), stamp = stamp),
            file.path(out_dir, "dates", sprintf("%s_%d.csv", site, year)),
            row.names = FALSE)

  s  <- as.data.frame(do.call(rbind, lapply(got, `[[`, "stats")))
  fv <- vapply(got, function(g) g$frac_valid %||% NA_real_, numeric(1))
  write.csv(data.frame(
    stamp = stamp,
    SPACECRAFT_NAME = sub("^sentinel-2([abc])$", "Sentinel-2\\U\\1", it$platform,
                          perl = TRUE),
    MEAN_SOLAR_ZENITH_ANGLE          = s$sza,
    MEAN_INCIDENCE_ZENITH_ANGLE_B8A  = s$vza,
    CLOUDY_PIXEL_PERCENTAGE          = it$cloud,
    PRODUCT_ID   = it$id,
    rel_azimuth  = s$raa,
    frac_lit     = s$frac_lit,      # 1 = unshadowed by terrain
    frac_clear   = s$frac_clear,    # fmask == 1
    frac_cloud_dea  = s$frac_cloud_dea,   # DEA's own s2cloudless-mask == 2
    frac_cloud_prob = s$frac_cloud_prob,  # prob > cloud_prob_thresh
    frac_contig  = s$frac_contig,   # 1 = all bands present at that pixel
    frac_valid   = fv,              # AOI fraction with data and no qc_valid_bits
    tile         = it$tile,
    Site         = site,
    layer        = seq_along(stamp),
    stringsAsFactors = FALSE),
    file.path(out_dir, "meta", sprintf("%s_%d.csv", site, year)), row.names = FALSE)
  file.create(done)
  sprintf("%s %d: %d dates", site, year, nrow(it))
}

## ---- run ---------------------------------------------------------------------
jobs <- expand.grid(site = site_list, year = years, stringsAsFactors = FALSE)
message(nrow(jobs), " site-years on ", n_workers, " workers")
plan(multisession, workers = n_workers)
# One bad scene must not cancel the other 239 site-years: future_lapply aborts
# every iteration on an unhandled error, so each job catches its own and no
# .done marker is written, leaving it to be retried on the next run.
log <- future_lapply(seq_len(nrow(jobs)), function(k)
  tryCatch(do_site_year(jobs$site[k], jobs$year[k], site_aoi(jobs$site[k])),
           error = function(e) sprintf("%s %d FAILED: %s", jobs$site[k],
                                       jobs$year[k], conditionMessage(e))),
  future.seed = TRUE, future.packages = c("terra", "jsonlite"))
plan(sequential)
writeLines(unlist(log))
