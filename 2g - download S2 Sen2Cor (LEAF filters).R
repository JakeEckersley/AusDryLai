################################################################################
#### 2g - Sentinel-2 Sen2Cor reflectance with the LEAF input filters          ####
################################################################################
# One annual stack per site-year of Sen2Cor L2A reflectance (COPERNICUS/S2_SR_HARMONIZED),
# by Earth Engine batch export to Cloud Storage. The images are screened as the LEAF
# toolbox screens its Sentinel-2 inputs before the networks (script 2f), with the same
# rules written directly in Earth Engine:
#   scene cloud   CLOUDY_PIXEL_PERCENTAGE below max_cloud
#   s2cloudless   the scene has a COPERNICUS/S2_CLOUD_PROBABILITY image (LEAF joins to it, so a
#                 scene without one drops out); the probability itself is not used
#   QA60          opaque cloud (bit 10) and cirrus (bit 11) clear
#   SCL           scene classes 4 (vegetation) and 5 (not vegetated) only
#   date gate     an image is kept when min_valid_frac of the site survives the masks
# Masked pixels are NA in the stacks and every tile delivery is kept. Script 4 reads the
# stacks as its raw arm: it computes the vegetation indices and takes one image per date.
#
# layout, read by R/stack_io.R:
#   <out_dir>/stacks/<SITE>_<year>.tif   n_dates x n_bands layers, image-major
#   <out_dir>/dates/<SITE>_<year>.csv    acquisition times, one row per image
#   <out_dir>/meta/<SITE>_<year>.csv     per-image solar and view geometry
#   <out_dir>/done/<SITE>_<year>.done    completion marker
#
# .done markers make a restart safe: finished site-years are skipped.

library(terra)
library(reticulate)
library(jsonlite)

## ---- config -----------------------------------------------------------------
py_ver      <- "C:/Users/eck025/.conda/envs/s2/python.exe"
ee_project  <- "sent-dryland-validation"
serv_email  <- "jakeeckersley@sent-dryland-validation.iam.gserviceaccount.com"
serv_key    <- "D:/PilbaraLaiProject/serv-key/sent-dryland-validation-e2e627d1d6e2.json"
bucket_name <- "s2-bucket"

out_dir     <- "D:/PilbaraLaiProject/Sen2Cor_LEAFfilter_stacks"
sites_path  <- "data/SiteBoundsForSatDownloads"
sites       <- NULL                     # NULL = every site with a bounds shapefile, or a character vector
years       <- 2015:2024                # clipped to start_date and end_date
start_date  <- "2015-03-01"
end_date    <- "2024-09-01"
out_scale   <- 10L                      # native Sentinel-2 grid
out_crs     <- "EPSG:32750"
max_cloud   <- 50                       # as LEAF
min_valid_frac <- 0.80                  # as script 2f
bands       <- c("B2", "B3", "B4", "B8")   # the vegetation indices and the 10 m SL2P network
angle_band  <- "B8A"                    # SL2P keys its view geometry off B8A
max_running <- 25L                      # tasks in flight at once
max_imgs    <- 4500L %/% length(bands)  # images per export before the window is halved:
                                        # toBands() makes n_images x n_bands layers and
                                        # GEE refuses a request past ~5000 bands

poll_secs   <- 30L
overwrite   <- FALSE

for (s in c("stacks", "dates", "meta", "done"))
  dir.create(file.path(out_dir, s), showWarnings = FALSE, recursive = TRUE)

## ---- python side ------------------------------------------------------------
use_python(py_ver, required = TRUE)
ee <- import("ee")
stopifnot(file.exists(serv_key))
ee$Initialize(credentials = ee$ServiceAccountCredentials(serv_email, serv_key),
              project = ee_project)

py_run_string(paste0("
import ee
from google.cloud import storage
from google.oauth2 import service_account as sa

MIN_VALID_FRAC = ", min_valid_frac, "
_key = r'", serv_key, "'
_gcs = storage.Client(credentials=sa.Credentials.from_service_account_file(_key),
                      project='", ee_project, "')

def mask_clear(img):
    # QA60 opaque cloud (bit 10) and cirrus (bit 11) clear
    qa = img.select('QA60')
    return img.updateMask(qa.bitwiseAnd(1 << 10).eq(0).And(qa.bitwiseAnd(1 << 11).eq(0)))

def mask_land(img):
    # Sen2Cor scene classification 4 (vegetation) or 5 (not vegetated)
    scl = img.select('SCL')
    return img.updateMask(scl.eq(4).Or(scl.eq(5)))

def add_valid_frac(img, aoi, scale):
    # fraction of the AOI still unmasked after the masks above. (specific band doesn't matter)
    v = img.select('B4').mask().rename('v')
    f = v.reduceRegion(reducer=ee.Reducer.mean(), geometry=aoi, scale=scale,
                       maxPixels=int(1e9), bestEffort=True).get('v')
    return img.set('valid_frac', f)

def build(aoi, d0, d1, maxcloud, scale):
    s2 = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
            .filterBounds(aoi).filterDate(d0, d1)
            .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', maxcloud))
    cl = (ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
            .filterBounds(aoi).filterDate(d0, d1))
    s2 = ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply(
        primary=s2, secondary=cl,
        condition=ee.Filter.equals(leftField='system:index', rightField='system:index')))
    return (s2.map(mask_clear).map(mask_land)
              .map(lambda im: add_valid_frac(im, aoi, scale))
              .filter(ee.Filter.gte('valid_frac', MIN_VALID_FRAC))
              .sort('system:time_start'))

def n_images(ic):
    return ic.size().getInfo()

def submit(ic, bands, aoi, scale, crs, bucket, prefix):
    t = ee.batch.Export.image.toCloudStorage(**{
        'image': ic.select(bands).toBands().toFloat(),
        'description': prefix[:100], 'bucket': bucket, 'fileNamePrefix': prefix,
        'scale': scale, 'region': aoi, 'crs': crs,
        'fileFormat': 'GeoTIFF', 'maxPixels': int(1e13)})
    t.start()
    return t

def state(t):
    return t.status().get('state')

def meta(ic, props):
    out = {'stamp': ic.aggregate_array('system:time_start').map(
              lambda x: ee.Date(x).format('YYYY_MM_dd_HH_mm_ss')).getInfo()}
    for p in props:
        out[p] = ic.aggregate_array(p).getInfo()
    return out

def blob_exists(bucket, name):
    return _gcs.bucket(bucket).blob(name).exists()

def blob_download(bucket, name, dest):
    _gcs.bucket(bucket).blob(name).download_to_filename(dest)

def blob_delete(bucket, name):
    b = _gcs.bucket(bucket).blob(name)
    if b.exists(): b.delete()
"))

## ---------------------------------------------------------------------
## ---- RUN ------------------------------------------------------------
## ---------------------------------------------------------------------

site_list <- sort(tools::file_path_sans_ext(
  basename(Sys.glob(file.path(sites_path, "*.shp")))))
if (!is.null(sites)) site_list <- intersect(site_list, sites)

props <- c("SPACECRAFT_NAME", "MEAN_SOLAR_ZENITH_ANGLE", "MEAN_SOLAR_AZIMUTH_ANGLE",
           paste0("MEAN_INCIDENCE_ZENITH_ANGLE_",  angle_band),
           paste0("MEAN_INCIDENCE_AZIMUTH_ANGLE_", angle_band),
           "CLOUDY_PIXEL_PERCENTAGE", "PRODUCT_ID")

aoi_of <- function(site) {
  tmp <- tempfile(fileext = ".geojson"); on.exit(unlink(tmp), add = TRUE)
  writeVector(project(vect(file.path(sites_path, paste0(site, ".shp"))), "EPSG:4326"),
              tmp, filetype = "GeoJSON", overwrite = TRUE)
  ee$FeatureCollection(fromJSON(tmp, simplifyVector = FALSE))$geometry()$bounds()
}

jobs <- expand.grid(site = site_list, year = years, stringsAsFactors = FALSE)
jobs$key  <- sprintf("%s_%d", jobs$site, jobs$year)
jobs$done <- file.path(out_dir, "done", paste0(jobs$key, ".done"))
jobs <- jobs[!file.exists(jobs$done) | overwrite, ]
message(nrow(jobs), " site-years to export")

## ---- submit, poll, download -------------------------------------------------
# Tasks are kept to max_running in flight: submits, waits for space, then downloads the
# finished ones. GEE does the work concurrently.
inflight <- list()

harvest <- function(force = FALSE) {
  keep <- list()
  for (h in inflight) {
    st <- tryCatch(py_to_r(py$state(h$task)), error = function(e) "UNKNOWN")
    if (st %in% c("COMPLETED", "FAILED", "CANCELLED")) {
      if (st == "COMPLETED") {
        blob <- paste0(h$key, ".tif")
        if (py_to_r(py$blob_exists(bucket_name, blob))) {
          dest <- file.path(out_dir, "stacks", blob)
          py$blob_download(bucket_name, blob, dest)
          py$blob_delete(bucket_name, blob)          # keep the bucket clean
          write.csv(h$dates, file.path(out_dir, "dates", paste0(h$key, ".csv")),
                    row.names = FALSE)
          write.csv(h$meta,  file.path(out_dir, "meta",  paste0(h$key, ".csv")),
                    row.names = FALSE)
          # a split year is only finished when every part has landed
          left <- get0(h$done_path, envir = parts_left, ifnotfound = 1L) - 1L
          assign(h$done_path, left, envir = parts_left)
          if (left <= 0L) file.create(h$done_path)
          message("  downloaded ", blob)
        } else message("  ", h$key, ": completed but no blob")
      } else message("  ", h$key, ": ", st)
    } else keep[[length(keep) + 1L]] <- h
  }
  inflight <<- keep
  invisible(NULL)
}

# One site-year is usually one export, but a year with too many images is split or GEE
# rejects it. Parts are named <key>a, <key>b ... and are stitched back together by date
# downstream. The site-year .done marker is only written once every part has landed.
parts_left <- new.env(parent = emptyenv())

submit_window <- function(key, site, d0, d1, done_path, depth = 0L) {
  aoi <- aoi_of(site)
  ic  <- ee$ImageCollection(py$build(aoi, format(d0), format(d1), max_cloud, out_scale))
  n   <- tryCatch(py_to_r(py$n_images(ic)), error = function(e) NA_integer_)
  if (is.na(n)) { message('  ', key, ': collection failed'); return(0L) }
  if (!n) return(0L)

  if (n > max_imgs && depth < 3L) {
    mid <- d0 + floor(as.numeric(d1 - d0) / 2)
    if (mid > d0 && mid < d1) {
      message(sprintf('  %s: %d images, splitting', key, n))
      a <- submit_window(paste0(key, 'a'), site, d0, mid, done_path, depth + 1L)
      b <- submit_window(paste0(key, 'b'), site, mid, d1, done_path, depth + 1L)
      return(a + b)
    }
  }

  md <- try(py_to_r(py$meta(ic, as.list(props))), silent = TRUE)
  if (inherits(md, 'try-error')) { message('  ', key, ': metadata failed'); return(0L) }
  mdf <- as.data.frame(lapply(md, unlist), stringsAsFactors = FALSE)
  saa <- mdf[['MEAN_SOLAR_AZIMUTH_ANGLE']]
  vaa <- mdf[[paste0('MEAN_INCIDENCE_AZIMUTH_ANGLE_', angle_band)]]
  mdf$rel_azimuth <- abs((saa - vaa + 180) %% 360 - 180)
  mdf$Site <- site; mdf$layer <- seq_len(nrow(mdf))

  task <- py$submit(ic, as.list(bands), aoi, out_scale, out_crs,
                    bucket_name, key)
  inflight[[length(inflight) + 1L]] <<-
    list(task = task, key = key, done_path = done_path,
         dates = mdf[, c('layer', 'stamp')], meta = mdf)

  message(sprintf('submitted %-14s %3d images x %d bands', key, n, length(bands)))
  1L
}

for (i in seq_len(nrow(jobs))) {
  while (length(inflight) >= max_running) { Sys.sleep(poll_secs); harvest() }
  j   <- jobs[i, ]
  d0  <- max(as.Date(start_date), as.Date(sprintf('%d-01-01', j$year)))
  d1  <- min(as.Date(end_date),   as.Date(sprintf('%d-01-01', j$year + 1)))
  if (d1 <= d0) { file.create(j$done); next }
  n_parts <- submit_window(j$key, j$site, d0, d1, j$done)
  if (n_parts == 0L) { file.create(j$done); next }   # nothing to fetch
  assign(j$done, n_parts, envir = parts_left)
}

while (length(inflight)) { Sys.sleep(poll_secs); harvest() }
message("written to ", normalizePath(out_dir))
