################################################################################
#### 3 - Sentinel-2 surface reflectance, annual stacks via GCS batch export   ####
################################################################################
# R equivalent of the older "3 - DownloadSentinel2viaGCS.ipynb", but one stack per
# site-year instead of one file per date, and using Earth Engine batch export to
# Cloud Storage rather than getDownloadURL. Significantly faster.
#
# layout matching the LEAF downloads:
#   <out_dir>/stacks/<SITE>_<year>.tif annual stack
#   <out_dir>/dates/<SITE>_<year>.csv  acquisition times, one row per date
#   <out_dir>/meta/<SITE>_<year>.csv   per-granule metadata (see below)
#   <out_dir>/done/<SITE>_<year>.done  completion marker
#
# metadata info:
#   SPACECRAFT_NAME  - used to select the appropriate spectral response curve from the workbook
#   solar/view angles  - inputs to the network inversion (script 3e)
#
# GEE and GCS downloads are finicky .done markers are used for safe restarts
# tasks already in the bucket are downloaded rather than re-exported.

library(terra)
library(reticulate)
library(jsonlite)

## ---- config -----------------------------------------------------------------
py_ver      <- "C:/Users/eck025/.conda/envs/s2/python.exe"
ee_project  <- "sent-dryland-validation"
serv_email  <- "jakeeckersley@sent-dryland-validation.iam.gserviceaccount.com"
serv_key    <- "D:/PilbaraLaiProject/serv-key/sent-dryland-validation-e2e627d1d6e2.json"
bucket_name <- "s2-bucket"

out_dir     <- "D:/PilbaraLaiProject/S2_annual_stacks"
sites_path  <- "data/SiteBoundsForSatDownloads"
start_date  <- "2015-03-01"
end_date    <- "2024-09-01"
out_scale   <- 10L                      # native Sentinel-2 grid
out_crs     <- "EPSG:32750"
max_cloud   <- 50                       # CLOUD_FILTER
cld_prb        <- 10   # CLD_PRB_THRESH s2cloudless probability above which a pixel
                       #                is called cloud, lowering it masks more,
                       #                10 is conservative - default was 15-20  
                       #                but visually confirmed in gee that it's needed to catch smoke
nir_drk_thresh <- 0.11 # NIR_DRK_THRESH lower = less shadow sensitivity, but less
                       #                likely to mask genuinely dark bare ground
cld_prj_dist   <- 2    # CLD_PRJ_DIST   distance from cloud to search for cloud shadow, km
buffer_m       <- 50   # BUFFER         dilation around the cloud mask, m
min_valid_frac <- 0.80 # MIN_VALID_FRAC keep an image only if at least this much
                       #                of the site AOI survives the whole mask
                       #                chain. One number covering cloud, shadow,
                       #                SCL, saturation and low reflectance -
                       #                partial scenes were a major source of
                       #                noise in site-mean timeseries.
aot_max        <- 0.30 # AOT_MAX        mean aerosol optical thickness over the AOI
                       #                above which the scene is dropped. Catches
                       #                smoke and dust haze that's ignored by the cloud mask (which looks for bright cold patches)
                       #                The L2A band is scaled by 0.001.

angle_band  <- "B8A"                    # SL2P keys its view geometry off B8A
max_running <- 25L                      # tasks in flight at once
max_imgs    <- 180L                     # images per export before the window is
                                        # halved. toBands() flattens the collection
                                        # to n_images x n_bands, and past
                                        # ~ 5000 bands GEE refuses the request (double orbits)
one_pass_a_day <- FALSE  # keeps only the nearer-nadir granule when a site is imaged
                                        # twice in a day. doesn't make much of a difference.

poll_secs   <- 30L
overwrite   <- FALSE

refl_bands <- c("B1","B2","B3","B4","B5","B6","B7","B8","B8A","B9","B11","B12")
vi_bands   <- c("NDVI","EVI","MSAVI2","NIRv","SR","DVI","EVI2","MSR","kNDVI", "GCVI","NDGI")
keep_bands <- c(refl_bands, vi_bands, "timestamp")

for (s in c("stacks", "dates", "meta", "done"))
  dir.create(file.path(out_dir, s), showWarnings = FALSE, recursive = TRUE)

## ---- python side ------------------------------------------------------------
use_python(py_ver, required = TRUE)
ee <- import("ee")
stopifnot(file.exists(serv_key))
ee$Initialize(credentials = ee$ServiceAccountCredentials(serv_email, serv_key),
              project = ee_project)

# the full masking chain, ported verbatim from the notebook. The earlier version
# of this script had only the s2cloudless cloud test - no shadow projection, no
# SCL mask, no saturation or low-reflectance masks.
py_run_string(paste0("
import ee
from google.cloud import storage
from google.oauth2 import service_account as sa

CLD_PRB_THRESH = ", cld_prb, "
NIR_DRK_THRESH = ", nir_drk_thresh, "
CLD_PRJ_DIST   = ", cld_prj_dist, "
BUFFER         = ", buffer_m, "
MIN_VALID_FRAC = ", min_valid_frac, "
AOT_MAX        = ", aot_max, "
ONE_PASS_A_DAY = ", ifelse(one_pass_a_day, 'True', 'False'), "
SR_BAND_SCALE  = 1e4
_key = r'", serv_key, "'
_gcs = storage.Client(credentials=sa.Credentials.from_service_account_file(_key),
                      project='", ee_project, "')

def add_cloud_bands(img):
    cld_prb = ee.Image(img.get('s2cloudless')).select('probability')
    is_cloud = cld_prb.gt(CLD_PRB_THRESH).rename('clouds')
    return img.addBands(ee.Image([cld_prb, is_cloud]))

def add_shadow_bands(img):
    not_water = img.select('SCL').neq(6)
    dark_pixels = (img.select('B8').lt(NIR_DRK_THRESH*SR_BAND_SCALE)
                     .multiply(not_water).rename('dark_pixels'))
    shadow_azimuth = ee.Number(90).subtract(ee.Number(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')))
    cld_proj = (img.select('clouds')
                  .directionalDistanceTransform(shadow_azimuth, CLD_PRJ_DIST*10)
                  .reproject(**{'crs': img.select(0).projection(), 'scale': 100})
                  .select('distance').mask().rename('cloud_transform'))
    shadows = cld_proj.multiply(dark_pixels).rename('shadows')
    return img.addBands(ee.Image([dark_pixels, cld_proj, shadows]))

def add_cld_shdw_mask(img):
    img_cloud = add_cloud_bands(img)
    img_cloud_shadow = add_shadow_bands(img_cloud)
    is_cld_shdw = img_cloud_shadow.select('clouds').add(
                  img_cloud_shadow.select('shadows')).gt(0)
    is_cld_shdw = (is_cld_shdw.focalMin(2).focalMax(BUFFER*2/20)
                     .reproject(**{'crs': img.select([0]).projection(), 'scale': 20})
                     .rename('cloudmask'))
    return img.addBands(is_cld_shdw)

def apply_cld_shdw_mask(img):
    not_cld_shdw = img.select('cloudmask').Not()
    scl = img.select('SCL')
    # Scene Classification Layer codes: 1 saturated, 2 dark area, 3 cloud shadow,
    # 8 cloud medium probability, 9 cloud high probability, 10 thin cirrus,
    # 11 snow/ice. Already there and can help filter cloud edges.
    scl_mask = (scl.neq(1).And(scl.neq(2)).And(scl.neq(3))
                  .And(scl.neq(8)).And(scl.neq(9))
                  .And(scl.neq(10)).And(scl.neq(11)))
    return img.select('B.*').updateMask(not_cld_shdw).updateMask(scl_mask)

def mask_saturated(img):
    saturated = img.select(['B2','B3','B4','B8']).lt(65535).reduce(ee.Reducer.min())
    return img.updateMask(saturated)

def mask_low_reflectance(img):
    return (img.updateMask(img.select('B8').gt(100))
               .updateMask(img.select('B4').gt(100))
               .updateMask(img.select('B2').gt(30))
               .updateMask(img.select('B3').gt(100)))

def add_indices(img):
    NIR   = img.select('B8').divide(SR_BAND_SCALE)
    RED   = img.select('B4').divide(SR_BAND_SCALE)
    BLUE  = img.select('B2').divide(SR_BAND_SCALE)
    GREEN = img.select('B3').divide(SR_BAND_SCALE)
    NDVI = NIR.subtract(RED).divide(NIR.add(RED)).rename('NDVI')
    EVI  = (NIR.subtract(RED)
              .divide(NIR.add(RED.multiply(6)).subtract(BLUE.multiply(7.5)).add(1))
              .multiply(2.5).rename('EVI'))
    MSAVI2 = (NIR.multiply(2).add(1)
                .subtract(NIR.multiply(2).add(1).pow(2)
                          .subtract(NIR.subtract(RED).multiply(8)).sqrt())
                .divide(2).rename('MSAVI2'))
    NIRv  = NDVI.multiply(NIR).rename('NIRv')
    SR    = NIR.divide(RED).rename('SR')
    DVI   = NIR.subtract(RED).rename('DVI')
    EVI2  = (NIR.subtract(RED).divide(NIR.add(RED.multiply(2.4)).add(1))
               .multiply(2.5).rename('EVI2'))
    MSR   = SR.subtract(1).divide(SR.add(1).sqrt()).rename('MSR')
    kNDVI = NDVI.pow(2).tanh().rename('kNDVI')
    GCVI  = NIR.divide(GREEN).subtract(1).rename('GCVI')
    NDGI  = GREEN.subtract(RED).divide(GREEN.add(RED)).rename('NDGI')
    return img.addBands(ee.Image([NDVI,EVI,MSAVI2,NIRv,SR,DVI,EVI2,MSR,kNDVI,GCVI,NDGI]))

def add_time(img):
    return img.addBands(ee.Image(img.getNumber('system:time_start')).rename('timestamp'))

def add_atmos(img, aoi):
    # AOT and WVP are Sen2Cor outputs carried in L2A, both scaled by 0.001. Kept
    # as scene properties so haze can be filtered on without exporting the bands.
    a = img.select('AOT').multiply(0.001)
    w = img.select('WVP').multiply(0.001)
    st = a.addBands(w).reduceRegion(reducer=ee.Reducer.mean(), geometry=aoi,
                                    scale=60, maxPixels=int(1e9), bestEffort=True)
    return img.set('aot_mean', st.get('AOT'), 'wvp_mean', st.get('WVP'))

def one_per_day(ic):
    # Where a site is imaged twice in a day from adjacent orbits, keep the granule
    # closer to nadir: sort ascending by view zenith and take the first per date. Default = off
    ic = ic.map(lambda im: im.set('day', ee.Date(im.get('system:time_start'))
                                          .format('YYYY-MM-dd')))
    days = ic.aggregate_array('day').distinct()
    def pick(d):
        return (ic.filter(ee.Filter.eq('day', d))
                  .sort('MEAN_INCIDENCE_ZENITH_ANGLE_B8A').first())
    return ee.ImageCollection(days.map(pick))

def add_valid_frac(img, aoi):
    # fraction of the AOI still unmasked after the full chain above. (specific band doesn't matter)
    v = img.select('B4').mask().rename('v')
    f = v.reduceRegion(reducer=ee.Reducer.mean(), geometry=aoi, scale=20,
                       maxPixels=int(1e9), bestEffort=True).get('v')
    return img.set('valid_frac', f)

def build(aoi, d0, d1, maxcloud):
    s2 = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
            .filterBounds(aoi).filterDate(d0, d1)
            .filterMetadata('CLOUDY_PIXEL_PERCENTAGE','less_than', maxcloud))
    cl = (ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
            .filterBounds(aoi).filterDate(d0, d1))
    j = ee.ImageCollection(ee.Join.saveFirst('s2cloudless').apply(
        primary=s2, secondary=cl,
        condition=ee.Filter.equals(leftField='system:index', rightField='system:index')))
    
    # apply_cld_shdw_mask keeps only B.* , so the indices are computed after masking, on masked reflectance.
    # dropping by aerosol first saves running the whole mask chain on data that will not be used.
    out = (j.map(lambda im: add_atmos(im, aoi))
             .filter(ee.Filter.lte('aot_mean', AOT_MAX))
             .map(add_cld_shdw_mask).map(apply_cld_shdw_mask)
             .map(mask_saturated).map(mask_low_reflectance)
             .map(lambda im: add_valid_frac(im, aoi))
             .filter(ee.Filter.gte('valid_frac', MIN_VALID_FRAC))
             .map(add_indices).map(add_time).sort('system:time_start'))
    return one_per_day(out).sort('system:time_start') if ONE_PASS_A_DAY else out

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

years <- seq(as.integer(format(as.Date(start_date), "%Y")),
             as.integer(format(as.Date(end_date), "%Y")))

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

# One site-year is usually one export, but years with too many images have to be
# split or GEE rejects it. Parts are named <key>_p1, _p2 ... and are stitched
# back together downstream by date so it doesn't matter for scripts 3b-4
# The site-year .done marker is only written once every part has landed.
parts_left <- new.env(parent = emptyenv())

submit_window <- function(key, site, d0, d1, done_path, depth = 0L) {
  aoi <- aoi_of(site)
  ic  <- ee$ImageCollection(py$build(aoi, format(d0), format(d1), max_cloud))
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

  task <- py$submit(ic, as.list(keep_bands), aoi, out_scale, out_crs,
                    bucket_name, key)
  inflight[[length(inflight) + 1L]] <<-
    list(task = task, key = key, done_path = done_path,
         dates = mdf[, c('layer', 'stamp')], meta = mdf)
  
  message(sprintf('submitted %-14s %3d images x %d bands', key, n, length(keep_bands)))
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
