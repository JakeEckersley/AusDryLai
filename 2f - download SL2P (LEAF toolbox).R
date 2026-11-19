################################################################################
#### 7 - LAI from the CCRS LEAF-Toolbox (SL2P) via the Python Earth Engine ####
################################################################################
# Downloads one stack per site year (every date as bands).
# Script 3c assembles and filters these stacks (domain flag, one retrieval per
# date). The s2cloudless probability is exported as a fourth band but not used:
# the SL2P arm keeps the toolbox's own QA60 and SCL masks. layout:
#   stacks/<SITE>_<start>_<end>.tif   the annual stacks
#   dates/<SITE>_<start>_<end>.csv    acquisition times, image ids and view zenith
#                                     angles, one row per date
#   done/<SITE>_<year>.done           completion markers (for safe restarts)

# Both SL2P resolutions are downloaded to separate folders:
#   10 m networks use B2, B3, B4, B8        -> LeafDownloads_10m
#   20 m networks use B3-B7, B8A, B11, B12  -> LeafDownloads_20m
# Both are Sentinel-2; the Landsat networks are separate dictionary entries.

# Stacks download to <name>.part and are renamed on success
# site-year is skipped once the .done marker exists. Deleting a .done marker
# forces that site-year to be re-checked; deleting the stack forces a re-download.

# LOCAL FIX TO THE VENDORED TOOLBOX - re-apply it after any re-clone.
# vendor/LEAF-Toolbox/Source-Python/production/LEAF.py, makeProductCollection:
# setDefaultProjection() took its CRS from image.bandNames().slice(0,1), which is
# band 0 of the Sentinel-2 image - B1, the 60 m aerosol band. That fixed the
# information content at 60 m and the following .reproject() merely laid it onto
# the requested grid, so exported SL2P had a 60 m change period at BOTH
# resolutions (6 px at 10 m, 3 px at 20 m, measured on the raw downloads).
# Changed to use `projection`, the network's first reflectance band, which is
# what the line above it already computes. Anything downloaded before
# 2026-10-13 carries the 60 m artefact and needs re-downloading.
#
# needed to run:
#   git clone --depth 1 https://github.com/rfernand387/LEAF-Toolbox vendor/LEAF-Toolbox
#   "C:/Users/eck025/.conda/envs/s2/python.exe" -m pip install tqdm

library(terra)
library(reticulate)
library(jsonlite)
library(future.apply)

## ---- config -----------------------------------------------------------------
leaf_dir    <- normalizePath("vendor/LEAF-Toolbox/Source-Python/production",
                             winslash = "/", mustWork = FALSE)
py_ver      <- "C:/Users/eck025/.conda/envs/s2/python.exe"
ee_project  <- "sent-dryland-validation"
serv_email  <- "jakeeckersley@sent-dryland-validation.iam.gserviceaccount.com"
serv_key    <- "D:/PilbaraLaiProject/serv-key/sent-dryland-validation-e2e627d1d6e2.json"

end_path    <- "outputs"
sites_path  <- "data/SiteBoundsForSatDownloads"
mvg_path    <- "data/Vegetation NVIS classification 15-3-23.shp"

start_date  <- "2015-03-01"                # full Sentinel-2 record
end_date    <- "2024-09-01"
variable    <- "LAI"                       # or fAPAR, fCOVER, Albedo
max_cloud   <- 50                          # CLOUD_FILTER
out_crs     <- "EPSG:32750"                # same as the Sentinel-2 downloads
min_valid_frac <- 0.80                     # keep a date only if at least this much of
                                           # the site AOI survives LEAF's masking.
                                           # Every tile delivery is kept; script 3c
                                           # collapses them and applies the QC band.
n_workers   <- 6                           # parallel EE sessions

runs <- list(
  list(tag = "10m", scale = 10L, dir = "D:/PilbaraLaiProject/LeafDownloads_10m",
       collection = "COPERNICUS/S2_SR_HARMONIZED_10m"),
  list(tag = "20m", scale = 20L, dir = "D:/PilbaraLaiProject/LeafDownloads_20m",
       collection = "COPERNICUS/S2_SR_HARMONIZED"))

# The 10 m entry needed two bugs in LEAF's dictionary patched. Both are
# not in the networks:
#   -  "tools" is wired to toolsHLS, which selects an Fmask band Sentinel-2 does
#      not have - should be toolsS2.
#   -  inputOffset has 8 entries where inputBands and inputScaling have 7, so
#      scaleBands() does a 7-band .add() against an 8-element list. GEE
#      returned "Images must contain the same number of bands. Got 7 and 8".
#      offsets were 0, so trimming was identical. Affects
#      all seven network variables at 10 m, not just LAI.

## ---- what to run ------------------------------------------------------------
paired    <- readRDS(file.path(end_path, "paired_ref_and_sat_dataset.rds"))
recs      <- unique(paired[, c("Site", "MVG", "Date")])
recs$Date <- as.Date(recs$Date)

# Downloads every site that has a bounds shapefile
# paired_ref_and_sat_dataset.rds covers only 15 sites - calibration/validation only uses the paired ones.
site_list <- sort(tools::file_path_sans_ext(
  basename(Sys.glob(file.path(sites_path, "*.shp")))))
years     <- seq(as.integer(format(as.Date(start_date), "%Y")),
                 as.integer(format(as.Date(end_date),   "%Y")))

# passes bounding box as GeoJSON text, same as "3 - DownloadSentinel2viaGCS.ipynb":
# FeatureCollection -> geometry() -> bounds()... gee freaks out with too many vertices
# needed to be text not SpatVector for safe multi-worker processing.
site_json <- vapply(site_list, function(s) {
  shp <- file.path(sites_path, paste0(s, ".shp"))
  if (!file.exists(shp)) return(NA_character_)
  tmp <- tempfile(fileext = ".geojson"); on.exit(unlink(tmp), add = TRUE)
  writeVector(project(vect(shp), "EPSG:4326"), tmp, filetype = "GeoJSON",
              overwrite = TRUE)
  paste(readLines(tmp, warn = FALSE), collapse = "")
}, character(1))
site_list <- site_list[!is.na(site_json[site_list])]

for (r in runs) {
  dir.create(r$dir, showWarnings = FALSE, recursive = TRUE)
  for (sub in c("stacks", "dates", "done")) dir.create(file.path(r$dir, sub), showWarnings = FALSE, recursive = TRUE)
}

## ---- download the site-year stacks -----------------------------------------
# One job per site: the worker builds its own Python session (reticulate objects
# cannot cross process boundaries) and then loops years x resolutions, so the
# ~15 s Earth Engine start-up is paid once per site rather than once per stack.
do_site <- function(k, site_list, site_json, runs, years, cfg) {
  library(reticulate); library(terra)
  reticulate::use_python(cfg$py_ver, required = TRUE)
  reticulate::py_run_string("import sys")
  reticulate::py_run_string(paste0("sys.path.insert(0, r'", cfg$leaf_dir, "')"))
  ee <- reticulate::import("ee")
  ee$Initialize(credentials = ee$ServiceAccountCredentials(cfg$serv_email, cfg$serv_key),
                project = cfg$ee_project)
  LEAF   <- reticulate::import("LEAF")
  dicts  <- reticulate::import("dictionariesSL2P")
  SL2PV1 <- reticulate::import("SL2PV1")
  reticulate::py_run_string("
import ee
def stamps(ic):
    return ic.aggregate_array('system:time_start').map(
        lambda t: ee.Date(t).format('YYYY_MM_dd_HH_mm_ss')).getInfo()
def geometry(ic):
    return {'id': ic.aggregate_array('system:index').getInfo(),
            'vza': ic.aggregate_array('MEAN_INCIDENCE_ZENITH_ANGLE_B8A').getInfo()}
def year_url(ic, bands, aoi, scale, crs):
    return (ic.select(bands).toBands().toFloat()
              .getDownloadURL({'scale': scale, 'region': aoi, 'crs': crs,
                               'format': 'GEO_TIFF', 'filePerBand': False}))
def keep_full_scenes(ic, aoi, band, scale, minfrac):
    # Drops dates where LEAF's masking has removed more than 1 - minfrac of the sites
    def add(img):
        v = img.select(band).mask().rename('v')
        fr = v.reduceRegion(reducer=ee.Reducer.mean(), geometry=aoi, scale=scale,
                            maxPixels=int(1e9), bestEffort=True).get('v')
        return img.set('valid_frac', fr)
    return ic.map(add).filter(ee.Filter.gte('valid_frac', minfrac))
def band_names(ic, bands):
    return ic.select(bands).toBands().bandNames().getInfo()
")
  site  <- site_list[k]
  aoi   <- ee$FeatureCollection(jsonlite::fromJSON(site_json[[site]],
                                  simplifyVector = FALSE))$geometry()$bounds()
  bands <- list(paste0("estimate", cfg$variable), paste0("error", cfg$variable), "QC",
                "s2cloudless_probability")
  msg   <- character()

  for (r in runs) {
    co <- dicts$make_collection_options(SL2PV1)[[r$collection]]
    no <- dicts$make_net_options()[[cfg$variable]][[r$collection]]
    if (grepl("_10m$", r$collection)) {                      # patch the 10 m entry
      co[["tools"]] <- reticulate::import("toolsS2")
      nb <- length(reticulate::py_to_r(no[["inputBands"]]))
      no[["inputOffset"]] <- reticulate::r_to_py(as.list(rep(0, nb)))
    }
    # Fetch one date window: download the stack and cut it into per-date tiles.
    # Returns TRUE on success. Sites in a double orbit overlap get ~270 images a
    # year instead of ~135, and toBands() exceeds a GEE request limit
    # Window is halved on failure and retried, down to about six weeks.
    fetch <- function(d0, d1, depth = 0L) {
      lab <- sprintf("%s %s %s..%s", site, r$tag, format(d0), format(d1))
      ic <- try(ee$ImageCollection(LEAF$makeProductCollection(
                  co, no, cfg$variable, aoi, format(d0), format(d1),
                  cfg$max_cloud, r$scale)), silent = TRUE)
      if (!inherits(ic, "try-error") && min_valid_frac > 0)
        ic <- try(ee$ImageCollection(reticulate::py$keep_full_scenes(
                    ic, aoi, paste0("estimate", cfg$variable), r$scale,
                    min_valid_frac)), silent = TRUE)
      n <- if (inherits(ic, "try-error")) NULL else
        try(reticulate::py_to_r(ic$size()$getInfo()), silent = TRUE)
      if (is.null(n) || inherits(n, "try-error")) {
        if (depth >= 3L) { msg <<- c(msg, paste(lab, "collection failed")); return(FALSE) }
      } else {
        if (!n) return(TRUE)                                  # nothing in window

        f <- file.path(r$dir, "stacks",
                       sprintf("%s_%s_%s.tif", site, format(d0, "%Y%m%d"),
                               format(d1, "%Y%m%d")))
        ok <- file.exists(f) && file.size(f) > 0
        if (!ok) {
          
          # transient 503s happen under load, so retries with a backoff. Downloads
          # to .part and renamed, to stop crash leaving truncated stacks
          p <- paste0(f, ".part")
          for (attempt in 1:3) {
            u <- try(reticulate::py$year_url(ic, bands, aoi, r$scale, cfg$out_crs),
                     silent = TRUE)
            if (!inherits(u, "try-error")) {
              ok <- tryCatch({ download.file(u, p, mode = "wb", quiet = TRUE); TRUE },
                             error = function(e) FALSE)
              if (ok && file.exists(p) && file.size(p) > 0) break
            }
            unlink(p); ok <- FALSE
            Sys.sleep(10 * attempt)
          }
          if (ok) file.rename(p, f)
        }
        if (ok) {
          # download strips band names, so the acquisition times are written
          # beside the stack as a dates .csv. toBands() emits the four bands per
          # image in collection order - estimateLAI, errorLAI, QC, s2cloudless
          # probability - so row i of the sidecar is stack layers 4*(i-1) + 1:4.
          
          d <- file.path(r$dir, "dates", sub(".tif", ".csv", basename(f), fixed = TRUE))
          if (!file.exists(d)) {
            st <- try(reticulate::py_to_r(reticulate::py$stamps(ic)), silent = TRUE)
            g  <- try(reticulate::py_to_r(reticulate::py$geometry(ic)), silent = TRUE)
            x  <- try(terra::rast(f), silent = TRUE)
            if (inherits(st, "try-error") || inherits(g, "try-error") || inherits(x, "try-error") ||
                terra::nlyr(x) != length(st) * length(bands)) {
              unlink(f)                                       # unusable, redo
              if (depth >= 3L) { msg <<- c(msg, paste(lab, "metadata failed")); return(FALSE) }
            } else {
              utils::write.csv(data.frame(layer = seq_along(st), stamp = st,
                                          id = unlist(g$id), vza = unlist(g$vza)),
                               d, row.names = FALSE)
              return(TRUE)
            }
          } else return(TRUE)
        }
        if (depth >= 3L) { msg <<- c(msg, paste(lab, "failed after halving")); return(FALSE) }
      }
      mid <- d0 + floor(as.numeric(d1 - d0) / 2)
      if (mid <= d0 || mid >= d1) { msg <<- c(msg, paste(lab, "cannot split")); return(FALSE) }
      a <- fetch(d0, mid, depth + 1L)
      b <- fetch(mid, d1, depth + 1L)
      a && b
    }

    for (y in years) {
      done <- file.path(r$dir, "done", sprintf("%s_%d.done", site, y))
      if (file.exists(done)) next
      y0 <- max(as.Date(cfg$start_date), as.Date(sprintf("%d-01-01", y)))
      y1 <- min(as.Date(cfg$end_date),   as.Date(sprintf("%d-01-01", y + 1)))
      if (y1 <= y0) { file.create(done); next }
      if (fetch(y0, y1)) file.create(done)
    }
  }
  c(sprintf("%s done", site), msg)
}

cfg <- list(py_ver = py_ver, leaf_dir = leaf_dir, serv_email = serv_email,
            serv_key = serv_key, ee_project = ee_project, variable = variable,
            start_date = start_date, end_date = end_date, max_cloud = max_cloud,
            out_crs = out_crs)

message(length(site_list), " sites x ", length(years), " years x ", length(runs),
        " resolutions on ", n_workers, " workers")
plan(multisession, workers = n_workers)
log <- future_lapply(seq_along(site_list), do_site, site_list = site_list,
                     site_json = site_json, runs = runs, years = years, cfg = cfg,
                     future.seed = TRUE)
plan(sequential)
writeLines(unlist(log))


################################################################################
# NOTES
#
# Output. One annual stack per site per resolution under <dir>/stacks/, named
# SITE_<start>_<end>.tif. Each has 4 bands per date - estimateLAI, errorLAI, QC and
# the s2cloudless probability (0-100, unused) - and the band names carry the
# source image id, e.g.
#   20210101T021349_20210101T022139_T50KQV_estimateLAI
# so the acquisition date of every layer is recoverable from the name alone. No
# per-date tiles are written: script 3c assembles the per-site temporal stacks
# straight from these annual stacks.
#
# SL2PV1 is the current production algorithm module; SL2PV0 is the
# original SNAP-equivalent and SL2PV1_1 a later revision. Swapping the import
# swaps the network, which is a cheap sensitivity test worth running.
################################################################################
