################################################################################
#### 5b - Daily climate grids and site timeseries, from AGCD               ####
################################################################################
# Australian Gridded Climate Data at 0.05 deg daily: rainfall, maximum and
# minimum temperature. AGCD v1 is the Australian Water Availability Project
# (AWAP) analysis renamed and continued, so this is the same product the earlier
# BOM grid downloads were drawing on, from its maintained source.
#
# Release v1-0-4 covers 1900-01-01 to 2025-12-31. Note the other AGCD releases
# on NCI do not: v1 stops at 2020-05-31, and v2-0-4 is monthly only.
#
# Pulled from the NCI THREDDS NetCDF Subset Service, which crops to the study
# bounding box server-side - a Pilbara year is ~7 MB against ~50 MB for the
# national file. Open access, no account or key.
#   https://geonetwork.nci.org.au/geonetwork/srv/eng/catalog.search#/metadata/f4742_3902_2481_1714
#
# Grids are cached in data/climGrids, one netCDF per variable per year, and the
# per-site daily series written to outputs/site_climate_AGCD.rds.

library(terra)

grid_dir  <- "data/climGrids"
MVG_shp   <- "data/Vegetation NVIS classification 15-3-23.shp"
end_path  <- "outputs"
out_file  <- file.path(end_path, "site_climate_AGCD.rds")
years     <- 2015:as.integer(format(Sys.Date(), "%Y"))
bbox      <- c(north = -21.0, south = -23.8, west = 116.5, east = 120.5)
release   <- "v1-0-4"
overwrite <- FALSE

# variable -> its folder under zv2/agcd/<release>, which is also its file stem
vars <- c(precip = "precip/total", tmax = "tmax/mean", tmin = "tmin/mean")

dir.create(grid_dir, showWarnings = FALSE, recursive = TRUE)

agcd_get <- function(v, y) {
  stem <- sub("/", "_", vars[[v]])
  f    <- file.path(grid_dir, sprintf("agcd_%s_%d_pilbara.nc", stem, y))
  if (file.exists(f) && !overwrite) return(f)
  u <- sprintf(paste0("https://thredds.nci.org.au/thredds/ncss/grid/zv2/agcd/%s/",
                      "%s/r005/01day/agcd_v1_%s_r005_daily_%d.nc",
                      "?var=%s&north=%s&south=%s&west=%s&east=%s",
                      "&horizStride=1&time=all&accept=netcdf"),
               release, vars[[v]], stem, y, v,
               bbox["north"], bbox["south"], bbox["west"], bbox["east"])
  # years past the end of the release 404, which is expected and reported below
  ok <- suppressWarnings(try(utils::download.file(u, f, mode = "wb", quiet = TRUE),
                             silent = TRUE))
  if (inherits(ok, "try-error") || !file.exists(f) || file.size(f) < 1e4) {
    if (file.exists(f)) unlink(f)
    return(NULL)
  }
  f
}

# site points, one per NVIS site centroid
pts <- centroids(terra::aggregate(project(vect(MVG_shp), "EPSG:4326"), by = "Site"))

series <- lapply(names(vars), function(v) {
  out <- do.call(rbind, lapply(years, function(y) {
    f <- agcd_get(v, y)
    if (is.null(f)) { message("  ", v, " ", y, " unavailable"); return(NULL) }
    r <- rast(f); crs(r) <- "EPSG:4326"
    m <- terra::extract(r, pts, ID = FALSE)
    message("  ", v, " ", y, ": ", nlyr(r), " days")
    data.frame(Site = rep(pts$Site, each = nlyr(r)),
               Date = rep(as.Date(time(r)), times = nrow(m)),
               value = as.numeric(t(as.matrix(m))), stringsAsFactors = FALSE)
  }))
  if (!is.null(out)) names(out)[names(out) == "value"] <- v
  out
})
names(series) <- names(vars)
series <- Filter(Negate(is.null), series)
if (!length(series)) stop("no AGCD data retrieved - check the network")

clim <- Reduce(function(a, b) merge(a, b, by = c("Site", "Date"), all = TRUE), series)
clim <- clim[order(clim$Site, clim$Date), ]
names(clim)[names(clim) == "precip"] <- "Rainfall"
names(clim)[names(clim) == "tmax"]   <- "Tmax_dC"
names(clim)[names(clim) == "tmin"]   <- "Tmin_dC"

dir.create(end_path, showWarnings = FALSE)
saveRDS(split(clim, clim$Site), out_file)
message(length(unique(clim$Site)), " sites, ", as.character(min(clim$Date)), " to ",
        as.character(max(clim$Date)), ", written to ", out_file)
