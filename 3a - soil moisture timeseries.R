library(terra)
library(future)
library(future.apply)
library(progressr)

smap_files <- list.files("data/spl4smau", full.names = TRUE)
site_maps  <- Sys.glob(file.path("data/SiteBoundsForSatDownloads", "*.shp"))
site_codes <- tools::file_path_sans_ext(basename(site_maps))
n_cores    <- max(1L, parallelly::availableCores() - 5L)

# timestamps from filenames; Etc/GMT-8 is AWST
stamp     <- sub("^.*_", "", tools::file_path_sans_ext(basename(smap_files)))
TimeStamp <- as.POSIXct(stamp, format = "%Y%m%dT%H%M%S", tz = "UTC")
attr(TimeStamp, "tzone") <- "Etc/GMT-8"
LocalDate <- as.Date(TimeStamp, tz = "Etc/GMT-8")   # as.Date() defaults to UTC

# site cells (one shared grid across all granules)
ref           <- rast(smap_files[1], lyrs = "sm_surface_analysis")
cells_by_site <- lapply(site_maps, function(s) cells(ref, ext(project(vect(s), ref))))
need <- sort(unique(unlist(cells_by_site)))
pick <- lapply(cells_by_site, function(ci) match(ci, need))

read_chunk <- function(ii, files, need, pick, p) {
  out <- t(vapply(files[ii], function(f) {
    v <- terra::values(terra::rast(f, lyrs = "sm_surface_analysis"), mat = FALSE)[need]
    vapply(pick, function(k) mean(v[k], na.rm = TRUE), numeric(1))
  }, numeric(length(pick)), USE.NAMES = FALSE))
  p(amount = length(ii))   # workers buffer ticks until the chunk resolves
  out
}

handlers("progress")
plan(multisession, workers = n_cores)
chunks <- split(seq_along(smap_files), cut(seq_along(smap_files), n_cores * 4, labels = FALSE))
mat <- with_progress({
  p <- progressor(steps = length(smap_files))
  do.call(rbind, future_lapply(chunks, read_chunk,
                               files = smap_files, need = need, pick = pick, p = p))
})
plan(sequential)

n <- nrow(mat)
SubDaily_ts <- data.frame(
  Site_code    = rep(site_codes, each = n),
  TimeStamp    = rep(TimeStamp,  times = ncol(mat)),
  Date         = rep(LocalDate,  times = ncol(mat)),
  SoilMoisture = as.vector(mat)
)

Daily_ts <- aggregate(SoilMoisture ~ Site_code + Date, data = SubDaily_ts, FUN = mean, na.rm = TRUE)
Daily_ts <- Daily_ts[order(Daily_ts$Site_code, Daily_ts$Date), ]
Daily_ts$prev_day    <- Daily_ts$SoilMoisture[match(paste(Daily_ts$Site_code, Daily_ts$Date - 1),
                                                    paste(Daily_ts$Site_code, Daily_ts$Date))]
Daily_ts$sm_dif      <- Daily_ts$SoilMoisture - Daily_ts$prev_day
# wetting day: the daily rise is more than 4 median absolute deviations of the
# site's daily changes over its whole record (6.9% of site-days)
Daily_ts$norm_change <- ave(Daily_ts$sm_dif, Daily_ts$Site_code,
                            FUN = function(z) (z - median(z, na.rm = TRUE)) / mad(z, na.rm = TRUE))
Daily_ts$high_change <- !is.na(Daily_ts$norm_change) & Daily_ts$norm_change > 4

saveRDS(SubDaily_ts, "data/SmSubDailyTimeseries.RDS")
saveRDS(Daily_ts,    "data/SmDailyTimeseries.RDS")

site_code <- "BN2"
sub <- SubDaily_ts[SubDaily_ts$Site_code == site_code, ]
day <- Daily_ts[Daily_ts$Site_code == site_code, ]
with(sub, plot(TimeStamp, SoilMoisture, main = site_code))
with(day, plot(Date, SoilMoisture, type = "l", main = site_code))
with(day[day$high_change, ], points(Date, SoilMoisture, col = "red", pch = 19))

day <- Daily_ts[Daily_ts$Site_code == site_code & Daily_ts$Date > '2022-01-01', ]
with(sub, plot(TimeStamp, SoilMoisture, main = site_code))
with(day, plot(Date, SoilMoisture, type = "l", main = site_code))
with(day[day$high_change, ], points(Date, SoilMoisture, col = "red", pch = 19))
