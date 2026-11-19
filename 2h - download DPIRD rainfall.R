################################################################################
#### 2h - Rainfall from the DPIRD Weather API v2                            ####
################################################################################
# The two DPIRD-network stations near the study sites, both DBCA's: MILL,
# Millstream Chichester National Park (beside DC, CS, DR and MS) and KARI,
# Karijini (30 to 100 km from the Hamersley sites; it stopped reporting in
# December 2023). Daily rainfall is the 24 hours to 9 am. Hourly summaries run
# from the station start; 15 minute summaries only from the start of the previous
# year, which is all the API keeps. Data are DPIRD's, CC BY 3.0 AU. The key is
# read from data/bom_data and sent as a header, never written out.

library(httr2)

key_file <- "data/bom_data/dpirdApiKey.txt"
out_dir  <- "data/bom_data/dpird"
stations <- c("MILL", "KARI")
base     <- "https://api.agric.wa.gov.au/v2/weather"
fields   <- "period,rainfall,observations,errors,stationAvailability"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
key <- trimws(readLines(key_file, warn = FALSE)[1])

# 10 requests a second and 4000 an hour; a page holds 1000 records at most
dpird <- function(path, query = list()) {
  Sys.sleep(0.15)
  request(paste0(base, path)) |>
    req_headers(`API-KEY` = key, Accept = "application/json") |>
    req_url_query(!!!query) |>
    req_retry(max_tries = 6, backoff = function(i) 2^i,
              is_transient = function(r) resp_status(r) %in% c(429, 500, 502, 503, 504)) |>
    req_perform() |> resp_body_json()
}

num <- function(x) if (is.null(x) || is.list(x)) NA_real_ else as.numeric(x)

# one window of summaries as a data frame; availability is a number for hourly
# and 15 minute records and three for daily ones, of which the 24 hours to 9 am
# is kept
fetch_window <- function(station, kind, from, to) {
  q <- if (kind == "daily") list(startDate = format(from), endDate = format(to))
       else list(startDateTime = format(from, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
                 endDateTime   = format(to,   "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
  rec <- list(); offset <- 0L
  repeat {
    r <- dpird(sprintf("/stations/%s/summaries/%s", station, kind),
               c(q, list(select = fields, limit = 1000, offset = offset)))
    s <- r$data$summaries
    rec <- c(rec, s)
    offset <- offset + length(s)
    if (!length(s) || offset >= r$metadata$collection$count) break
  }
  if (!length(rec)) return(NULL)
  data.frame(
    from_utc     = vapply(rec, function(s) s$period$from, ""),
    to_utc       = vapply(rec, function(s) s$period$to, ""),
    year         = vapply(rec, function(s) num(s$period$year), 0),
    month        = vapply(rec, function(s) num(s$period$month), 0),
    day          = vapply(rec, function(s) num(s$period$day), 0),
    rainfall_mm  = vapply(rec, function(s) num(s$rainfall), 0),
    obs_count    = vapply(rec, function(s) num(s$observations$count), 0),
    obs_pct      = vapply(rec, function(s) num(s$observations$percentage), 0),
    errors       = vapply(rec, function(s) num(s$errors), 0),
    availability = vapply(rec, function(s)
      num(if (is.list(s$stationAvailability)) s$stationAvailability$to9AM else s$stationAvailability), 0),
    stringsAsFactors = FALSE)
}

# windows of `days` days, end points inclusive, each small enough for one page
pull <- function(station, kind, from, to, days) {
  step   <- if (kind == "daily") days else days * 86400
  starts <- seq(from, to, by = step)
  ends   <- pmin(starts + (step - 1), to)
  out <- do.call(rbind, lapply(seq_along(starts), function(i) fetch_window(station, kind, starts[i], ends[i])))
  if (is.null(out)) return(NULL)
  out <- out[!duplicated(out$from_utc), ]
  out[order(out$from_utc), ]
}

t_now <- Sys.time(); attr(t_now, "tzone") <- "UTC"
t_15  <- as.POSIXct(sprintf("%d-01-01", as.integer(format(t_now, "%Y")) - 1), tz = "UTC")
cols  <- c("from_utc", "to_utc", "rainfall_mm", "obs_count", "obs_pct", "errors", "availability")
write_out <- function(x, station, kind) {
  if (is.null(x)) return(message(station, " ", kind, ": no records"))
  write.csv(x, file.path(out_dir, sprintf("%s_%s_rainfall.csv", station, kind)), row.names = FALSE)
  message(station, " ", kind, ": ", nrow(x), " records, ", min(x$from_utc), " to ", max(x$from_utc))
}

for (station in stations) {
  st <- dpird(paste0("/stations/", station))$data
  message(st$stationCode, " ", st$stationName, " (", st$latitude, ", ", st$longitude,
          "), record from ", st$startDate)
  write.csv(data.frame(code = st$stationCode, name = st$stationName, lat = st$latitude,
                       lon = st$longitude, altitude_m = st$altitude, owner = st$ownerCode,
                       start = st$startDate, end = if (is.null(st$endDate)) NA else st$endDate,
                       gauge_height_m = num(st$rainGaugeHeight), logger = st$model),
            file.path(out_dir, paste0(station, "_station.csv")), row.names = FALSE)

  t_start <- as.POSIXct(paste(st$startDate, "00:00:00"), tz = "UTC") - 86400

  daily <- pull(station, "daily", as.Date(t_start), as.Date(t_now), days = 1000)
  daily$date <- as.Date(sprintf("%04d-%02d-%02d", daily$year, daily$month, daily$day))
  write_out(daily[, c("date", cols)], station, "daily")
  write_out(pull(station, "hourly", t_start, t_now, days = 40)[, cols], station, "hourly")
  q15 <- pull(station, "15min", t_15, t_now, days = 10)
  write_out(if (!is.null(q15)) q15[, cols], station, "15min")
}
