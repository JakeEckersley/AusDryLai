################################################################################
#### 3g - NDVI flag against the band and wetting-day flags, with rainfall    ####
################################################################################
# One figure per site compares the flagging processes against its daily rainfall:
#   NDVI flag     script 3b: each pixel's NDVI against the median of its dates within
#                 screen_days either side (limits mad_k)
#   band flags    highs in B4 and B8 on any date, lows in B4 and B11 within days_buffer
#                 days of a wetting day, and a date dropped when more than scene_frac of
#                 its valid pixels spike
#   wetting days  SMAP L4 daily rises above the threshold in script 3a
# Rainfall is the gauge nearest the site where one was downloaded (a BOM Climate Data
# Online daily rainfall file in data/bom_data, or a DPIRD gauge of script 2h: MILL, at
# Millstream, for DC, CS, DR and MS, and KARI, at Karijini, for JD, MCK2 and FRS),
# otherwise the AGCD 5 km series at the site.
# Daily rainfall is the 24 h to 9 am, so rain on the afternoon of day d is recorded on
# d + 1: a flag on day d is matched to rain recorded on d or d + 1. The mask, the NDVI
# screen and the date gate are read from script 3b so the two cannot drift apart.
#
# The NDVI series is also smoothed with a Savitzky-Golay (SG) filter on a 2-day grid.
# Two tables, written to outputs/, give the share of a 5 to 30 d pulse that each SG
# window keeps and the error in predicting a held-out observation, by how many
# neighbours the observation has within 15 d.

library(terra)
library(ggplot2)
library(patchwork)
library(future.apply)

sites_plot  <- c("FRW", "CW", "FRE", "KC")   # FRW is rainfall driven, 43 km from Karijini North; CW is 6 km from Coolawanyah; FRE and KC are 25 km from Newman Aero
rain_thr    <- c(1, 5)                # mm, for the agreement table
days_buffer <- 3                      # days after a wetting day (the day itself and the next n) in which lows are flagged
circle_min  <- 0.05                   # circle a date when at least this share of its valid pixels is NDVI-flagged
grey_min    <- 0.10                   # draw the raw site mean on dates with at least this share of the pixels valid
grid_days   <- 2; max_fill <- 30      # SG runs on a 2-day grid; gaps over max_fill days stay empty
sg_variants <- list("22 d" = c(11, 11), "30 d" = c(15, 15), "38 d" = c(19, 19), "46 d" = c(23, 23), "widened" = c(15, 31))   # cubic SG windows, in grid points (shortest, longest)
wide_obs    <- 6                      # the widened window is the shortest that holds this many observations
pulse_days  <- c(5, 10, 15, 20, 30)   # pulse widths (full width at half maximum), days
n_workers   <- min(8L, parallelly::availableCores(omit = 2))
refresh     <- TRUE                   # FALSE re-plots from the cached flags and SG tables
bom_dir     <- "data/bom_data/unzipped"
gauge_files <- c(CW  = "coolawanyah_IDCJAC0009_005001_1800/IDCJAC0009_005001_1800_Data.csv",
                 FRW = "karijiniNorth_IDCJAC0009_005098_1800/IDCJAC0009_005098_1800_Data.csv",
                 KC  = "newmanAero_IDCJAC0009_007176_1800/IDCJAC0009_007176_1800_Data.csv",    # 25 km
                 FRE = "newmanAero_IDCJAC0009_007176_1800/IDCJAC0009_007176_1800_Data.csv")
gauge_files <- setNames(file.path(bom_dir, gauge_files), names(gauge_files))
gauge_names <- c(CW = "Coolawanyah 005001", FRW = "Karijini North 005098", KC = "Newman Aero 007176", FRE = "Newman Aero 007176")
mill_sites  <- c("DC", "CS", "DR", "MS")      # within 5 km of MILL
kari_sites  <- c("JD", "MCK2", "FRS")         # 29, 34 and 51 km from KARI, which reported to December 2023
gauge_files <- c(gauge_files,
                 setNames(rep("data/bom_data/dpird/MILL_daily_rainfall.csv", length(mill_sites)), mill_sites),
                 setNames(rep("data/bom_data/dpird/KARI_daily_rainfall.csv", length(kari_sites)), kari_sites))
gauge_names <- c(gauge_names,
                 setNames(rep("Millstream MILL (DPIRD)", length(mill_sites)), mill_sites),
                 setNames(rep("Karijini KARI (DPIRD)", length(kari_sites)), kari_sites))
panels <- c(B4 = "B4 red", B8 = "B8 NIR", B11 = "B11 SWIR1")      # the bands the band flags use
fig_path <- "plots/filteringComparisons"; end_path <- "outputs"
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

# rules from 3b, dates analysed in 5b
take <- function(file, nms) for (e in parse(file, keep.source = FALSE))
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% nms) eval(e, globalenv())
take("3b - VI timeseries.R", c("mask_apply", "mask_dilate_m", "one_pass_a_day", "min_valid_frac", "vi_valid", "sat_root",
                               "mad_k", "screen_days", "screen_side", "zdev"))
take("5b - phenology timeseries analysis.R", "obs_range")
scene_frac <- 1 - min_valid_frac      # a date is dropped when more than this share of its valid pixels spike
source("R/stack_io.R")
stack_idx <- stack_index(sat_root)
DailySM   <- readRDS("data/SmDailyTimeseries.RDS")
agcd      <- readRDS(file.path(end_path, "site_climate_AGCD.rds"))
# the wetting threshold in script 3a, for the titles
wet_mad <- as.numeric(sub(".*> *([0-9.]+).*", "\\1", grep("high_change *<-", readLines("3a - soil moisture timeseries.R"), value = TRUE)[1]))

# wet soil flag for pixel timeseries
# Spike detector for pixel timeseries (ts), applied rowwise over pixel x dates matrix
# points flagged when it departs from the nearest valid dates either side in the same direction
# aimed at finding isolated cloud shadow or soil darkening after rain (low) and bright
# haze or cloud edges (high); fire and green-up persist across dates, so are safe.
#   dates    acquisition dates (same length and order as ts - s2md_sub$Date
#   mag      flag is based on MAD (median dif between dates), diff has to exceed (mad * mag) to be flagged
#   max_gap  (days) from the point to each neighbouring valid observation; flagging is skipped beyond it
qc_flag <- function(ts, dates = NULL, direction = "low", mag = 1.5, max_gap = 10) {
  n <- length(ts); out <- rep(FALSE, n)
  v <- which(!is.na(ts)); m <- length(v)
  if (m < 3) return(out)
  x <- ts[v]
  d_prev <- x - c(NA, x[-m])
  d_next <- x - c(x[-1], NA)
  s <- mad(c(d_prev, d_next), na.rm = TRUE)
  ok <- !is.na(d_prev) & !is.na(d_next) & is.finite(s) & s > 0

  # skip flag step if either neighbouring valid observation is further than max_gap days away
  if (!is.null(max_gap)) {
    if (is.null(dates)) stop("qc_flag: dates are required when max_gap is set")
    if (length(dates) != n) stop("qc_flag: dates and ts differ in length")
    d <- as.numeric(as.Date(dates))[v]
    ok <- ok & c(NA, diff(d)) <= max_gap & c(diff(d), NA) <= max_gap
  }
  ok[is.na(ok)] <- FALSE

  low_qc  <- ok & d_prev / s < -mag & d_next / s < -mag
  high_qc <- ok & d_prev / s >  mag & d_next / s >  mag

  if (direction == "low") out[v] <- low_qc
  else if (direction == "high") out[v] <- high_qc
  else if (direction == "both") out[v] <- high_qc | low_qc
  else stop('direction must be "low", "high", or "both"')
  out
}

# used to identify if two bands are spiking the direction, as a pixels x dates matrix
spike_flag <- function(a_mat, b_mat, dates, direction) {
  t(apply(a_mat, 1, qc_flag, dates = dates, direction = direction)) &
    t(apply(b_mat, 1, qc_flag, dates = dates, direction = direction))
}

# Climate Data Online daily rainfall. A row that closes an accumulation over several
# days is set to NA with the days it covers, since the daily timing is unknown. A DPIRD
# day (script 2h) is set to NA when under 95% of its minutes were recorded.
read_gauge <- function(f) {
  g <- read.csv(f, check.names = FALSE)
  if ("rainfall_mm" %in% names(g))
    return(data.frame(Date = as.Date(g$date), Rainfall = ifelse(g$obs_pct >= 95, g$rainfall_mm, NA)))
  r <- data.frame(Date = as.Date(sprintf("%s-%s-%s", g$Year, g$Month, g$Day)), Rainfall = g[["Rainfall amount (millimetres)"]])
  per <- g[["Period over which rainfall was measured (days)"]]
  for (i in which(!is.na(per) & per > 1)) r$Rainfall[max(1, i - per[i] + 1):i] <- NA
  r
}
site_rain <- function(site, src = c("gauge", "AGCD")) {
  if (match.arg(src) == "gauge" && site %in% names(gauge_files))
    return(list(rain = read_gauge(gauge_files[[site]]), label = gauge_names[[site]]))
  list(rain = agcd[[site]][, c("Date", "Rainfall")], label = "AGCD")
}

## ---- agreement of the wetting flag with rainfall, every site ------------------
# precision: flagged days with rain >= thr recorded on that day or the next
# recall:    rain days >= thr recorded on r that carry a flag on r - 1 or r
sm_by_site <- split(DailySM, DailySM$Site_code)
agree <- do.call(rbind, lapply(names(sm_by_site), function(s) {
  sm <- sm_by_site[[s]]
  do.call(rbind, lapply(c("AGCD", if (s %in% names(gauge_files)) "gauge"), function(src) {
    r <- site_rain(s, src)$rain; r <- r[!is.na(r$Rainfall), ]
    d <- data.frame(Date = seq(max(min(sm$Date), min(r$Date)), min(max(sm$Date), max(r$Date)), by = 1))
    d$flag <- sm$high_change[match(d$Date, sm$Date)]; d$rain <- r$Rainfall[match(d$Date, r$Date)]
    d$rain_next <- c(d$rain[-1], NA); d$flag_prev <- c(NA, d$flag[-nrow(d)])
    do.call(rbind, lapply(rain_thr, function(th) {
      wet <- d$rain >= th | d$rain_next >= th
      data.frame(Site = s, source = src, rain_mm = th, flag_days = sum(d$flag, na.rm = TRUE), flag_pct = 100 * mean(d$flag, na.rm = TRUE),
                 precision = 100 * mean(wet[which(d$flag)], na.rm = TRUE), base_rate = 100 * mean(wet[!is.na(d$flag)], na.rm = TRUE),
                 recall = 100 * mean((d$flag | d$flag_prev)[which(d$rain >= th)], na.rm = TRUE))
    }))
  }))
}))
write.csv(agree, file.path(end_path, "wetting_flag_vs_rain.csv"), row.names = FALSE)
cat("median over all sites, AGCD rainfall:\n")
print(aggregate(cbind(flag_pct, precision, base_rate, recall) ~ rain_mm, agree[agree$source == "AGCD", ], median), digits = 3)
cat("\nsites with a gauge (flag_pct is of days with a gauge record):\n")
g_sites <- unique(agree$Site[agree$source == "gauge"])
print(agree[agree$Site %in% g_sites, ], digits = 3, row.names = FALSE)

## ---- the flags for one site ----------------------------------------------------
# The pixel-dates the NDVI screen flags, below and above, for a pixels x dates matrix.
ndvi_screen <- function(V, dates) {
  dn <- as.numeric(dates); plan(multisession, workers = n_workers); on.exit(plan(sequential))
  Z <- do.call(rbind, future_lapply(lapply(seq_len(nrow(V)), function(i) V[i, ]), zdev, d = dn))
  list(below = !is.na(Z) & Z < -mad_k[["below"]], above = !is.na(Z) & Z > mad_k[["above"]])
}

# per date, band flags: the share of valid pixels spiking low, high, and flagged (a high, or a low
# inside a wetting window), and whether the date is dropped as a scene-level spike; and, for
# each of the band panels, the pixel mean and standard deviation before the flags (dates with
# at least grey_min of the pixels valid), and the mean after the flagged pixel-dates are
# removed (dates the band flags keep: min_valid_frac valid and unflagged, not a scene-level spike).
# per date, NDVI flag (nd_): the share of valid pixels flagged below and above, the pixel mean
# and standard deviation before the flag (dates with at least grey_min of the pixels valid), and
# the mean after it (dates kept by the min_valid_frac gate).
# Returned with the pixel x date NDVI matrix after the masks, and its flags.
site_flags <- function(site) {
  md <- site_meta(site, stack_idx)
  band <- function(b) site_band_series(site, b, stack_idx)
  B4 <- band("B4"); B8 <- band("B8"); B11 <- band("B11")
  stopifnot(nlyr(B4) == nrow(md), nlyr(B8) == nlyr(B4), nlyr(B11) == nlyr(B4))
  q    <- site_qc_series(site, stack_idx)
  bits <- c(c(cloud_dea = 1L, cloud_prob = 2L, fmask = 4L, contig = 8L, terrain = 16L)[names(mask_apply)[mask_apply]], nodata = 32L)
  cloudy <- c("cloud_dea", "cloud_prob", "fmask")
  soft <- sum(bits[names(bits) %in% cloudy]); hard <- sum(bits[!names(bits) %in% cloudy])
  bad  <- if (soft > 0) bits_set(q, soft) else NULL
  if (!is.null(bad) && mask_dilate_m > 0) bad <- focal(bad, w = 2 * ceiling(mask_dilate_m / res(B4)[1]) + 1, fun = "max", na.rm = TRUE) > 0
  if (hard > 0) { h <- bits_set(q, hard); bad <- if (is.null(bad)) h else (bad | h) }
  B4 <- mask(B4, bad, maskvalues = 1); B8 <- mask(B8, bad, maskvalues = 1); B11 <- mask(B11, bad, maskvalues = 1)
  sel <- seq_len(nrow(md))
  if (one_pass_a_day) {
    cov <- 1 - colMeans(is.na(values(B4)))
    sel <- sort(unname(vapply(split(seq_along(cov), md$Date), function(j) j[order(-cov[j], md$MEAN_INCIDENCE_ZENITH_ANGLE_B8A[j])][1], integer(1))))
  }
  mats <- list(B4 = values(B4[[sel]]) / 1e4, B8 = values(B8[[sel]]) / 1e4, B11 = values(B11[[sel]]) / 1e4)   # DEA scale 10000
  dates <- md$Date[sel]
  lo <- spike_flag(mats$B4, mats$B11, dates, "low"); hi <- spike_flag(mats$B4, mats$B8, dates, "high"); spike <- lo | hi
  sm  <- DailySM[DailySM$Site_code == site, ]
  win <- vapply(dates, function(d) any(sm$high_change[sm$Date %in% (d - 0:days_buffer)]), logical(1))
  flagged <- hi | sweep(lo, 2, win, "&")
  nv  <- pmax(colSums(!is.na(mats$B4)), 1)
  fl  <- data.frame(Date = dates, n_pix = nrow(mats$B4), n_valid = colSums(!is.na(mats$B4)), low = colSums(lo) / nv, high = colSums(hi) / nv,
                    spike = colSums(spike) / nv, flagged = colSums(flagged) / nv, in_window = win)
  fl$dir <- ifelse(fl$low >= fl$high, "low", "high")
  fl$scene_drop <- fl$spike > scene_frac
  V <- values(mask(site_vi_series(site, "NDVI", stack_idx, lyrs = sel), bad[[sel]], maskvalues = 1))
  V[!is.na(V) & (V < vi_valid$NDVI[1] | V > vi_valid$NDVI[2])] <- NA
  mats$NDVI <- V
  sc <- ndvi_screen(V, dates); nflag <- sc$below | sc$above; nvv <- colSums(!is.na(V)); grey <- nvv / nrow(V) >= grey_min
  Vn <- V; Vn[nflag] <- NA
  fl$nd_valid <- nvv; fl$nd_low <- colSums(sc$below) / pmax(nvv, 1); fl$nd_high <- colSums(sc$above) / pmax(nvv, 1); fl$nd_share <- fl$nd_low + fl$nd_high
  fl$nd_raw <- ifelse(grey, colMeans(V, na.rm = TRUE), NA); fl$nd_sd <- ifelse(grey, apply(V, 2, sd, na.rm = TRUE), NA)
  fl$nd_flt <- ifelse(colMeans(!is.na(V) & !nflag) >= min_valid_frac, colMeans(Vn, na.rm = TRUE), NA)
  grey_b <- fl$n_valid / fl$n_pix >= grey_min
  fl$kept_band <- fl$n_valid * (1 - fl$flagged) / fl$n_pix >= min_valid_frac & !fl$scene_drop
  for (v in names(mats)) {
    V <- mats[[v]]; Vf <- V; Vf[flagged] <- NA
    fl[[paste0(v, "_mean")]] <- ifelse(grey_b, colMeans(V, na.rm = TRUE), NA)
    fl[[paste0(v, "_sd")]]   <- ifelse(grey_b, apply(V, 2, sd, na.rm = TRUE), NA)
    fl[[paste0(v, "_flt")]]  <- ifelse(fl$kept_band, colMeans(Vf, na.rm = TRUE), NA)
  }
  list(fl = fl, ndvi = mats$NDVI, screen = sc)
}

## ---- SG windows on the screened NDVI, one site -----------------------------------------
# Screen: the NDVI flags above, with a date kept when min_valid_frac of the site's pixels are valid and
# unflagged. SG: the kept points are filled linearly to a 2-day grid (gaps over max_fill days stay empty),
# smoothed by a cubic SG of the given window and read back at the observation dates; "widened" takes, at
# each grid point, the shortest window of its range that holds wide_obs observations. Pulse test: a
# Gaussian pulse of 0.05 NDVI added to real, screened pixel series in the dry and the wet season; a linear
# smoother keeps the same share at any amplitude, and observations the screen removed count as lost.
# Held out: 60 observations per pixel, each predicted from the rest of its series.
sg_n <- seq(min(unlist(sg_variants)), max(unlist(sg_variants)), by = 2)
Fs <- lapply(setNames(sg_n, sg_n), function(n) { k <- (n - 1) / 2; t <- -k:k; X <- outer(t, 0:3, `^`); X %*% solve(crossprod(X), t(X)) })   # row i: the cubic fit at position i of an n-point window
sg <- function(x, F) {   # centred weights inside, end windows fitted at the edges (as signal::sgolayfilt)
  n <- nrow(F); k <- (n - 1) / 2; m <- length(x)
  y <- as.numeric(stats::filter(x, F[k + 1, ], sides = 2))
  y[1:k] <- F[1:k, , drop = FALSE] %*% x[1:n]; y[(m - k + 1):m] <- F[(k + 2):n, , drop = FALSE] %*% x[(m - n + 1):m]
  y
}
sg_fit <- function(x, d, n_min, n_max = n_min, at = d) {   # SG of the series x (days d) on a 2-day grid, read at `at`
  v <- which(!is.na(x)); out <- rep(NA_real_, length(at)); K <- length(v); if (K < max(3, wide_obs)) return(out)
  dv <- d[v]; g <- seq(dv[1], by = grid_days, length.out = ceiling((dv[K] - dv[1]) / grid_days) + 1); N <- length(g); if (N < n_max) return(out)
  y <- approx(dv, x[v], xout = g, rule = 2)$y
  j <- pmin(findInterval(g, dv, rightmost.closed = TRUE), K - 1L); gap <- (dv[j + 1L] - dv[j]) > max_fill & g > dv[j] & g < dv[j + 1L]
  nn <- seq(n_min, n_max, by = 2); S <- vapply(nn, function(n) sg(y, Fs[[as.character(n)]]), numeric(N))
  if (length(nn) == 1) ys <- S[, 1] else {   # per grid point, blend the two window lengths either side of the one that holds wide_obs observations
    cand <- outer(findInterval(g, dv), (-(wide_obs - 1)):wide_obs, `+`); ok <- cand >= 1 & cand <= K; cand[!ok] <- 1L
    D <- abs(matrix(dv[cand], nrow = N) - g); D[!ok] <- Inf; H <- apply(D, 1, function(r) sort(r, partial = wide_obs)[wide_obs])
    pos <- (pmin(pmax(2 * H / grid_days + 1, n_min), n_max) - n_min) / 2 + 1; lo <- pmin(floor(pos), length(nn) - 1L); w <- pos - lo
    ys <- S[cbind(seq_len(N), lo)] * (1 - w) + S[cbind(seq_len(N), lo + 1L)] * w
  }
  ys[gap] <- NA
  i <- pmin(pmax(findInterval(at, g, rightmost.closed = TRUE), 1L), N - 1L); u <- (at - g[i]) / grid_days
  out <- ys[i] * (1 - u) + ys[i + 1L] * u; out[at < g[1] | at > g[N]] <- NA; out
}
pulse_keep <- function(b, dn) {   # b: a screened pixel series y and a pulse centre c0; by pulse width, the pulse energy and the part of it each fit keeps
  f0 <- lapply(sg_variants, function(s) sg_fit(b$y, dn, s[1], s[2]))
  t(sapply(pulse_days, function(w) {
    pulse <- 0.05 * exp(-4 * log(2) * (dn - b$c0)^2 / w^2); y1 <- b$y + pulse; u <- which(!is.na(b$y) & abs(dn - b$c0) <= max(2 * w, 10))
    c(den = sum(pulse[u]^2), mapply(function(s, f) { f1 <- sg_fit(y1, dn, s[1], s[2]); k <- u[!is.na(f1[u]) & !is.na(f[u])]; sum((f1[k] - f[k]) * pulse[k]) }, sg_variants, f0))
  }))
}
loo_pixel <- function(x2, dn) {   # hold out 60 observations in turn and predict each from the rest of the series
  v <- which(!is.na(x2)); k <- length(v); if (k < 60) return(NULL)
  do.call(rbind, lapply(v[sample(12:(k - 11), min(60, k - 22))], function(j) { xl <- x2; xl[j] <- NA
    c(dens = sum(!is.na(x2) & abs(dn - dn[j]) <= 15) - 1, truth = x2[j], vapply(sg_variants, function(s) sg_fit(xl, dn, s[1], s[2], at = dn[j]), numeric(1))) }))
}
sg_window_check <- function(V, dates, sc) {   # V: pixels x dates NDVI and sc its flags; the flag shares, and the pulse and held-out results
  ix <- dates >= obs_range[1] & dates <= obs_range[2]; V <- V[, ix, drop = FALSE]; dimnames(V) <- NULL; dn <- as.numeric(dates[ix])
  below <- sc$below[, ix, drop = FALSE]; above <- sc$above[, ix, drop = FALSE]
  set.seed(8); plan(multisession, workers = n_workers); on.exit(plan(sequential))
  kept <- colMeans(!is.na(V) & !(below | above)) >= min_valid_frac
  X2 <- V; X2[below | above] <- NA; X2[, !kept] <- NA
  pxs <- which(rowSums(!is.na(V)) >= 300)      # pixels with a long record
  centres <- as.numeric(as.Date(c(paste0(2019:2023, "-06-01"), paste0(2018:2022, "-12-01"), paste0(2019:2023, "-02-01"), paste0(2019:2023, "-04-01"))))   # mid dry season; Dec, Feb and Apr of the wet season
  b <- data.frame(p = sample(pxs, 210, replace = TRUE), c0 = sample(centres, 210, replace = TRUE) + runif(210, -3, 3))
  pk <- Reduce(`+`, future_lapply(lapply(seq_len(nrow(b)), function(i) list(y = X2[b$p[i], ], c0 = b$c0[i])), pulse_keep, dn = dn)); rownames(pk) <- pulse_days
  cv <- do.call(rbind, future_lapply(lapply(sample(pxs, 200), function(p) X2[p, ]), loo_pixel, dn = dn, future.seed = TRUE))
  nv <- sum(!is.na(V))
  list(screen = c(pixels = nrow(V), dates = ncol(V), flagged = 100 * sum(below | above) / nv, below = 100 * sum(below) / nv, above = 100 * sum(above) / nv,
                  dates_valid = sum(colMeans(!is.na(V)) >= min_valid_frac), dates_kept = sum(kept)), pulse = pk, loo = cv)
}

## ---- figure ----------------------------------------------------------------------
check_plot <- function(site, fl, xlim) {
  sm <- DailySM[DailySM$Site_code == site, ]; rr <- site_rain(site)
  r <- rr$rain[rr$rain$Date >= xlim[1] & rr$rain$Date <= xlim[2], ]
  s <- sm[sm$Date >= xlim[1] & sm$Date <= xlim[2], ]; f <- fl[fl$Date >= xlim[1] & fl$Date <= xlim[2], ]
  wet <- s[s$high_change, ]
  ax <- list(scale_x_date(limits = xlim, expand = c(0.005, 0)), theme_bw(base_size = 11), theme(axis.title.x = element_blank()))
  rings <- function(nm) scale_size_area(nm, max_size = 6, limits = c(0, 1), breaks = c(0.05, 0.1, 0.2, 0.4, 0.8))
  p1 <- ggplot(r, aes(Date, Rainfall)) + geom_col(width = 1, fill = "grey30") + ax +
    labs(y = sprintf("Rainfall (mm)\n%s", rr$label),
         title = sprintf("%s: NDVI flag against the band and wetting-day flags\nSMAP wetting days (red dots): daily rise > %s MAD (%d of %d days); blue bands = 0-%d d after one, the window in which band lows are flagged\nNDVI flag: more than %g MADs below or %g above the median of the pixel's dates within %d d either side (at least %d a side); a date is dropped when more than %.0f%% of its pixels are invalid or flagged\nGrey dots = raw site mean, black dots = mean after that panel's flags; red line (NDVI panel) = mean after the band flags; circles = dates with at least %.0f%% of the pixels NDVI-flagged\nBand flags: highs in B4 and B8 on any date; lows in B4 and B11 within 0-%d d of a wetting day (> %s MAD from the nearest clear dates within %s d); rings = dates dropped because more than %.0f%% of the valid pixels spike",
                         site, wet_mad, nrow(wet), nrow(s), days_buffer, mad_k[["below"]], mad_k[["above"]], screen_days, screen_side, 100 * (1 - min_valid_frac),
                         100 * circle_min, days_buffer, formals(qc_flag)$mag, formals(qc_flag)$max_gap, 100 * scene_frac)) + theme(plot.title = element_text(size = 9))
  win <- geom_rect(data = data.frame(xmin = wet$Date, xmax = wet$Date + days_buffer + 1), inherit.aes = FALSE,
                   aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "steelblue", alpha = 0.12)
  p2 <- ggplot(s, aes(Date, SoilMoisture)) + win + geom_line(colour = "steelblue4", linewidth = 0.5) + ax +
    geom_point(data = wet, colour = "red", size = 1.5) + labs(y = "SMAP surface soil\nmoisture (m3/m3)")
  nd <- f[!is.na(f$nd_raw), ]; big_nd <- nd[nd$nd_share >= circle_min, ]
  p3 <- ggplot(nd, aes(Date)) +
    geom_ribbon(aes(ymin = nd_raw - nd_sd, ymax = nd_raw + nd_sd), fill = "grey70", alpha = 0.6) +
    geom_line(aes(y = NDVI_flt), colour = "red", linewidth = 0.4) +
    geom_point(aes(y = nd_raw), colour = "grey60", size = 1.5) + geom_point(aes(y = nd_flt), colour = "black", size = 0.9, na.rm = TRUE) +
    geom_point(data = big_nd, aes(y = nd_raw, size = nd_share), shape = 1, colour = "red", stroke = 0.9) + rings("share of pixels\nNDVI-flagged") + ax + labs(y = "NDVI")
  p4 <- ggplot(f) + geom_col(aes(Date, nd_share), width = 1, fill = "darkorange") + geom_col(aes(Date, nd_low), width = 1, fill = "steelblue4") +
    geom_hline(yintercept = 1 - min_valid_frac, linetype = "dashed") + ax +
    labs(y = "NDVI flag: share\nof pixels below (blue)\nabove (orange)") + theme(axis.title.y = element_text(size = 8))
  big <- f[f$scene_drop, ]
  ts_panel <- function(v) {
    m <- paste0(v, "_mean"); d <- f[!is.na(f[[m]]), ]
    g <- ggplot(d, aes(Date)) +
      geom_ribbon(aes(ymin = .data[[m]] - .data[[paste0(v, "_sd")]], ymax = .data[[m]] + .data[[paste0(v, "_sd")]]), fill = "grey70", alpha = 0.6) +
      geom_point(aes(y = .data[[m]]), colour = "grey60", size = 1.5) + geom_point(aes(y = .data[[paste0(v, "_flt")]]), colour = "black", size = 0.9, na.rm = TRUE) +
      ax + labs(y = panels[[v]])
    bd <- big[!is.na(big[[m]]), ]
    if (nrow(bd)) g <- g + geom_point(data = bd, aes(y = .data[[m]], shape = dir, size = spike), colour = "red", fill = NA, stroke = 0.8) +
      scale_shape_manual("dominant spike", values = c(low = 6, high = 2)) + rings("share of pixels\nspiking")
    g
  }
  pend <- ggplot(f) + geom_col(aes(Date, spike), width = 1, fill = "darkorange") + geom_col(aes(Date, flagged), width = 1, fill = "red") +
    geom_hline(yintercept = scene_frac, linetype = "dashed") + ax +
    labs(y = "Band flags: share\nof pixels spiking\n(orange) flagged (red)") + theme(axis.title.y = element_text(size = 8))
  Reduce(`/`, c(list(p1, p2, p3, p4), lapply(names(panels), ts_panel), list(pend))) +
    plot_layout(heights = c(1.7, 1.3, 2, 0.8, 1.1, 1.1, 1.1, 0.9), guides = "collect")
}

for (site in sites_plot) {
  cache <- file.path(end_path, paste0("wetting_flag_check_", site, ".rds")); sg_cache <- file.path(end_path, paste0("sg_window_check_", site, ".rds"))
  if (refresh || !file.exists(cache) || !file.exists(sg_cache)) { sf <- site_flags(site); saveRDS(sf$fl, cache); saveRDS(sg_window_check(sf$ndvi, sf$fl$Date, sf$screen), sg_cache) }
  fl <- readRDS(cache)
  inside <- function(x) sprintf("%.0f", 100 * sum(x[fl$in_window] * fl$n_valid[fl$in_window]) / sum(x * fl$n_valid))
  message(site, ": ", nrow(fl), " dates; ", sprintf("%.1f", 100 * mean(fl$in_window)), "% inside a post-wetting window; spike pixel-dates inside one: ",
          inside(fl$low), "% of low, ", inside(fl$high), "% of high; share of all pixel-dates, low ", sprintf("%.2f", 100 * sum(fl$low * fl$n_valid) / sum(fl$n_valid)),
          "%, high ", sprintf("%.2f", 100 * sum(fl$high * fl$n_valid) / sum(fl$n_valid)), "%, flagged ", sprintf("%.2f", 100 * sum(fl$flagged * fl$n_valid) / sum(fl$n_valid)), "%")
  big <- fl$scene_drop       # scene-level dates, dropped whole by the band flags
  message("  ", sum(big), " of ", nrow(fl), " dates are dropped as scene-level spikes (> ", 100 * scene_frac, "% of pixels; ", sum(big & fl$in_window),
          " inside a window); they hold ", sprintf("%.0f", 100 * sum(fl$spike[big] * fl$n_valid[big]) / sum(fl$spike * fl$n_valid)),
          "% of the spiking pixel-dates; pixel-dates flagged on the other dates: ",
          sprintf("%.2f", 100 * sum((fl$flagged * fl$n_valid)[!big]) / sum(fl$n_valid[!big])), "%")
  message("  NDVI flag: ", sprintf("%.2f", 100 * sum(fl$nd_share * fl$nd_valid) / sum(fl$nd_valid)), "% of the valid pixel-dates (",
          sprintf("%.2f", 100 * sum(fl$nd_low * fl$nd_valid) / sum(fl$nd_valid)), "% below, ", sprintf("%.2f", 100 * sum(fl$nd_high * fl$nd_valid) / sum(fl$nd_valid)), "% above)")
  # dates with imagery flags, soil moisture and rainfall all present
  ok   <- site_rain(site)$rain$Date[!is.na(site_rain(site)$rain$Rainfall)]
  span <- c(max(min(fl$Date), min(ok)), min(max(fl$Date), max(ok)))
  g <- check_plot(site, fl, span)
  nm <- paste0("Filtering comparison ", site)
  ggsave(g, filename = file.path(fig_path, paste0(nm, ".tiff")), width = 40, height = 34, units = "cm")
}

## ---- SG window tables, every site ------------------------------------------------------
sgc <- lapply(setNames(sites_plot, sites_plot), function(s) readRDS(file.path(end_path, paste0("sg_window_check_", s, ".rds"))))
cat(sprintf("\nNDVI screen: %g MADs below, %g above the median of the points within %d d either side (at least %d points a side); flagged, below and above are %% of the valid pixel-dates; dates_valid, dates_kept: dates with at least %.0f%% of the pixels valid before the screen, and valid and unflagged after it\n",
            mad_k[["below"]], mad_k[["above"]], screen_days, screen_side, 100 * min_valid_frac))
print(round(do.call(rbind, lapply(sgc, `[[`, "screen")), 2))
pulse_tab <- function(s) {
  p <- Reduce(`+`, lapply(sgc[s], `[[`, "pulse")); m <- t(round(100 * p[, names(sg_variants)] / p[, "den"], 1)); colnames(m) <- paste(pulse_days, "d")
  data.frame(site = if (length(s) > 1) "All" else s, window = names(sg_variants), m, row.names = NULL, check.names = FALSE)
}
loo_tab <- function(s) {
  cv <- do.call(rbind, lapply(sgc[s], `[[`, "loo")); cv <- cv[complete.cases(cv), ]
  cls <- cut(cv[, "dens"], c(-1, 2, 4, Inf), labels = c("sparse (0-2)", "medium (3-4)", "dense (5+)"))
  do.call(rbind, lapply(levels(cls), function(l) { e <- cv[cls == l, names(sg_variants), drop = FALSE] - cv[cls == l, "truth"]
    data.frame(site = if (length(s) > 1) "All" else s, neighbours = l, n = nrow(e), stat = c("MAE", "RMSE"), round(rbind(colMeans(abs(e)), sqrt(colMeans(e^2))), 4), check.names = FALSE) }))
}
pulse_all <- do.call(rbind, c(lapply(sites_plot, pulse_tab), list(pulse_tab(sites_plot))))
loo_all   <- do.call(rbind, c(lapply(sites_plot, loo_tab), list(loo_tab(sites_plot))))
write.csv(pulse_all, file.path(end_path, "sg_window_pulse_kept.csv"), row.names = FALSE)
write.csv(loo_all, file.path(end_path, "sg_window_loo_error.csv"), row.names = FALSE)
cat("\nShare of a pulse kept (%) by each SG window, by pulse width (full width at half maximum): a 0.05 NDVI Gaussian pulse added to real screened pixel series, dry and wet season\n")
print(pulse_all, row.names = FALSE)
cat("\nMean absolute error (NDVI) in predicting a held-out observation, by how many other observations lie within 15 d of it (RMSE in the csv)\n")
print(loo_all[loo_all$stat == "MAE", names(loo_all) != "stat"], row.names = FALSE)
