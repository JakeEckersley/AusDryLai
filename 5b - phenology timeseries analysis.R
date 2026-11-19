################################################################################
#### 5b - LAI growth phases and climate predictors                          ####
################################################################################
# Three LAI sources, treated identically: the filtered index with the lowest
# held-out RMSE (NDVI or EVI, chosen in script 5a) through its fitted lm, SL2P
# 10 m, and SL2P 10 m calibrated against the in situ data. The VI and SL2P stacks
# share a grid, so all three describe the same pixels. Smoothing is applied to
# LAI rather than the index, since SL2P has no index. Climate windows are found
# by date arithmetic against the daily series. A rain event is a run of
# consecutive days of at least event_min_day mm, so lighter rain neither starts
# nor chains events.
#
# Writes outputs/LAI_phase_stats_by_source.rds, one row per phase per pixel per
# source, read by script 5c.

library(terra)
library(future.apply)
library(parallelly)

## ---------------------------------------------------------------- config ----
proj_crs   <- "EPSG:4326"
S2_path    <- "D:/PilbaraLaiProject/FliteredTimeseries"   # filtered VI stacks
raw_path   <- "D:/PilbaraLaiProject/RawTimeseries_Sen2Cor"  # unfiltered Sen2Cor VI stacks, script 4
leaf_path  <- "D:/PilbaraLaiProject/FilteredLeaf"      # assembled and QC-filtered by script 3c
MVG_shp    <- "data/Vegetation NVIS classification 15-3-23.shp"
rip_shp    <- "data/riparian vegetation boundary edited 2025.shp"
end_path   <- "outputs"
plot_path  <- "plots"
txt_size   <- 11
dir.create(plot_path, showWarnings = FALSE)

obs_range  <- as.Date(c("2018-01-01", "2024-09-01"))  # observations to read (S2 is sparse before 2018)
sg_days    <- 40                      # Savitzky-Golay window, days
sg_order   <- 3                       # polynomial order
sg_step    <- 2                       # the series is linearly filled to a grid this many days apart
max_fill   <- 30                      # gaps longer than this (days) stay empty
min_obs    <- 40                      # pixels with fewer clear dates are skipped
min_amp    <- 0.025                   # LAI. Turning points are strict local
                                      # extrema, so a ripple of 0.01 creates a
                                      # peak/trough pair and a spurious phase.
                                      # A phase must change LAI by more than
                                      # this, up or down; 0 keeps every turning
                                      # point. 0.025 is about 0.004 EVI, under the
                                      # 0.01 EVI amplitude threshold of Broich et
                                      # al. (2014, Biogeosciences 11, 5181-5198).
lags       <- seq(0, 56, by = 7)      # days before the phase, for the lagged vars
event_min_day <- 1                    # mm. A day under this belongs to no rain event: AGCD
                                      # spreads trace rain over many days the gauges
                                      # record dry. 0 counts any day with rain.
n_large    <- 3                       # events kept apart in each window: the largest (LRE),
                                      # RE2, ... RE<n_large>, and the rain left once the
                                      # 1, 2, ... n_large largest are taken out (MRE, MRE2, ...)
sites      <- NULL                   # NULL = every site; or a character vector

# the sites held out of the models fitted in script 5a, each with the one NVIS
# class used for its timeseries figure (woodland, shrubland, hummock grassland and
# Melaleuca), which has reference records in both campaigns
test_sites <- c("MCK", "CW", "YU", "DR")
ts_mvg     <- c(MCK = "Eucalypt Woodlands",
                CW  = "Acacia Shrublands",
                YU  = "Hummock Grasslands",
                DR  = "Melaleuca Forests and Woodlands")

# index -> LAI model chosen in script 5a, and the plausible range of each index.
# EVI's denominator can approach zero, so single pixels reach values that the
# log model turns into absurd LAI; masked before conversion, not after.
vi_choice <- readRDS(file.path(end_path, "Optimum_closest_VI_choice.rds"))
vi_valid  <- list(NDVI = c(-1, 1), EVI = c(-1, 2.5))

## ------------------------------------------------------------- the models ---
# keyed <index>_<form>; the unfiltered series has its own fit, key raw<index>_<form>
lin_mods <- readRDS(file.path(end_path, "Optimum_closest_VI_lin_reg_LAI_mod.rds"))
stopifnot(c(vi_choice$key, paste0("raw", vi_choice$key)) %in% names(lin_mods))
message("VI model: ", vi_choice$key, " (held-out RMSE ", vi_choice$test_RMSE, ")")

# reference LAI paired with the satellite values, filtered as in script 5a
paired <- readRDS(file.path(end_path, "paired_ref_and_sat_dataset.rds"))
paired <- paired[!(paired$n_records_understorey < 8 | paired$n_records_canopy < 8 |
                   paired$n_pixels < 8), ]

# SL2P calibration: reference LAI on the 10 m retrieval over the training sites.
# Mirrors the fit in script 5a - keep the two in step if either changes.
leaf_cal <- local({
  d <- paired[!paired$Site %in% test_sites, ]               # training sites only
  d <- d[is.finite(d$leaf_10m_closest) & is.finite(d$LAI_LX_total), ]
  lm(LAI_LX_total ~ leaf_10m_closest, data = d)
})
message("SL2P calibration: LAI = ", round(coef(leaf_cal)[1], 3), " + ",
        round(coef(leaf_cal)[2], 3), " x retrieval (n = ", nobs(leaf_cal), ")")

# what each source needs: where the stack lives and how it becomes LAI
sources <- c(
  setNames(list(list(kind = "vi", stem = vi_choice$stem, dir = S2_path,
                     model = vi_choice$key)), vi_choice$stem),
  list(`SL2P 10 m`            = list(kind = "lai", stem = "leaf_10m", dir = leaf_path),
       `SL2P 10 m calibrated` = list(kind = "cal", stem = "leaf_10m", dir = leaf_path)))

site_list <- if (is.null(sites))
  basename(list.dirs(S2_path, recursive = FALSE)) else sites

## ------------------------------------------------------------- extraction ---
# One site, one source -> long data frame of pixel x date LAI, carrying the
# vegetation class and landform of each pixel. Pixels outside the NVIS polygons
# are dropped, which is most of the raster.
extract_site <- function(site, src_name, srcs = sources) {
  src <- srcs[[src_name]]
  f <- file.path(src$dir, site, paste0(src$stem, "-ts.tif"))
  if (!file.exists(f)) return(NULL)
  r <- try(rast(f), silent = TRUE); if (inherits(r, "try-error")) return(NULL)

  d <- as.Date(substr(names(r), 1, 10))
  keep <- d > obs_range[1] & d < obs_range[2]
  if (!any(keep)) return(NULL)
  r <- r[[which(keep)]]; d <- d[keep]

  # QC already applied by script 3c

  # class and landform on the same grid
  ref <- r[[1]]
  cls <- function(shp, field) {
    v <- vect(shp); v <- v[v$Site == site, ]
    if (!nrow(v)) return(rep(NA_character_, ncell(ref)))
    g <- terra::rasterize(project(v, crs(ref)), ref, field = field, fun = "max")
    lv <- levels(g)[[1]]
    if (is.null(lv)) as.character(values(g)[, 1])
    else as.character(factor(values(g)[, 1], levels = lv[[1]], labels = lv[[2]]))
  }
  mvg <- cls(MVG_shp, "MVG_Nam")
  lf  <- cls(rip_shp, "Landform")

  v  <- values(r)
  ok <- which(!is.na(mvg) & rowSums(is.finite(v)) >= min_obs)
  if (!length(ok)) return(NULL)
  xy <- crds(ref, na.rm = FALSE)[ok, , drop = FALSE]

  out <- data.frame(
    Site = site, source = src_name, pixel = rep(ok, each = length(d)),
    x = rep(xy[, 1], each = length(d)), y = rep(xy[, 2], each = length(d)),
    MVG_Nam = rep(mvg[ok], each = length(d)), Landform = rep(lf[ok], each = length(d)),
    Date = rep(d, times = length(ok)),
    value = as.numeric(t(v[ok, , drop = FALSE])), stringsAsFactors = FALSE)
  out <- out[is.finite(out$value), ]

  # index -> LAI, after masking implausible index values
  if (src$kind == "vi") {
    rng <- vi_valid[[sub("^raw", "", src$stem)]]
    out <- out[out$value >= rng[1] & out$value <= rng[2], ]
    nd  <- setNames(data.frame(out$value), paste0(src$stem, "_closest"))
    out$LAI <- as.numeric(predict(lin_mods[[src$model]], newdata = nd))
  } else if (src$kind == "cal") {
    out$LAI <- as.numeric(predict(leaf_cal,
                 newdata = data.frame(leaf_10m_closest = out$value)))
  } else out$LAI <- out$value

  out$LAI[out$LAI < 0] <- 0 # LAI below zero is clipped
  out[, c("Site", "source", "pixel", "x", "y", "MVG_Nam", "Landform", "Date", "LAI")]
}

## -------------------------------------------------- smoothing and phases ----
# Savitzky-Golay filter of one pixel's LAI. The observations are linearly filled
# to a grid sg_step days apart (gaps over max_fill days stay empty) and each grid
# point takes the order-sg_order polynomial fitted over the sg_days around it; at
# the ends the nearest full window is fitted. Negative LAI is set to zero. Used by
# the turning points and the test-site figure, so the smooth is defined here only.
sg_n <- sg_days / sg_step + 1
sg_F <- local({ k <- (sg_n - 1) / 2; t <- -k:k; X <- outer(t, 0:sg_order, `^`); X %*% solve(crossprod(X), t(X)) })
sg_smooth <- function(day, y) {
  m <- length(day); g <- seq(day[1], day[m], by = sg_step); N <- length(g)
  if (N < sg_n) return(NULL)
  yg <- approx(day, y, xout = g, rule = 2, ties = mean)$y
  k <- (sg_n - 1) / 2
  out <- as.numeric(stats::filter(yg, sg_F[k + 1, ], sides = 2))
  out[1:k] <- sg_F[1:k, , drop = FALSE] %*% yg[1:sg_n]
  out[(N - k + 1):N] <- sg_F[(k + 2):sg_n, , drop = FALSE] %*% yg[(N - sg_n + 1):N]
  j <- pmin(findInterval(g, day, rightmost.closed = TRUE), m - 1L)
  out[(day[j + 1L] - day[j]) > max_fill & g > day[j] & g < day[j + 1L]] <- NA
  data.frame(day = g, LAI = pmax(out, 0))
}

# Smooth, then strict local maxima and minima; the third condition marks the end
# of a flat-at-zero stretch. full = TRUE also returns the smooth for the
# protocol figure below.
turning_points <- function(d, full = FALSE) {
  if (nrow(d) < min_obs) return(NULL)
  d  <- d[order(d$Date), ]
  sm <- sg_smooth(as.numeric(d$Date), d$LAI)
  if (is.null(sm)) return(NULL)
  grid <- as.Date(sm$day, origin = "1970-01-01"); lai <- sm$LAI
  n <- length(lai); i <- 2:(n - 1)
  peak   <- lai[i] > lai[i - 1] & lai[i] > lai[i + 1]
  trough <- lai[i] < lai[i - 1] & lai[i] < lai[i + 1]
  flat   <- lai[i] == 0 & lai[i - 1] == 0 & lai[i] < lai[i + 1]
  peak[is.na(peak)] <- FALSE; trough[is.na(trough)] <- FALSE; flat[is.na(flat)] <- FALSE
  turn <- rep(NA_integer_, n)
  turn[i[peak]] <- 1L; turn[i[trough | flat]] <- 0L
  k <- which(!is.na(turn))
  if (length(k) < 2) return(NULL)

  # A gap left empty can put two peaks, or two troughs, in a row: keep the higher
  # peak and the lower trough.
  repeat {
    s <- which(turn[k][-1] == turn[k][-length(k)])
    if (!length(s)) break
    a <- k[s[1]]; b <- k[s[1] + 1L]
    k <- k[k != if (turn[a] == 1L) { if (lai[a] >= lai[b]) b else a } else { if (lai[a] <= lai[b]) b else a }]
  }
  if (length(k) < 2) return(NULL)

  # Drop swings smaller than min_amp. Turning points alternate peak/trough by
  # construction, so removing an adjacent PAIR preserves that alternation -
  # repeatedly take out the smallest swing until none is below the threshold.
  if (min_amp > 0) {
    repeat {
      if (length(k) < 2) break
      a <- abs(diff(lai[k]))
      i <- which.min(a)
      if (a[i] >= min_amp) break
      k <- k[-c(i, i + 1L)]
    }
    if (length(k) < 2) return(NULL)
  }

  tp <- data.frame(d[1, c("Site", "source", "pixel", "x", "y", "MVG_Nam", "Landform")],
                   Date = grid[k], LAI = lai[k], turn = turn[k],
                   LAI_glob_mean = mean(lai, na.rm = TRUE), LAI_glob_sd = sd(lai, na.rm = TRUE),
                   row.names = NULL, stringsAsFactors = FALSE)
  if (!full) return(tp)
  list(tp = tp, smooth = data.frame(Date = grid, LAI = lai))
}

# Climate over a window, by date arithmetic on a complete daily series. Returns
# NA rather than a wrong answer when the window falls outside the record. An event
# is a run of consecutive days of at least event_min_day mm; LRE is the biggest
# event total (0 if there is none), RF all the rain, n_wet the days with any rain.
# RE2 ... RE<n_large> are the next largest events (0 if there are fewer) and
# MRE, MRE2, ... the rain left once the 1, 2, ... n_large largest events are out.
re_extra <- if (n_large > 1) paste0("RE", 2:n_large)
mre_cols <- c("MRE", if (n_large > 1) paste0("MRE", 2:n_large))
clim_window <- function(cl, from, to) {
  i1 <- as.integer(from - cl$Date[1]) + 1L
  i2 <- as.integer(to   - cl$Date[1]) + 1L
  if (is.na(i1) || is.na(i2) || i1 < 1L || i2 > nrow(cl) || i2 < i1)
    return(setNames(rep(NA_real_, 5 + length(re_extra) + length(mre_cols)),
                    c("RF", "LRE", "Tmax", "Tmax_sd", "n_wet", re_extra, mre_cols)))
  rf <- cl$Rainfall[i1:i2]; tx <- cl$Tmax_dC[i1:i2]
  wet <- which(rf > 0)
  ev  <- which(rf > 0 & rf >= event_min_day)
  et  <- if (length(ev)) as.numeric(tapply(rf[ev], cumsum(c(1, diff(ev) != 1)), sum)) else numeric(0)
  re  <- c(sort(et, decreasing = TRUE), rep(0, n_large))[seq_len(n_large)]   # largest first
  c(RF = sum(rf), LRE = re[1], Tmax = mean(tx), Tmax_sd = sd(tx), n_wet = length(wet),
    setNames(re[-1], re_extra), setNames(pmax(sum(rf) - cumsum(re), 0), mre_cols))
}

# One pixel's turning points -> one row per phase, with the response and the
# climate of the phase and of each lagged window before it.
phase_stats <- function(tp, cl) {
  tp <- tp[order(tp$Date), ]
  if (nrow(tp) < 2) return(NULL)
  do.call(rbind, lapply(2:nrow(tp), function(p) {
    prev <- tp$Date[p - 1]; cur <- tp$Date[p]
    dL <- tp$LAI[p] - tp$LAI[p - 1]
    w  <- clim_window(cl, prev, cur)
    out <- data.frame(
      tp[p, c("Site", "source", "pixel", "x", "y", "MVG_Nam", "Landform")],
      Date = cur, turn = tp$turn[p],
      phase_length = as.numeric(cur - prev),
      LAI = tp$LAI[p],
      LAI_prev = tp$LAI[p - 1],
      LAI_change = dL,
      LAI_prop_change = dL / max(tp$LAI[p - 1], 0.005),
      LAI_norm_change = dL / tp$LAI_glob_sd[p],
      RF_phase = w[["RF"]],
      LRE_phase = w[["LRE"]],
      Tmax_mean_phase = w[["Tmax"]],
      Tmax_sd_phase = w[["Tmax_sd"]],
      n_RF_days_phase = w[["n_wet"]],
      length_prev_phase = if (p > 2) as.numeric(prev - tp$Date[p - 2]) else NA,
      row.names = NULL, stringsAsFactors = FALSE)
    if (p > 2) {
      wp <- clim_window(cl, tp$Date[p - 2], prev)
      out$Tmax_mean_prev_phase <- wp[["Tmax"]]
      out$Tmax_sd_prev_phase   <- wp[["Tmax_sd"]]
    } else { out$Tmax_mean_prev_phase <- NA; out$Tmax_sd_prev_phase <- NA }
    for (L in lags) {
      wl <- clim_window(cl, prev - L, cur - L)
      out[[paste0("RF_", L, "_Dlag")]]      <- wl[["RF"]]
      out[[paste0("LRE_", L, "_Dlag")]]     <- wl[["LRE"]]
      for (nm in c(re_extra, mre_cols)) out[[paste0(nm, "_", L, "_Dlag")]] <- wl[[nm]]
      out[[paste0("Tmax_", L, "_Dlag")]]    <- wl[["Tmax"]]
      out[[paste0("Tmax_sd_", L, "_Dlag")]] <- wl[["Tmax_sd"]]
    }
    out
  }))
}

## -------------------------------------------------------------------- run ---
clim_DF <- readRDS(file.path(end_path, "site_climate_AGCD.rds"))

run_site <- function(site) {
  cl <- clim_DF[[site]]
  if (is.null(cl)) return(NULL)
  cl <- cl[order(cl$Date), ]
  if (any(diff(cl$Date) != 1)) stop("climate for ", site, " is not a complete daily series")
  out <- lapply(names(sources), function(sn) {
    d <- extract_site(site, sn)
    if (is.null(d)) return(NULL)
    tp <- do.call(rbind, lapply(split(d, d$pixel), turning_points))
    if (is.null(tp)) return(NULL)
    do.call(rbind, lapply(split(tp, tp$pixel), phase_stats, cl = cl))
  })
  out <- do.call(rbind, Filter(Negate(is.null), out))
  if (!is.null(out)) message("  ", site, ": ", nrow(out), " phases across ",
                             length(unique(out$pixel)), " pixels")
  out
}

plan(multisession, workers = pmin(availableCores(omit = 2), length(site_list)))
phase_df <- do.call(rbind, Filter(Negate(is.null),
  future_lapply(site_list, run_site, future.seed = TRUE,
                future.packages = "terra")))
plan(sequential)

saveRDS(phase_df, file.path(end_path, "LAI_phase_stats_by_source.rds"))
message(nrow(phase_df), " phases, ", length(unique(phase_df$Site)), " sites, ",
        length(unique(phase_df$source)), " sources -> ",
        file.path(end_path, "LAI_phase_stats_by_source.rds"))
print(table(phase_df$source))

################################################################################
# example of the protocol for a single pixel
################################################################################
# One pixel end to end: the filtered index, the LAI modelled from it, the smooth,
# the turning points, the phases they bound, and the rainfall over each phase.
# The index is in the top panel. Conversion happens per observation and the
# smoother runs on LAI, so the LAI points are already in LAI units.
library(ggplot2)
library(patchwork)

protocol_figure <- function(site = "MCK", src_name = vi_choice$stem, pixel = NULL,
                            lag_days = 14,
                            file = file.path(plot_path, "example_pixel_LAI_timeseries.tif")) {
  d <- extract_site(site, src_name)
  if (is.null(d)) stop("no data for ", site, " / ", src_name)
  if (is.null(pixel)) {                       # the best-observed pixel
    n <- table(d$pixel)
    pixel <- as.integer(names(n)[which.max(n)])
  }
  d <- d[d$pixel == pixel, ]
  z <- turning_points(d, full = TRUE)
  if (is.null(z)) stop("no turning points for pixel ", pixel)
  tp <- z$tp; sm <- z$smooth
  d  <- d[d$Date >= min(sm$Date), ]           # observations within the drawn window

  # the index itself at those dates, before conversion to LAI
  src <- sources[[src_name]]
  idx <- NULL
  if (src$kind == "vi") {
    r   <- rast(file.path(src$dir, site, paste0(src$stem, "-ts.tif")))
    idx <- data.frame(Date = as.Date(substr(names(r), 1, 10)), value = as.numeric(unlist(r[pixel])))
    idx <- idx[idx$Date %in% d$Date, ]
  }
  lab_obs <- sprintf("LAI modelled from %s (NBART corrected)", vi_choice$stem)

  # The bracket beneath each growth phase marks the window the lagged climate
  # predictors are summed over: (phase start - lag) to (phase end - lag). It is
  # drawn at the height of the trough that opens the phase, so it reads as
  # "this rain, shifted back by lag days, drove that rise".
  st <- which(tp$turn == 0 & seq_len(nrow(tp)) < nrow(tp))
  ph <- data.frame(x    = tp$Date[st] - lag_days,
                   xend = tp$Date[st + 1] - lag_days,
                   y    = tp$LAI[st])
  tick <- diff(range(sm$LAI, na.rm = TRUE)) * 0.02

  cl <- clim_DF[[site]]
  cl <- cl[cl$Date >= min(sm$Date) & cl$Date <= max(sm$Date) & cl$Rainfall > 0.1, ]
  sf <- 0.4 * max(sm$LAI, na.rm = TRUE) / max(cl$Rainfall, na.rm = TRUE)

  ttl <- sprintf("%s  pixel %d  (%s, %s)", site, pixel, src_name, d$MVG_Nam[1])
  xl  <- range(sm$Date)
  fs  <- 1.2 * txt_size
  txt <- theme(panel.grid = element_blank(), legend.title = element_blank(),
               plot.title = element_text(size = fs, face = "plain"),
               axis.text = element_text(size = fs, colour = "black"),
               axis.title = element_text(size = fs),
               legend.text = element_text(size = fs))

  p <- ggplot() +
    geom_segment(data = cl, aes(x = Date, xend = Date, y = 0, yend = Rainfall * sf,
                                colour = "Rainfall"), linewidth = 0.4) +
    geom_point(data = d, aes(x = Date, y = LAI, colour = lab_obs), size = 0.9,
               alpha = 0.55) +
    geom_line(data = sm, aes(x = Date, y = LAI, colour = "smoothed LAI"),
              linewidth = 0.8) +
    # climate window for each growth phase, lagged back by lag_days
    geom_segment(data = ph, aes(x = x, xend = xend, y = y, yend = y),
                 colour = "grey35", linewidth = 0.4) +
    geom_segment(data = ph, aes(x = x, xend = x, y = y - tick, yend = y + tick),
                 colour = "grey35", linewidth = 0.4) +
    geom_segment(data = ph, aes(x = xend, xend = xend, y = y - tick, yend = y + tick),
                 colour = "grey35", linewidth = 0.4) +
    geom_point(data = tp, aes(x = Date, y = LAI,
                              colour = ifelse(turn == 1, "peak", "trough")), size = 2.4) +
    scale_colour_manual(values = c(setNames("grey30", lab_obs),
                                   `smoothed LAI` = "darkgreen",
                                   trough = "blue3", peak = "red3",
                                   Rainfall = "lightblue3"),
                        breaks = c(lab_obs, "smoothed LAI", "peak", "trough",
                                   "Rainfall")) +
    guides(colour = guide_legend(override.aes = list(linewidth = 0.8))) +   # rainfall key as wide as the smooth
    scale_x_date(name = "", limits = xl, date_breaks = "1 year", date_labels = "%Y") +
    scale_y_continuous(name = expression(LAI~(m^2~m^-2)),
                       sec.axis = sec_axis(~ . / sf, name = "daily rainfall (mm)")) +
    labs(colour = "") +
    theme_bw() + txt +
    theme(legend.position = "bottom")
  if (is.null(idx)) {
    fig <- p + labs(title = ttl)
  } else {
    p_idx <- ggplot(idx, aes(x = Date, y = value)) +
      geom_point(colour = "darkorange2", size = 0.9, alpha = 0.7) +
      scale_x_date(name = "", limits = xl, date_breaks = "1 year", date_labels = "%Y") +
      scale_y_continuous(name = sprintf("%s (NBART corrected)", src$stem)) +
      labs(title = ttl) +
      theme_bw() + txt +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    fig <- p_idx / p + plot_layout(heights = c(1, 2.2))
  }
  ggsave(file, fig, width = 25, height = if (is.null(idx)) 13 else 18, units = "cm", dpi = 300)
  message("protocol figure: ", site, " pixel ", pixel, " -> ", file)
  fig
}

protocol_figure()

################################################################################
# LAI timeseries at the test sites
################################################################################
# Mean LAI of one NVIS class at each test site from four datasets: the index
# after the DEA filtering of script 3b, the same index from unfiltered Sen2Cor,
# SL2P 10 m, and SL2P 10 m calibrated. Each goes to LAI through its own
# regression (the unfiltered index through its own fit, as the raw panels of
# script 5a do). The points, the smoother used for the turning points applied to
# the class mean, and a ribbon of one SD across the class's pixels on each date.
# The triangles are the field campaigns.
min_valid_frac <- 0.80    # same name and value as scripts 2a, 2f and 3b

raw_src    <- paste0("raw", vi_choice$stem)
ts_sources <- c(sources, setNames(list(list(kind = "vi", stem = raw_src, dir = raw_path,
                                            model = paste0("raw", vi_choice$key))), raw_src))
ts_names   <- c(raw_src, vi_choice$stem, "SL2P 10 m", "SL2P 10 m calibrated")
ts_cols    <- setNames(c("grey50", "darkgreen", "darkorange2", "purple3"),
                       c(sprintf("%s (Sen2Cor)", vi_choice$stem),
                         sprintf("%s (NBART)", vi_choice$stem),
                         "SL2P", "SL2P calibrated"))

# Class mean per date. A date is dropped when fewer than min_valid_frac of the
# best-observed date's pixels are clear - the mean of a small remnant is what spikes.
class_series <- function(site, src_name, srcs = ts_sources, mvg = ts_mvg[[site]]) {
  d <- extract_site(site, src_name, srcs)
  if (is.null(d)) return(NULL)
  d <- d[d$MVG_Nam == mvg, ]
  if (!nrow(d)) return(NULL)
  n <- tapply(d$LAI, d$Date, length)
  m <- tapply(d$LAI, d$Date, mean)
  s <- tapply(d$LAI, d$Date, sd)
  ok <- n >= min_valid_frac * max(n)
  data.frame(Site = site, source = src_name, Date = as.Date(names(m)[ok]),
             LAI = as.numeric(m[ok]), sd = as.numeric(s[ok]), row.names = NULL)
}

jobs <- expand.grid(site = test_sites, src = ts_names, stringsAsFactors = FALSE)
plan(multisession, workers = pmin(availableCores(omit = 2), nrow(jobs)))
ts_pts <- do.call(rbind, future_lapply(seq_len(nrow(jobs)),
  function(i) class_series(jobs$site[i], jobs$src[i]), future.seed = TRUE,
  future.packages = "terra"))
plan(sequential)

# the ribbon is broken where the series is, i.e. across gaps the smoother leaves empty
ts_pts <- ts_pts[order(ts_pts$Site, ts_pts$source, ts_pts$Date), ]
ts_pts$seg <- ave(as.numeric(ts_pts$Date), ts_pts$Site, ts_pts$source,
                  FUN = function(d) cumsum(c(0, diff(d) > max_fill)))

# the smoother of the turning points, applied to each class-mean series
class_smooth <- function(p) {
  out <- do.call(rbind, lapply(split(p, list(p$Site, p$source), drop = TRUE), function(s) {
    z <- sg_smooth(as.numeric(s$Date), s$LAI)
    if (!is.null(z)) data.frame(Site = s$Site[1], source = s$source[1],
                                Date = as.Date(z$day, origin = "1970-01-01"), LAI = z$LAI)
  }))
  rownames(out) <- NULL
  out
}
ts_line <- class_smooth(ts_pts)
print(table(ts_pts$Site, ts_pts$source))

ref_pts <- paired[paired$Site %in% test_sites & paired$MVG == unname(ts_mvg[paired$Site]), ]
ref_pts <- data.frame(Site = ref_pts$Site, Date = as.Date(ref_pts$Date), LAI = ref_pts$LAI_LX_total)

rain <- do.call(rbind, clim_DF[test_sites])
rain <- rain[!is.na(rain$Rainfall) & rain$Rainfall > 0.1 &
             rain$Date >= min(ts_pts$Date) & rain$Date <= max(ts_pts$Date), ]

panel_lab <- function(s) factor(paste0(s, " - ", ts_mvg[s]),
                                levels = paste0(test_sites, " - ", ts_mvg[test_sites]))
dataset_lab <- function(s) factor(s, levels = ts_names, labels = names(ts_cols))
ts_pts$dataset  <- dataset_lab(ts_pts$source);  ts_pts$panel  <- panel_lab(ts_pts$Site)
ts_line$dataset <- dataset_lab(ts_line$source); ts_line$panel <- panel_lab(ts_line$Site)
ref_pts$panel   <- panel_lab(ref_pts$Site); rain$panel    <- panel_lab(rain$Site)

# One LAI panel per site, each on a y axis of its own so a small spread (YU) can
# be seen, with the site's daily rainfall in a panel joined beneath it (the same
# rain axis at every site)
ts_figure <- function(file) {
  xlim <- range(ts_pts$Date)
  fs   <- 1.2 * txt_size
  th   <- theme_bw() +
    theme(legend.text = element_text(size = fs), axis.text = element_text(size = fs, colour = "black"),
          axis.title = element_text(size = fs), strip.text = element_text(size = fs),
          panel.grid = element_blank())
  no_x <- theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  x_sc <- scale_x_date(limits = xlim, date_breaks = "1 year", date_labels = "%Y", expand = expansion(0.01),
                       minor_breaks = seq(as.Date(cut(xlim[1], "year")), xlim[2], by = "6 months"),
                       guide = guide_axis(minor.ticks = TRUE))
  site_plots <- function(s, last) {
    pts   <- ts_pts[ts_pts$Site == s, ]
    y_max <- max(c(pts$LAI + pts$sd, ref_pts$LAI[ref_pts$Site == s]), na.rm = TRUE)
    p_lai <- ggplot(mapping = aes(x = Date)) +
      geom_ribbon(data = pts, alpha = 0.2,
                  aes(ymin = pmax(LAI - sd, 0), ymax = LAI + sd, fill = dataset, group = interaction(dataset, seg))) +
      geom_point(data = pts, aes(y = LAI, colour = dataset), size = 0.9, alpha = 0.3) +
      geom_line(data = ts_line[ts_line$Site == s, ], aes(y = LAI, colour = dataset),
                linewidth = 0.7, na.rm = TRUE) +
      geom_point(data = ref_pts[ref_pts$Site == s, ], aes(y = LAI, shape = "reference total LAI"),
                 size = 3, stroke = 0.9, colour = "black") +
      facet_wrap(~panel) +
      scale_colour_manual(values = ts_cols) +
      scale_fill_manual(values = ts_cols) +
      scale_shape_manual(values = c(`reference total LAI` = 2)) +
      x_sc +
      scale_y_continuous(name = expression(LAI~(m^2~m^-2)), limits = c(0, y_max)) +
      labs(colour = "", shape = "") +
      guides(colour = guide_legend(nrow = 2, order = 1,
                                   override.aes = list(alpha = 1, size = 2, linewidth = 1.2)),
             shape  = guide_legend(nrow = 1, order = 2), fill = "none") +
      th + no_x + theme(axis.ticks.x = element_blank(), plot.margin = margin(5.5, 5.5, 0, 5.5))
    p_rain <- ggplot(rain[rain$Site == s, ], aes(Date, ymin = 0, ymax = Rainfall)) +
      geom_linerange(colour = "steelblue", linewidth = 0.4) +
      x_sc +
      scale_y_continuous(name = "Rainfall\n(mm)", limits = c(0, 200), breaks = c(0, 100, 200), minor_breaks = c(50, 150),
                         expand = expansion(0), position = "right", guide = guide_axis(minor.ticks = TRUE)) +
      labs(x = NULL) + th + theme(plot.margin = margin(0, 5.5, 5.5, 5.5))
    if (!last) p_rain <- p_rain + no_x
    list(p_lai, p_rain)
  }
  ps  <- unlist(lapply(seq_along(test_sites), function(i) site_plots(test_sites[i], i == length(test_sites))),
                recursive = FALSE)
  fig <- wrap_plots(ps, ncol = 1, guides = "collect", heights = rep(c(1, 0.4), length(test_sites))) &
    theme(legend.position = "bottom")
  ggsave(file, fig, width = 25, height = 32, units = "cm", dpi = 300)
  message("timeseries figure -> ", file)
  invisible(fig)
}

ts_figure(file.path(plot_path, "LAI_timeseries_test_sites_sd.tif"))

################################################################################
# smoothed class mean against the reference LAI
################################################################################
# Every index model (both indices, filtered NBART and unfiltered Sen2Cor) and
# SL2P with and without calibration: the smoothed class mean of the figure above
# on each field visit, against the reference LAI. The other index is extracted
# here; the rest are the figure's series.
oth     <- setdiff(c("NDVI", "EVI"), vi_choice$stem)
acc_src <- c(ts_sources, setNames(list(
  list(kind = "vi", stem = oth, dir = S2_path, model = paste0(oth, "_", vi_choice$form)),
  list(kind = "vi", stem = paste0("raw", oth), dir = raw_path, model = paste0("raw", oth, "_", vi_choice$form))),
  c(oth, paste0("raw", oth))))
oth_jobs <- expand.grid(site = test_sites, src = c(oth, paste0("raw", oth)), stringsAsFactors = FALSE)
plan(multisession, workers = pmin(availableCores(omit = 2), nrow(oth_jobs)))
oth_pts <- do.call(rbind, future_lapply(seq_len(nrow(oth_jobs)),
  function(i) class_series(oth_jobs$site[i], oth_jobs$src[i], acc_src), future.seed = TRUE,
  future.packages = "terra"))
plan(sequential)

acc_lab <- c(NDVI = "NDVI (NBART corrected)", EVI = "EVI (NBART corrected)",
             rawNDVI = "NDVI (Sen2Cor L2A)", rawEVI = "EVI (Sen2Cor L2A)",
             `SL2P 10 m` = "SL2P 10 m (uncalibrated)", `SL2P 10 m calibrated` = "SL2P 10 m (calibrated)")
oth_line <- class_smooth(oth_pts)
sm_all   <- rbind(ts_line[names(oth_line)], oth_line)
acc <- do.call(rbind, lapply(split(sm_all, list(sm_all$Site, sm_all$source), drop = TRUE), function(s) {
  r   <- ref_pts[ref_pts$Site == s$Site[1], ]
  est <- approx(as.numeric(s$Date), s$LAI, xout = as.numeric(r$Date), na.rm = FALSE)$y   # NA in a gap
  data.frame(Model = acc_lab[[s$source[1]]], Site = r$Site, Date = r$Date, reference = r$LAI,
             smoothed = est, diff = est - r$LAI)
}))
acc$Model <- factor(acc$Model, levels = acc_lab)
acc <- acc[order(acc$Model, acc$Site, acc$Date), ]
message("smoothed class mean against the reference LAI at the test sites (NRMSE = RMSE / mean reference)")
print(acc, digits = 2, row.names = FALSE)
print(do.call(rbind, lapply(split(acc, acc$Model), function(a) { a <- a[is.finite(a$diff), ]
  data.frame(Model = a$Model[1], n = nrow(a), max_abs_diff = round(max(abs(a$diff)), 2),
             RMSE = round(sqrt(mean(a$diff^2)), 2), NRMSE_pct = round(100 * sqrt(mean(a$diff^2)) / mean(a$reference))) })),
  row.names = FALSE)

################################################################################
# change between the two campaigns
################################################################################
# Every test-site vegetation type with a reference record in both the late dry (October
# 2021) and early wet (March 2022) campaigns: the change in reference LAI and in
# the smoothed class mean of each figure series between the two visits.
chg_cls <- aggregate(Date ~ Site + MVG, unique(paired[paired$Site %in% test_sites, c("Site", "MVG", "Date")]), length)
chg_jobs <- merge(chg_cls[chg_cls$Date == 2, c("Site", "MVG")], data.frame(src = ts_names), by = NULL)
plan(multisession, workers = pmin(availableCores(omit = 2), nrow(chg_jobs)))
chg_pts <- do.call(rbind, future_lapply(seq_len(nrow(chg_jobs)), function(i) {
  p <- class_series(chg_jobs$Site[i], chg_jobs$src[i], mvg = chg_jobs$MVG[i])
  if (!is.null(p)) cbind(p, MVG = chg_jobs$MVG[i]) }, future.seed = TRUE, future.packages = "terra"))
plan(sequential)

chg <- do.call(rbind, lapply(split(chg_pts, list(chg_pts$Site, chg_pts$MVG, chg_pts$source), drop = TRUE), function(s) {
  r   <- paired[paired$Site == s$Site[1] & paired$MVG == s$MVG[1], ]
  r   <- r[order(r$Date), ]
  z   <- sg_smooth(as.numeric(s$Date), s$LAI)
  est <- if (is.null(z)) c(NA, NA) else approx(z$day, z$LAI, xout = as.numeric(as.Date(r$Date)), na.rm = FALSE)$y
  data.frame(Site = s$Site[1], MVG = s$MVG[1], reference = diff(r$LAI_LX_total),
             Model = acc_lab[[s$source[1]]], smoothed = diff(est))
}))
chg$Model <- factor(chg$Model, levels = acc_lab)
chg <- chg[order(chg$Model, chg$Site, chg$MVG), ]
chg_w <- reshape(chg, idvar = c("Site", "MVG", "reference"), timevar = "Model", direction = "wide")
names(chg_w) <- sub("^smoothed\\.", "", names(chg_w))
message("change in LAI from the late dry to the early wet campaign, test-site vegetation types visited in both (reference ",
        paste(sprintf("%+.2f", range(chg_w$reference)), collapse = " to "), ")")
print(cbind(chg_w[1:2], round(chg_w[-(1:2)], 2)), row.names = FALSE)
print(do.call(rbind, lapply(split(chg, chg$Model, drop = TRUE), function(a) { a <- a[is.finite(a$smoothed), ]
  data.frame(Model = a$Model[1], n = nrow(a), increased = sum(a$smoothed > 0),
             same_direction = sum(sign(a$smoothed) == sign(a$reference)),
             same_direction_pct = round(100 * mean(sign(a$smoothed) == sign(a$reference)))) })),
  row.names = FALSE)

