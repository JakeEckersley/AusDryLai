################################################################################
#### 5e - rainfall assumptions                                              ####
################################################################################
# Checks behind the rainfall terms of scripts 5b and 5c:
#   1. which lag lines the rain window up best with growth: Pearson r between the
#      standardised LAI change of growth phases (dLAI over the pixel's SD of
#      smoothed LAI, 5b's LAI_norm_change) and the rain of the window lagged by
#      0-56 d (t1 - lag to t2 - lag), all of it (RF) and its largest event (RE1),
#      for each source of 5c. 95% intervals from n_boot resamples of sites.
#   2. AGCD rainfall against the DPIRD gauges. MILL (Millstream) and KARI
#      (Karijini), script 2h, are DBCA gauges, not BOM stations AGCD is built
#      from, so each is an independent check on the AGCD series that script 5b
#      uses. AGCD is read at the 0.05 deg cell that holds the gauge, for the days
#      the gauge recorded at least min_obs % of its minutes, from its first such
#      day to its last or the end of the AGCD grids (KARI logged every 15 minutes
#      until mid 2017, so its earlier days fall short of min_obs). Both are the 24
#      hours to 9 am. One figure per gauge: daily rainfall of the two (gauge up,
#      AGCD down, on a signed log axis, days without a usable gauge record
#      shaded), their daily scatter, and their monthly totals.

library(terra)
library(ggplot2)
library(patchwork)
library(scales)

fig_path <- file.path("plots", "RainfallAssumptions")
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

#############################################################################
# 1. rain over lagged windows against growth
#############################################################################
end_path    <- "outputs"
mod_sources <- c(readRDS(file.path(end_path, "Optimum_closest_VI_choice.rds"))$stem, "SL2P 10 m calibrated")
src_short   <- setNames(c(mod_sources[1], "SL2Pcal"), mod_sources)
lags   <- seq(0, 56, by = 7)       # as stored by 5b
n_boot <- 1000

# 5c's growth phases
ph <- readRDS(file.path(end_path, "LAI_phase_stats_by_source.rds"))
ph <- ph[ph$source %in% mod_sources & ph$turn %in% 1 & !is.na(ph$Landform) & ph$Landform != "Bare channel", ]

# r from per-site sums, so each resample of sites is a weighted sum of them
pearson <- function(t) (t[, 1] * t[, 6] - t[, 2] * t[, 3]) /
  sqrt((t[, 1] * t[, 4] - t[, 2]^2) * (t[, 1] * t[, 5] - t[, 3]^2))
terms <- c(RF = "RF", RE1 = "LRE")  # 5b's column prefixes
set.seed(123)
lag_cor <- do.call(rbind, lapply(mod_sources, function(s) {
  d <- ph[ph$source == s, ]; site <- factor(d$Site)
  W <- t(rmultinom(n_boot, nlevels(site), rep(1, nlevels(site))))   # times each site is drawn
  do.call(rbind, lapply(names(terms), function(tm) do.call(rbind, lapply(lags, function(L) {
    x  <- d$LAI_norm_change; y <- d[[paste0(terms[[tm]], "_", L, "_Dlag")]]
    ok <- is.finite(x) & is.finite(y)
    S  <- rowsum(cbind(1, x, y, x^2, y^2, x * y)[ok, ], site[ok])
    rb <- pearson(W[, match(rownames(S), levels(site)), drop = FALSE] %*% S)
    data.frame(source = src_short[[s]], term = tm, lag = L, n = sum(ok), r = pearson(matrix(colSums(S), 1)),
               lower = unname(quantile(rb, 0.025)), upper = unname(quantile(rb, 0.975)))
  }))))
}))
rm(ph); invisible(gc())
print(lag_cor, digits = 3, row.names = FALSE)
print(do.call(rbind, lapply(split(lag_cor, list(lag_cor$source, lag_cor$term)), function(a) a[which.max(a$r), ])),
      digits = 3, row.names = FALSE)

# one panel per source and term, in the style of the earlier version of this figure
lag_panel <- function(src, tm) {
  a <- lag_cor[lag_cor$source == src & lag_cor$term == tm, ]
  x <- if (tm == "RF") "RF[t[1] - delay %->% t[2] - delay]" else "RE[list(1, t[1] - delay %->% t[2] - delay)]"
  ggplot(a, aes(x = lag, y = r)) +
    geom_point(size = 2.5, color = "steelblue") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, color = "steelblue") +
    theme_bw() +
    labs(title = parse(text = sprintf("Delta * LAI[t[1] %%->%% t[2]] ~ '~' ~ %s ~ '(%s)'", x, src)),
         x = expression(paste("Growth onset delay (", italic(n), " days)")),
         y = "Pearson\u2019s R (\u00b1 95 % CI)") +
    scale_x_continuous(breaks = lags) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text = element_text(size = 15, colour = 'black'), axis.title = element_text(size = 15),
          axis.line = element_line(colour = "black"))
}
lag_fig <- wrap_plots(lapply(seq_len(2 * length(src_short)), function(i)
  lag_panel(src_short[[ceiling(i / 2)]], names(terms)[2 - i %% 2])), ncol = 2)
ggsave(file.path(fig_path, "Supp_rainfall_lag_correlations.tif"), lag_fig, width = 25, height = 28, units = "cm", dpi = 300)

#############################################################################
# 2. AGCD rainfall against the DPIRD gauges
#############################################################################
gauges    <- c("MILL", "KARI")
grid_dir  <- "data/climGrids"
dpird_dir <- "data/bom_data/dpird"
min_obs   <- 95            # % of a day's minutes recorded
min_month <- 0.90          # share of a month's days needed for its total
txt_size  <- 11
col_g <- "steelblue4"; col_a <- "darkorange2"

grids <- Sys.glob(file.path(grid_dir, "agcd_precip_total_*_pilbara.nc"))
r1 <- rast(grids[1]); crs(r1) <- "EPSG:4326"

base_theme <- theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size = txt_size, colour = "black"),
        axis.title = element_text(size = txt_size), plot.subtitle = element_text(size = txt_size - 1),
        legend.text = element_text(size = txt_size), legend.title = element_blank())

stats <- function(x, y) c(n = length(x), r = cor(x, y), rmse = sqrt(mean((y - x)^2)), bias = mean(y - x),
                          total_ratio = sum(y) / sum(x))

# n, bias and r in the top left of a scatter, one line per annotate() as in 5a's
add_stats <- function(g, s, n_unit, unit, dp, lim, tr = identity_trans()) {
  ann <- c(sprintf("n = %s %s", format(s[["n"]], big.mark = ","), n_unit),
           sprintf(paste0("bias = %+.", dp, "f %s"), round(s[["bias"]], dp) + 0, unit),
           sprintf("italic(r) == '%.2f'", s[["r"]]))
  y <- tr$inverse(tr$transform(lim) * (1 - 0.07 * (seq_along(ann) - 1)))
  for (k in seq_along(ann))
    g <- g + annotate("text", x = 0, y = y[k], label = ann[k], parse = k == 3, hjust = 0, vjust = 1, size = txt_size / 3.2)
  g
}

# axis limit 5 % past the data, raised where a point would sit under those lines
ax_lim <- function(x, y, tr = identity_trans()) {
  lim <- max(x, y) * 1.05
  top <- tr$transform(y[tr$transform(x) < 0.5 * tr$transform(lim)])
  max(lim, tr$inverse(max(top, 0) / 0.75))
}
plog <- pseudo_log_trans(sigma = 1, base = 10)

compare <- function(code) {
  st <- read.csv(file.path(dpird_dir, paste0(code, "_station.csv")))
  g  <- read.csv(file.path(dpird_dir, paste0(code, "_daily_rainfall.csv")))
  g$date <- as.Date(g$date)
  g  <- data.frame(date = g$date, gauge = ifelse(g$obs_pct >= min_obs, g$rainfall_mm, NA))

  # AGCD at the gauge: the cell that holds it
  xy <- cbind(st$lon, st$lat); pt <- vect(xy, crs = "EPSG:4326")
  cell_xy <- xyFromCell(r1, cellFromXY(r1, xy))
  if (anyNA(cell_xy)) stop(code, " lies outside the cached AGCD grids; widen bbox in script 2d")
  cell_km <- round(distance(xy, cell_xy, lonlat = TRUE) / 1000, 1)
  a <- do.call(rbind, lapply(grids, function(f) {
    r <- rast(f); crs(r) <- "EPSG:4326"
    data.frame(date = as.Date(time(r)), agcd = as.numeric(unlist(terra::extract(r, pt, ID = FALSE))))
  }))
  message(code, ": AGCD cell centre ", paste(round(cell_xy, 3), collapse = ", "), ", ", cell_km, " km from the gauge")

  d <- merge(g, a, by = "date")
  d <- d[is.finite(d$gauge) & is.finite(d$agcd), ]
  d <- d[order(d$date), ]
  win <- range(d$date)                                    # first to last usable gauge day
  daily_s <- stats(d$gauge, d$agcd)

  # monthly totals over the days both have, for months with enough of them
  d$month <- format(d$date, "%Y-%m")
  m <- do.call(rbind, lapply(split(d, d$month), function(z) {
    n_days <- as.integer(format(seq(as.Date(paste0(z$month[1], "-01")), by = "month", length.out = 2)[2] - 1, "%d"))
    if (nrow(z) < min_month * n_days) return(NULL)
    data.frame(month = z$month[1], gauge = sum(z$gauge), agcd = sum(z$agcd))
  }))
  monthly_s <- stats(m$gauge, m$agcd)

  wet <- d$gauge >= 1; awet <- d$agcd >= 1
  n_win <- as.integer(diff(win)) + 1L
  message(sprintf("%s daily: n %d of %d days in %s to %s, r %.2f, RMSE %.2f mm, bias %+.2f mm/day, AGCD/gauge total %.2f; days >= 1 mm: both %d, gauge only %d, AGCD only %d",
                  code, daily_s[["n"]], n_win, win[1], win[2], daily_s[["r"]], daily_s[["rmse"]], daily_s[["bias"]], daily_s[["total_ratio"]],
                  sum(wet & awet), sum(wet & !awet), sum(!wet & awet)))
  message(sprintf("%s monthly: n %d, r %.2f, RMSE %.1f mm, bias %+.1f mm/month, AGCD/gauge total %.2f",
                  code, monthly_s[["n"]], monthly_s[["r"]], monthly_s[["rmse"]], monthly_s[["bias"]], monthly_s[["total_ratio"]]))

  ## ---- the figure
  # every day without a gauge record, shaded as runs
  miss <- setdiff(seq(win[1], win[2], by = "day"), d$date)
  gaps <- if (length(miss)) {
    b <- cumsum(c(1, diff(miss) != 1))
    data.frame(start = as.Date(tapply(miss, b, min), origin = "1970-01-01"), end = as.Date(tapply(miss, b, max), origin = "1970-01-01"))
  }

  ts <- rbind(data.frame(date = d$date, rain =  d$gauge, src = paste(code, "gauge")),
              data.frame(date = d$date, rain = -d$agcd,  src = "AGCD cell"))
  ts$src <- factor(ts$src, levels = c(paste(code, "gauge"), "AGCD cell"))
  brk <- c(-100, -30, -10, -3, -1, 0, 1, 3, 10, 30, 100)
  p_ts <- ggplot(ts, aes(date, ymin = 0, ymax = rain, colour = src)) +
    { if (!is.null(gaps)) geom_rect(data = gaps, aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -Inf, ymax = Inf),
                                    inherit.aes = FALSE, fill = "grey60", alpha = 0.35) } +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_linerange(linewidth = 0.3) +
    scale_colour_manual(values = setNames(c(col_g, col_a), levels(ts$src))) +
    scale_y_continuous(trans = pseudo_log_trans(sigma = 1, base = 10), breaks = brk, labels = abs(brk)) +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    coord_cartesian(xlim = win + c(-20, 20)) +
    labs(x = NULL, y = "Daily rainfall (mm)") +
    base_theme + theme(legend.position = "bottom")

  lim <- ax_lim(d$gauge, d$agcd, plog)
  p_day <- ggplot(d, aes(gauge, agcd)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(shape = 21, size = 1.8, stroke = 0.4, colour = "grey20", fill = alpha("grey40", 0.25)) +
    scale_x_continuous(trans = plog, breaks = c(0, 1, 3, 10, 30, 100), limits = c(0, lim)) +
    scale_y_continuous(trans = plog, breaks = c(0, 1, 3, 10, 30, 100), limits = c(0, lim)) +
    coord_equal() +
    labs(x = paste(code, "gauge (mm/day)"), y = "AGCD (mm/day)") +
    base_theme
  p_day <- add_stats(p_day, daily_s, "days", "mm/day", 2, lim, plog)

  mlim <- ax_lim(m$gauge, m$agcd)
  p_mon <- ggplot(m, aes(gauge, agcd)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    geom_point(shape = 21, size = 2.4, stroke = 0.5, colour = "grey20", fill = alpha("grey40", 0.4)) +
    scale_x_continuous(limits = c(0, mlim)) + scale_y_continuous(limits = c(0, mlim)) +
    coord_equal() +
    labs(x = paste(code, "gauge (mm/month)"), y = "AGCD (mm/month)") +
    base_theme
  p_mon <- add_stats(p_mon, monthly_s, "months", "mm/month", 1, mlim)

  fig <- p_ts / (p_day | p_mon) +
    plot_layout(heights = c(1, 1.25)) +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
  f <- file.path(fig_path, sprintf("AGCD_vs_%s_gauge_rainfall.tif", code))
  ggsave(f, fig, width = 26, height = 24, units = "cm", dpi = 300)
  message("written to ", f)
}

for (code in gauges) compare(code)
