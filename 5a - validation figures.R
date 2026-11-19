################################################################################
#### validation figures: empirical VI models and LEAF-Toolbox SL2P      ####
################################################################################
# Compares two empirical index models (NDVI, EVI) against SL2P at 10 m, plus a
# locally calibrated 10 m product. Figures:
#   1. understorey / overstorey / total contribution (plain and coloured by MVG),
#      with the in situ LAI numbers by MVG
#   2. manuscript scatter panel
#   3. supplementary: SL2P 20 m against the in situ data
#
# Two kinds of model are handled:
#   "index"  a vegetation index, which needs a fitted lm to become LAI
#   "lai"    a model that already predicts LAI, where 1:1 agreement is the expected

library(ggplot2)
library(patchwork)
source("R/mvg_palette.R")   # mvg_cols, shared with 5a_kfold and 5g

end_path  <- "outputs"
plot_path <- "plots"
dir.create(plot_path, showWarnings = FALSE)
txt.size  <- 10
test_sites <- c('MCK', 'CW', 'YU', 'DR')
composites <- c('_closest', '_mean', '_median')
focus_composite <- '_closest'

total_LAI <- readRDS(paste0(end_path, "/paired_ref_and_sat_dataset.rds"))
total_LAI <- total_LAI[!(total_LAI$n_records_understorey < 8 |
                         total_LAI$n_records_canopy < 8 |
                         total_LAI$n_pixels < 8), ]
total_LAI$MVG[total_LAI$MVG == "Naturally bare - sand, rock, claypan, mudflat"] <- 'Naturally bare'
total_LAI$MVG[total_LAI$MVG == "Acacia Open Woodlands"] <- 'Acacia Forests and Woodlands'
stopifnot(all(unique(total_LAI$MVG) %in% names(mvg_cols)))

test_idx <- total_LAI$Site %in% test_sites
train <- total_LAI[!test_idx, ]
test  <- total_LAI[ test_idx, ]
total_LAI$.set <- ifelse(test_idx, "test", "train")   # used by the scatter panels

#############################################################################
# 1. understorey, overstorey and total contribution, and in situ LAI by MVG
#############################################################################
total_LAI$Sample_period <- factor(
  ifelse(total_LAI$Date > as.Date("2022-01-01"), 'March 22', 'October 21'),
  levels = c('October 21', 'March 22'))
total_LAI$Site_MVG <- paste0(total_LAI$Site, '_', total_LAI$MVG)

strata_cols <- c(Overstorey  = "LAI_LX_canopy",
                 Understorey = "LAI_LX_understorey",
                 Total       = "LAI_LX_total")
lvl <- names(strata_cols)

# a vegetation type (one MVG at one site) is paired when it appears in both sampling
# periods. Script 4 gives one reference LAI per vegetation type per date.
per_strata <- table(total_LAI$Site_MVG, total_LAI$Sample_period)
paired_strata <- rownames(per_strata)[rowSums(per_strata > 0) == 2]
message(length(paired_strata), " of ", nrow(per_strata), " ref points sampled in both periods")

p <- total_LAI[total_LAI$Site_MVG %in% paired_strata, ]

strata_long <- do.call(rbind, lapply(lvl, function(s) {
  a <- do.call(rbind, lapply(
    split(p, list(p$Site_MVG, p$Sample_period), drop = TRUE),
    function(x) data.frame(Site_MVG = x$Site_MVG[1], Site = x$Site[1], MVG = x$MVG[1],
                           Sample_period = x$Sample_period[1],
                           LAI = mean(x[[strata_cols[[s]]]], na.rm = TRUE))))
  a$Stratum <- s
  a
}))

strata_long$Stratum <- factor(strata_long$Stratum, levels = lvl)

# Per stratum, the mean paired change in LAI (March 22 - October 21) with its
# standard error, and the change in the two means as a percentage. prop_change is
# the ratio of the means, not the mean of the per-vegetation-type ratios

strata_p <- do.call(rbind, lapply(lvl, function(s) {
  a <- strata_long[strata_long$Stratum == s, ]
  w <- reshape(a[, c("Site_MVG", "Sample_period", "LAI")], direction = "wide",
               idvar = "Site_MVG", timevar = "Sample_period")
  oct <- w[["LAI.October 21"]]; mar <- w[["LAI.March 22"]]   # by name, not position
  ok  <- is.finite(oct) & is.finite(mar)
  d   <- mar[ok] - oct[ok]
  pc  <- mean(mar[ok]) / mean(oct[ok]) - 1
  se  <- sd(d) / sqrt(length(d))
  tt  <- suppressWarnings(wilcox.test(oct[ok], mar[ok], paired = TRUE))
  data.frame(Stratum = s, n = sum(ok), p = tt$p.value,
             mean_Oct21 = mean(oct[ok]), mean_Mar22 = mean(mar[ok]),
             mean_change = mean(d), se = se, prop_change = pc,
             label = sprintf("italic(Delta*LAI)~'='~'%.2f ± %.2f (%.0f %%)'",
                             round(mean(d), 2) + 0, se, round(100 * pc) + 0),
             p_label = if (tt$p.value < 0.001) "p < 0.001" else sprintf("p = %.3f", tt$p.value),
             stringsAsFactors = FALSE)
}))
strata_p$Stratum <- factor(strata_p$Stratum, levels = lvl)
print(strata_p[, c("Stratum", "n", "p", "mean_Oct21", "mean_Mar22", "mean_change", "se",
                   "prop_change")], digits = 3)

lai_max <- max(strata_long$LAI, na.rm = TRUE)
y_top   <- lai_max / 0.85

# by_mvg = TRUE is the same figure with the points and lines coloured by MVG
contribution_fig <- function(by_mvg = FALSE) {
  ln <- if (by_mvg) geom_line(aes(group = Site_MVG, colour = MVG), linewidth = 0.3) else
    geom_line(aes(group = Site_MVG), colour = 'grey70', linewidth = 0.3)
  pt <- if (by_mvg) geom_point(aes(colour = MVG), size = 1.3, alpha = 0.75) else
    geom_point(size = 1.3, colour = "black", alpha = 0.75)
  p <- ggplot(strata_long, aes(x = Sample_period, y = LAI)) +
    ln +
    geom_boxplot(aes(fill = Sample_period), alpha = 0.5, colour = "black",
                 width = 0.55, outlier.shape = NA) +
    pt +
    geom_text(data = strata_p, aes(x = 1.5, y = y_top, label = label), vjust = 1,
              parse = TRUE, inherit.aes = FALSE, size = txt.size / 3.6) +
    geom_text(data = strata_p, aes(x = 1.5, y = y_top - 0.075 * y_top, label = p_label), vjust = 1,
              inherit.aes = FALSE, size = txt.size / 3.6) +
    facet_wrap(~Stratum, nrow = 1) +
    scale_fill_manual(values = c("#E7B800", "#FC4E07")) +
    scale_y_continuous(limits = c(0, y_top), breaks = seq(0, lai_max, by = 0.5)) +
    labs(fill = "Sampling period", x = "", y = expression(LAI~(m^2~m^-2))) +
    theme_bw() +
    theme(axis.text = element_text(size = txt.size, colour = 'black'),
          axis.title = element_text(size = txt.size),
          legend.text = element_text(size = txt.size),
          legend.title = element_text(size = txt.size),
          legend.position = "bottom",
          strip.text = element_text(size = txt.size),
          panel.grid = element_blank())
  if (by_mvg) p <- p + scale_colour_manual(values = mvg_cols, name = NULL) +
    guides(fill = guide_legend(order = 1), colour = guide_legend(order = 2, nrow = 3)) +
    theme(legend.box = "vertical")
  p
}

ggsave(file.path(plot_path, "LAI_strata_contribution_paired.tif"),
       contribution_fig(), width = 20, height = 10, units = "cm", dpi = 300)
ggsave(file.path(plot_path, "LAI_strata_contribution_paired_MVG.tif"),
       contribution_fig(by_mvg = TRUE), width = 20, height = 12, units = "cm", dpi = 300)

# Numbers for the in situ results text. The overstorey share is overstorey /
# total LAI. Both models take site as a random intercept; MVGs are compared with
# Holm-adjusted pairwise contrasts.
se <- function(v) sd(v) / sqrt(length(v))

# late dry season by MVG, one record per vegetation type
oct <- total_LAI[total_LAI$Sample_period == "October 21", ]
oct$MVG <- factor(oct$MVG)
oct$over_pct <- 100 * oct$LAI_LX_canopy / oct$LAI_LX_total
message(sprintf("late dry total LAI: %.2f ± %.2f (mean ± se, n = %d vegetation types at %d sites)",
                mean(oct$LAI_LX_total), se(oct$LAI_LX_total), nrow(oct), length(unique(oct$Site))))

mvg_stats <- do.call(rbind, lapply(split(oct, oct$MVG), function(a) data.frame(
  MVG = a$MVG[1], n = nrow(a), sites = length(unique(a$Site)),
  over_min = min(a$LAI_LX_canopy), over_max = max(a$LAI_LX_canopy),
  over_ratio = max(a$LAI_LX_canopy) / min(a$LAI_LX_canopy),
  over_pct = mean(a$over_pct), over_pct_se = se(a$over_pct))))
print(mvg_stats, digits = 2, row.names = FALSE)

pct_mod <- lmerTest::lmer(over_pct ~ MVG + (1 | Site), data = oct)
print(anova(pct_mod))
b  <- names(lme4::fixef(pct_mod))
pw <- do.call(rbind, combn(levels(oct$MVG), 2, function(g) {
  L <- setNames((b == paste0("MVG", g[1])) - (b == paste0("MVG", g[2])), b)
  r <- lmerTest::contest(pct_mod, L, joint = FALSE)
  data.frame(MVG_a = g[1], MVG_b = g[2], diff_pct = r$Estimate, p = r$`Pr(>|t|)`)
}, simplify = FALSE))
pw$p_holm <- p.adjust(pw$p, "holm")
print(pw[order(pw$p), ], digits = 2, row.names = FALSE)

# total LAI against sampling period, MVG and site, every record
tot_mod <- lmerTest::lmer(LAI_LX_total ~ Sample_period * MVG + (1 | Site), data = total_LAI)
print(anova(tot_mod))
print(lmerTest::ranova(tot_mod))
print(lme4::VarCorr(tot_mod))

# the sampling period effect, from the model without the interaction
add_mod <- lmerTest::lmer(LAI_LX_total ~ Sample_period + MVG + (1 | Site), data = total_LAI)
pe <- summary(add_mod)$coefficients["Sample_periodMarch 22", ]
message(sprintf("March 22 - October 21 total LAI: %.2f ± %.2f (95%% CI %.2f to %.2f, p = %.2g), %.0f%% of the late dry mean",
                pe[["Estimate"]], pe[["Std. Error"]], pe[["Estimate"]] - 1.96 * pe[["Std. Error"]],
                pe[["Estimate"]] + 1.96 * pe[["Std. Error"]], pe[["Pr(>|t|)"]], 100 * pe[["Estimate"]] / mean(oct$LAI_LX_total)))

# r2 partitioned as in 5c: marginal r2 (Nakagawa) of every submodel, split
# individually (Shapley; the interaction enters only after both its terms, so the
# orders period-MVG and MVG-period are averaged) and uniquely (drop one). total
# adds the interaction to the term. Intervals from 50 parametric bootstraps of the
# full model, the same simulated responses for every submodel.
r2 <- function(m) { vf <- var(as.vector(lme4::getME(m, "X") %*% lme4::fixef(m))); vr <- lme4::VarCorr(m)$Site[1]
                    c(vf, vf + vr) / (vf + vr + sigma(m)^2) }
mu   <- as.vector(lme4::getME(tot_mod, "X") %*% lme4::fixef(tot_mod))
su   <- attr(lme4::VarCorr(tot_mod)$Site, "stddev")
site <- as.integer(factor(total_LAI$Site))
r2_b <- lapply(c(P = "Sample_period", M = "MVG", PM = "Sample_period + MVG", full = "Sample_period * MVG"), function(f) {
  m <- lme4::lmer(as.formula(paste("LAI_LX_total ~", f, "+ (1 | Site)")), data = total_LAI)
  rbind(r2(m), t(sapply(1:50, function(b) { set.seed(b)
    r2(suppressMessages(suppressWarnings(lme4::refit(m, newresp = mu + rnorm(max(site), 0, su)[site] +
                                                       rnorm(length(mu), 0, sigma(tot_mod)))))) })))
})
part <- function(i) {
  R   <- sapply(r2_b, `[`, i, 1); cond <- r2_b$full[i, 2]; int <- R[["full"]] - R[["PM"]]
  ind <- c((R[["P"]] + R[["PM"]] - R[["M"]]) / 2, (R[["M"]] + R[["PM"]] - R[["P"]]) / 2)
  uni <- c(R[["PM"]] - R[["M"]], R[["PM"]] - R[["P"]])
  data.frame(partition = rep(c("model", "individual", "unique"), c(3, 2, 2)),
             term = c("Full model", "Fixed effects", "Site effect", rep(c("Sampling period", "MVG"), 2)),
             solo = c(cond, R[["full"]], cond - R[["full"]], ind, uni), total = c(NA, NA, NA, ind + int, uni + int))
}
reps <- lapply(1:51, part)
ci   <- function(col, q) apply(sapply(reps[-1], `[[`, col), 1, quantile, q, na.rm = TRUE)
vp   <- cbind(reps[[1]], solo_lower = ci("solo", 0.025), solo_upper = ci("solo", 0.975),
              total_lower = ci("total", 0.025), total_upper = ci("total", 0.975))
print(cbind(vp[1:2], round(vp[-(1:2)], 3)), row.names = FALSE)

# change between the sampling periods in the vegetation types of the figure, and
# those whose understorey change exceeds the mean understorey LAI
chg <- reshape(strata_long[, c("Site_MVG", "Stratum", "Sample_period", "LAI")], direction = "wide",
               idvar = c("Site_MVG", "Stratum"), timevar = "Sample_period")
chg$change <- chg[["LAI.March 22"]] - chg[["LAI.October 21"]]
tot <- chg[chg$Stratum == "Total", ]
und <- chg[chg$Stratum == "Understorey", ]
ovr <- chg[chg$Stratum == "Overstorey", ][match(und$Site_MVG, chg$Site_MVG[chg$Stratum == "Overstorey"]), ]
message(sprintf("total LAI rose in %d of %d vegetation types sampled in both periods (change %.2f to %.2f)",
                sum(tot$change > 0), nrow(tot), min(tot$change), max(tot$change)))
message(sprintf("mean absolute change: understorey %.2f ± %.2f, overstorey %.2f ± %.2f (paired Wilcoxon p = %.2g)",
                mean(abs(und$change)), se(abs(und$change)), mean(abs(ovr$change)), se(abs(ovr$change)),
                wilcox.test(abs(und$change), abs(ovr$change), paired = TRUE)$p.value))
und_LAI <- strata_long$LAI[strata_long$Stratum == "Understorey"]
big <- abs(und$change) > mean(und_LAI)
message(sprintf("%d vegetation types changed by more than the mean understorey LAI (%.2f ± %.2f), %.0f%% of the summed understorey change: %s",
                sum(big), mean(und_LAI), se(und_LAI), 100 * sum(abs(und$change[big])) / sum(abs(und$change)),
                paste(und$Site_MVG[big], collapse = ", ")))

#############################################################################
# 2. helpers: gg_xy() and model_table()
#############################################################################
pull_p <- function(m) {
  f <- summary(m)$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = FALSE); attributes(p) <- NULL
  if (p < 0.001) '<0.001' else signif(p, 2)
}

# The candidate functional forms for turning veg indices into LAI. Shared by
# model_table() and the manuscript fig
vi_forms <- function(v, ref = "LAI_LX_total") c(
  linear       = paste(ref, "~", v),
  logx         = paste(ref, "~ log(", v, "+1)"),
  rootx        = paste(ref, "~ sqrt(", v, "+1)"),
  mixed_linear = paste(ref, "~", v, "+ MVG"),
  mixed_log    = paste(ref, "~ log(", v, "+1) + MVG"),
  mixed_root   = paste(ref, "~ sqrt(", v, "+1) + MVG"))
index_forms <- c("linear", "logx", "rootx")   # the forms that predict from the index alone
# the best form has the lowest validation RMSE to 3 significant figures; forms tied
# at that go in the order listed above
best_rank <- function(rmse) order(signif(rmse, 3))
# a model's name in tables and figures: "NDVI (filtered), linear", "SL2P 10 m (calibrated)"
model_label <- function(stem, call) ifelse(grepl("^leaf_", stem),
  sprintf("%s (%s)", sub("^leaf_(\\d+)m$", "SL2P \\1 m", stem), ifelse(call == "identity", "uncalibrated", "calibrated")),
  sprintf("%s (%s), %s", sub("^raw", "", stem), ifelse(grepl("^raw", stem), "Sen2Cor", "NBART"), call))

# table formatting: numbers to 2 dp with a true minus, signed bias, the fitted
# equation of any form, and a csv (with a BOM so Excel keeps the symbols) plus an
# xlsx with bold headers and rules, ready to paste into Word
num <- function(x, dp = 2) sub("^-", "−", formatC(x, format = "f", digits = dp))
sgn <- function(x) { x <- round(x, 2) + 0   # + 0 turns -0 into 0, so no "−0.00"
  ifelse(is.na(x), "", paste0(ifelse(x < 0, "", "+"), num(x))) }
eqn <- function(m, x, form = "linear") {
  b <- unname(coef(m))
  t <- switch(sub("^mixed_", "", form), logx = , log = sprintf("ln(%s + 1)", x), rootx = , root = sprintf("√(%s + 1)", x), x)
  if (grepl("^mixed_", form)) return(sprintf("LAI = %s × %s + MVG intercept", num(b[2]), t))
  sprintf("LAI = %s %s %s × %s", num(b[1]), if (b[2] < 0) "−" else "+", num(abs(b[2])), t)
}
# top: an optional header row above the column names, spanning top_cols
save_table <- function(x, file, top = NULL, top_cols = NULL) {
  out <- rbind(if (!is.null(top)) replace(rep("", ncol(x)), top_cols[1], top), names(x), as.matrix(x))
  con <- file(file.path(plot_path, paste0(file, ".csv")), "w", encoding = "UTF-8")
  writeLines("﻿", con, sep = ""); write.table(out, con, sep = ",", row.names = FALSE, col.names = FALSE); close(con)
  h  <- seq_len(nrow(out) - nrow(x)); cc <- seq_len(ncol(x))
  st <- function(...) openxlsx::createStyle(...)
  wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, "Table", gridLines = FALSE)
  openxlsx::writeData(wb, 1, out, colNames = FALSE)
  openxlsx::addStyle(wb, 1, st(textDecoration = "bold"), rows = h, cols = cc, gridExpand = TRUE)
  openxlsx::addStyle(wb, 1, st(border = "top"), rows = 1, cols = cc, stack = TRUE)
  openxlsx::addStyle(wb, 1, st(border = "bottom"), rows = max(h), cols = cc, stack = TRUE)
  openxlsx::addStyle(wb, 1, st(border = "bottom"), rows = nrow(out), cols = cc, stack = TRUE)
  if (!is.null(top)) {
    openxlsx::mergeCells(wb, 1, cols = top_cols, rows = 1)
    openxlsx::addStyle(wb, 1, st(halign = "center", border = "bottom"), rows = 1, cols = top_cols, stack = TRUE)
  }
  openxlsx::setColWidths(wb, 1, cc, "auto")
  openxlsx::saveWorkbook(wb, file.path(plot_path, paste0(file, ".xlsx")), overwrite = TRUE)
}

# Reference vs model scatter. Every point is drawn, but the statistics come from
# the held-out test sites only, so the panel and the summary table agree.
# Training points are hollow to keep that distinction visible.
gg_xy <- function(d, pred, label, ref = "LAI_LX_total", lims = NULL, fs = txt.size) {
  d <- d[is.finite(d[[pred]]) & is.finite(d[[ref]]), ]
  if (!nrow(d)) return(NULL)
  if (is.null(lims)) lims <- range(c(d[[ref]], d[[pred]]), na.rm = TRUE)
  te <- d[d$.set == "test", ]
  # one line per annotate(), so r^2 can be set as an expression
  ann <- if (nrow(te) >= 3) {
    c(sprintf("test n = %d (train %d)", nrow(te), nrow(d) - nrow(te)),
      sprintf("bias = %+.2f", mean(te[[pred]] - te[[ref]])),
      sprintf("RMSE = %.2f", sqrt(mean((te[[pred]] - te[[ref]])^2))),
      sprintf("italic(r)^2 == %.2f", suppressWarnings(cor(te[[pred]], te[[ref]]))^2))
  } else sprintf("n = %d (no test sites)", nrow(d))
  step <- diff(lims) * 0.07
  g <- ggplot(d, aes(x = .data[[ref]], y = .data[[pred]])) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = 'grey40') +
    geom_smooth(method = "lm", formula = y ~ x, colour = 'black', alpha = 0.4,
                linewidth = 0.5) +
    geom_point(aes(colour = MVG, shape = .set), size = 1.4) +
    scale_colour_manual(values = mvg_cols) +
    scale_shape_manual(values = c(train = 1, test = 19), guide = "none")
  for (k in seq_along(ann))
    g <- g + annotate("text", x = lims[1], y = lims[2] - (k - 1) * step, label = ann[k],
                      parse = (k == 4), hjust = 0, vjust = 1, size = fs / 3.2)
  g +
    coord_equal(xlim = lims, ylim = lims) +
    labs(x = expression(reference~LAI~(m^2~m^-2)),
         y = expression(modelled~LAI~(m^2~m^-2)), title = label) +
    theme_bw() +
    theme(plot.title = element_text(size = fs, face = "plain"),
          legend.position = 'bottom', legend.title = element_blank(),
          legend.text = element_text(size = fs),
          panel.grid = element_blank(),
          axis.text = element_text(size = fs, colour = 'black'),
          axis.title = element_text(size = fs)) +
    guides(colour = guide_legend(ncol = 3))
}

# Summary statistics for one predictor. kind = "index" fits the linear/log/root
# (+ MVG) forms; "lai" scores the prediction on the 1:1 line and adds a
# calibration lm, whose slope and intercept separate offset from scale error.
model_table <- function(stem, kind = c("index", "lai"), composite = focus_composite,
                        training_data = train, test_data = test,
                        ref = "LAI_LX_total") {
  kind <- match.arg(kind)
  v <- paste0(stem, composite)
  if (!v %in% names(training_data)) return(NULL)
  ok <- function(d) d[is.finite(d[[v]]) & is.finite(d[[ref]]), ]
  tr <- ok(training_data); te <- ok(test_data)
  if (nrow(tr) < 5 || !nrow(te)) return(NULL)

  x_lab <- if (kind == "lai") "SL2P" else sub("^raw", "", stem)
  score <- function(actual, predicted, call, mod = NULL) data.frame(
    Model      = paste0(stem, "_", call),
    Model.stem = stem,
    Model.call = call,
    composite  = sub("^_", "", composite),
    eq         = if (is.null(mod)) paste("LAI =", x_lab) else eqn(mod, x_lab, call),
    n_train    = nrow(tr),
    n_val      = length(actual),
    r2_train   = if (is.null(mod)) NA else summary(mod)$r.squared,
    p          = if (is.null(mod)) NA else as.character(pull_p(mod)),
    r2_val     = suppressWarnings(cor(predicted, actual))^2,
    bias       = mean(predicted - actual),
    RMSE       = sqrt(mean((predicted - actual)^2)),
    NRMSE      = sqrt(mean((predicted - actual)^2)) / mean(actual),
    stringsAsFactors = FALSE)

  if (kind == "lai") {
    out <- list(score(te[[ref]], te[[v]], "identity"))
    cal <- lm(as.formula(paste(ref, "~", v)), data = tr)
    out[[2]] <- score(te[[ref]], predict(cal, newdata = te), "calibrated", cal)
    return(do.call(rbind, out))
  }

  forms <- vi_forms(v, ref)
  do.call(rbind, lapply(names(forms), function(nm) {
    m <- try(lm(as.formula(forms[[nm]]), data = tr), silent = TRUE)
    if (inherits(m, "try-error")) return(NULL)
    p <- try(predict(m, newdata = te), silent = TRUE)
    if (inherits(p, "try-error")) return(NULL)
    score(te[[ref]], p, nm, m)
  }))
}

#############################################################################
# 3. what models are available
#############################################################################
nm <- names(total_LAI)
all_stems  <- unique(sub("_mean$", "", grep("_mean$", nm, value = TRUE)))
raw_stems  <- grep("^raw",  all_stems, value = TRUE)
leaf_stems <- grep("^leaf.*m$", all_stems, value = TRUE)   # not the _err uncertainty layers
# whatever is left is a filtered vegetation index. The merged dataset also
# carries LAI products from the other scripts, which end in _mean too; scoring
# one as an index would quietly fit an lm to it and report it twice
vi_stems   <- grep("^(raw|leaf|prosail)", all_stems, value = TRUE, invert = TRUE)

message("filtered VIs : ", length(vi_stems))
message("raw VIs      : ", length(raw_stems))
message("LEAF         : ", paste(leaf_stems, collapse = ", "))

#############################################################################
# 4. summary tables, one per composite, in plots/model_summary_stats
#############################################################################
# every model, and the best of each index and product: the index form (no MVG
# intercepts) with the lowest validation RMSE, and SL2P both uncalibrated and calibrated
stats_dir <- file.path(plot_path, "model_summary_stats")
dir.create(stats_dir, showWarnings = FALSE)
all_stats <- do.call(rbind, lapply(composites, function(comp) rbind(
  do.call(rbind, lapply(c(vi_stems, raw_stems), model_table, kind = "index", composite = comp)),
  do.call(rbind, lapply(leaf_stems,             model_table, kind = "lai",   composite = comp)))))
best_of <- function(s) do.call(rbind, lapply(unique(s$Model.stem), function(st) {
  a <- s[s$Model.stem == st, ]
  if (all(a$Model.call %in% c("identity", "calibrated"))) return(a)
  a <- a[a$Model.call %in% index_forms, ]; a[best_rank(a$RMSE)[1], ] }))
fmt_stats <- function(s) data.frame(
  Model = model_label(s$Model.stem, s$Model.call), "Regression model" = s$eq, "n (training)" = s$n_train, "n (validation)" = s$n_val,
  "r² (training)" = ifelse(is.na(s$r2_train), "–", num(s$r2_train)),
  p = ifelse(is.na(s$p), "–", sub("<", "< ", s$p)), "r² (validation)" = num(s$r2_val),
  Bias = sgn(s$bias), RMSE = num(s$RMSE), NRMSE = num(s$NRMSE), check.names = FALSE)

for (comp in unique(all_stats$composite)) {
  s <- all_stats[all_stats$composite == comp, ]
  save_table(fmt_stats(s),          file.path("model_summary_stats", paste0("model_summary_stats_", comp)))
  save_table(fmt_stats(best_of(s)), file.path("model_summary_stats", paste0("model_summary_stats_best_", comp)))
}
save_table(cbind(Composite = all_stats$composite, fmt_stats(all_stats)),
           file.path("model_summary_stats", "model_summary_stats_all_composites"))

#############################################################################
# 5. inspection plots: the best model of each index and product, _closest
#############################################################################
# reference LAI against the predictor, hollow = training, filled = validation; red =
# the fitted model (the calibration line for SL2P), dashed = 1:1 for uncalibrated SL2P
best_closest <- best_of(all_stats[all_stats$composite == sub("^_", "", focus_composite), ])
pdf(file.path(stats_dir, "model_inspection_best_closest.pdf"), width = 11, height = 8)
# titles: validation RMSE; r of predictor and reference LAI over every record; p of the fit
op <- par(mfrow = c(3, 4), mar = c(4, 4, 4, 1))
for (i in seq_len(nrow(best_closest))) {
  b  <- best_closest[i, ]; v <- paste0(b$Model.stem, focus_composite)
  d  <- total_LAI[is.finite(total_LAI[[v]]) & is.finite(total_LAI$LAI_LX_total), ]
  tr <- d$.set == "train"
  p  <- if (is.na(b$p)) "" else paste(",", if (grepl("<", b$p)) sub("<", "p < ", b$p) else paste("p =", b$p))
  x  <- if (grepl("^leaf_", b$Model.stem)) sub("^leaf_(\\d+)m$", "SL2P \\1 m", b$Model.stem) else sub("^raw", "", b$Model.stem)
  plot(d[[v]], d$LAI_LX_total, pch = ifelse(tr, 1, 19), cex = 0.8, col = "grey20", xlab = x, ylab = "reference LAI",
       main = sprintf("%s\nvalidation RMSE %.2f\nr = %.2f%s", model_label(b$Model.stem, b$Model.call), b$RMSE,
                      cor(d[[v]], d$LAI_LX_total), p), cex.main = 0.9)
  if (b$Model.call == "identity") { abline(0, 1, lty = 2); next }
  form <- if (b$Model.call == "calibrated") "linear" else b$Model.call
  m  <- lm(as.formula(vi_forms(v)[[form]]), data = d[tr, ])
  xs <- seq(min(d[[v]]), max(d[[v]]), length.out = 100)
  lines(xs, predict(m, newdata = setNames(data.frame(xs), v)), col = "red", lwd = 1.5)
}
par(op)
dev.off()

#############################################################################
# 6. manuscript scatter panel
#############################################################################
# NDVI and EVI only, raw before filtered: both are standard and directly
# comparable between the raw and filtered runs. MSAVI2 edged them on RMSE but is
# a less useful thing to report.
panel_index <- c("NDVI", "EVI")

# Locally calibrated SL2P: reference LAI regressed on the retrieval over the
# training sites, then applied everywhere. Fitted on training sites only, so the
# held-out statistics stay honest. One fit per resolution.
fit_cal <- function(stem) {
  v <- paste0(stem, focus_composite)
  if (!v %in% names(train)) return(NULL)
  d <- train[is.finite(train[[v]]) & is.finite(train$LAI_LX_total), ]
  if (nrow(d) >= 5) lm(as.formula(paste("LAI_LX_total ~", v)), data = d)
}
leaf_cals <- list(leaf_10m = fit_cal("leaf_10m"), leaf_20m = fit_cal("leaf_20m"))
leaf_cal  <- leaf_cals$leaf_10m
apply_cal <- function(x, stem = "leaf_10m") {
  m <- leaf_cals[[stem]]
  if (is.null(m)) return(x * NA)
  as.numeric(predict(m, newdata = setNames(data.frame(x), paste0(stem, focus_composite))))
}

# index -> LAI models for script 5b, fitted on the training sites over the
# _closest composite and keyed <index>_<form>. Only the forms that predict from
# the index alone are kept; the mixed_* forms also need MVG, which 5b does not
# carry into predict(). The unfiltered series is fitted the same way and keyed
# raw<index>_<form>.
lin_mods <- list()
for (v in c(panel_index, paste0("raw", panel_index))) {
  col <- paste0(v, "_closest")
  if (!col %in% names(train)) next
  tr <- train[is.finite(train[[col]]) & is.finite(train$LAI_LX_total), ]
  if (nrow(tr) < 5) next
  for (nm in index_forms) {
    m <- try(lm(as.formula(vi_forms(col)[[nm]]), data = tr), silent = TRUE)
    if (!inherits(m, "try-error")) lin_mods[[paste0(v, "_", nm)]] <- m
  }
}
saveRDS(lin_mods, file.path(end_path, "Optimum_closest_VI_lin_reg_LAI_mod.rds"))
message("VI -> LAI models saved: ", paste(names(lin_mods), collapse = ", "))

# The one index model used wherever a filtered index becomes LAI: lowest
# held-out RMSE over NDVI and EVI. Script 5b reads it, so the validation
# timeseries and the phase analysis share the same model.
vi_ts <- local({
  s <- do.call(rbind, lapply(panel_index, model_table, kind = "index",
                             composite = focus_composite))
  s <- s[s$Model.call %in% index_forms, ]
  s <- s[best_rank(s$RMSE)[1], c("Model.stem", "Model.call", "RMSE")]
  s$RMSE <- signif(s$RMSE, 3); s
})
vi_ts$key <- paste0(vi_ts$Model.stem, "_", vi_ts$Model.call)
stopifnot(vi_ts$key %in% names(lin_mods))
saveRDS(list(stem = vi_ts$Model.stem, form = vi_ts$Model.call, key = vi_ts$key,
             test_RMSE = vi_ts$RMSE), file.path(end_path, "Optimum_closest_VI_choice.rds"))
message("chosen VI model: ", vi_ts$key, " (held-out RMSE ", vi_ts$RMSE, ")")

# Every panel's prediction is worked out first so they can share one axis range.
# Per-panel ranges make the models look more alike than they are - a model biased
# high just gets a taller axis and still fills the box.
panel_spec <- function(stem, label, kind = "lai") {
  v <- paste0(stem, focus_composite)
  if (!v %in% names(total_LAI)) return(NULL)
  d <- total_LAI; m <- NULL
  if (kind == "index") {
    s  <- model_table(stem, kind = "index", composite = focus_composite)
    if (is.null(s)) return(NULL)
    s  <- s[s$Model.call %in% index_forms, ]      # no vegetation class intercepts
    cl <- s$Model.call[best_rank(s$RMSE)][1]      # best form on held-out RMSE
    m  <- try(lm(as.formula(vi_forms(v)[[cl]]),
                 data = train[is.finite(train[[v]]), ]), silent = TRUE)
    if (inherits(m, "try-error")) return(NULL)
    d$.pred <- suppressWarnings(predict(m, newdata = d))
    label <- sprintf("%s (%s, %s)", sub("^raw", "", stem), label, cl)
  } else if (kind == "cal") {
    m <- leaf_cals[[stem]]
    d$.pred <- apply_cal(d[[v]], stem)
  } else d$.pred <- d[[v]]
  if (all(!is.finite(d$.pred))) return(NULL)
  list(d = d, label = label, m = m, form = if (kind == "index") cl else "linear")
}

# pad_r (cm) is extra right margin for a title running past its panel, so the edges stay even
draw_panels <- function(specs, file, ncol, width, height, fs = txt.size, pad_r = 0) {
  specs <- Filter(Negate(is.null), specs)
  lims  <- range(unlist(lapply(specs, function(s) c(s$d$.pred, s$d$LAI_LX_total))),
                 na.rm = TRUE, finite = TRUE)
  p <- Filter(Negate(is.null),
              lapply(specs, function(s) gg_xy(s$d, ".pred", s$label, lims = lims, fs = fs)))
  if (!length(p)) return(invisible(NULL))
  ggsave(file.path(plot_path, file),
         Reduce(`+`, p) + plot_layout(ncol = ncol, guides = 'collect') +
           plot_annotation(theme = theme(plot.margin = margin(0.2, 0.2 + pad_r, 0.2, 0.2, "cm"))) &
           theme(legend.position = 'bottom', legend.title = element_blank()),
         width = width, height = height, units = "cm", dpi = 300)
  message(file, ": ", length(p), " panels")
}

# Laid out so each column is a pair: Sen2Cor above filtered for the two
# indices, uncalibrated SL2P above calibrated. patchwork fills by row at ncol = 3,
# so the top row is the three raw/uncalibrated panels and the bottom row their partners.
main_specs <- c(lapply(paste0("raw", panel_index), panel_spec, label = "Sen2Cor", kind = "index"),
                list(panel_spec("leaf_10m", "LEAF SL2P 10 m (uncalibrated)")),
                lapply(panel_index, panel_spec, label = "NBART", kind = "index"),
                list(panel_spec("leaf_10m", "LEAF SL2P 10 m (calibrated, linear)", kind = "cal")))
draw_panels(main_specs, "model_vs_reference_scatter.tif", ncol = 3, width = 22.2, height = 18, fs = 1.1 * txt.size,
            pad_r = 0.95)

# supplementary: the 20 m product against the same in situ data, as it comes and calibrated
draw_panels(list(panel_spec("leaf_20m", "LEAF SL2P 20 m (uncalibrated)"),
                 panel_spec("leaf_20m", "LEAF SL2P 20 m (calibrated, linear)", kind = "cal")),
            "supp_leaf20m_vs_reference.tif", ncol = 2, width = 20, height = 12.5)

for (s in names(leaf_cals)) if (!is.null(leaf_cals[[s]]))
  message("SL2P ", sub("leaf_", "", s), " calibration: LAI = ", round(coef(leaf_cals[[s]])[1], 3), " + ",
          round(coef(leaf_cals[[s]])[2], 3), " x retrieval  (n = ", nobs(leaf_cals[[s]]),
          ", R2 = ", round(summary(leaf_cals[[s]])$r.squared, 3), ")")

#############################################################################
# 7. manuscript tables
#############################################################################
# The numbers behind each scatter panel (validation sites unless marked
# training) and Table 2: the regression models with their validation bias by
# reference LAI class (save_table() in section 2).
names(main_specs) <- c("NDVI (Sen2Cor)", "EVI (Sen2Cor)", "SL2P 10 m (uncalibrated)",
                       "NDVI (NBART)", "EVI (NBART)", "SL2P 10 m (calibrated)")
x_name   <- c("NDVI", "EVI", "SL2P", "NDVI", "EVI", "SL2P")
lai_bins <- c(0, 0.5, 1, 1.5, 2)

val <- lapply(main_specs, function(s) {
  d <- s$d[is.finite(s$d$.pred) & is.finite(s$d$LAI_LX_total), ]
  list(n_train = sum(d$.set == "train"), te = d[d$.set == "test", ], m = s$m, form = s$form)
})
key_stats <- do.call(rbind, Map(function(v, x, nm) {
  e <- v$te$.pred - v$te$LAI_LX_total; fit <- !is.null(v$m)
  data.frame(Model = nm, "Regression model" = if (fit) eqn(v$m, x, v$form) else paste("LAI =", x),
             "n (training)" = v$n_train, "n (validation)" = nrow(v$te),
             "r² (training)" = if (fit) num(summary(v$m)$r.squared) else "–",
             p = if (fit) sub("<", "< ", pull_p(v$m)) else "–",
             "r² (validation)" = num(cor(v$te$.pred, v$te$LAI_LX_total)^2),
             Bias = sgn(mean(e)), RMSE = num(sqrt(mean(e^2))),
             NRMSE = num(sqrt(mean(e^2)) / mean(v$te$LAI_LX_total)), check.names = FALSE)
}, val, x_name, names(main_specs)))
save_table(key_stats, "validation_key_stats")

# bias by reference LAI class at the validation sites, n per class in the header
lai_cls  <- function(v) cut(v$te$LAI_LX_total, lai_bins, include.lowest = TRUE)
cls_n    <- table(lai_cls(val[[1]]))
bias_cls <- t(sapply(val, function(v) sgn(tapply(v$te$.pred - v$te$LAI_LX_total, lai_cls(v), mean))))
colnames(bias_cls) <- sprintf("%s%s – %s (n = %d)", c("", rep("> ", length(cls_n) - 1)),
                              head(lai_bins, -1), lai_bins[-1], cls_n)
table2 <- cbind(key_stats[, c("Model", "Regression model", "NRMSE", "p")], bias_cls, row.names = NULL)
save_table(table2[order(grepl("SL2P", table2$Model)), ], "Table2_validation_models",   # SL2P rows last
           top = "LAI bias", top_cols = 4 + seq_along(cls_n))
