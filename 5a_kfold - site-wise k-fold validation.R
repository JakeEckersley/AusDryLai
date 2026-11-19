################################################################################
#### 5a_kfold - validation by site-wise k-fold cross-validation           ####
################################################################################
# Script 5a's validation with every site held out in turn instead of four fixed
# test sites. The sites are split at random into k folds; each fold's records are
# predicted by models fitted on the other folds, so both campaigns at a site are
# always held out together. The split is repeated n_rep times. Validation
# statistics are the mean over repeats of those from the out-of-fold predictions;
# the regression model, training r2 and p are the fit to every record. Figures
# show each record's mean out-of-fold prediction. Written to plots/kfold.
#   k = 5: 15 sites gives 3 per fold and ~80% of the records to fit each model;
#   repeating the random split averages out which sites happen to share a fold.
# Models with MVG intercepts cannot predict a class missing from the training
# folds; those records are left out of that repeat's statistics.

library(ggplot2)
library(patchwork)
source("R/mvg_palette.R")   # mvg_cols, shared with 5a and 5g

end_path   <- "outputs"
plot_path  <- file.path("plots", "kfold")
dir.create(file.path(plot_path, "model_summary_stats"), showWarnings = FALSE, recursive = TRUE)
txt.size   <- 10
composites <- c('_closest', '_mean', '_median')
focus_composite <- '_closest'
k     <- 5
n_rep <- 100

total_LAI <- readRDS(file.path(end_path, "paired_ref_and_sat_dataset.rds"))
total_LAI <- total_LAI[!(total_LAI$n_records_understorey < 8 |
                         total_LAI$n_records_canopy < 8 |
                         total_LAI$n_pixels < 8), ]
total_LAI$MVG[total_LAI$MVG == "Naturally bare - sand, rock, claypan, mudflat"] <- 'Naturally bare'
total_LAI$MVG[total_LAI$MVG == "Acacia Open Woodlands"] <- 'Acacia Forests and Woodlands'

sites <- sort(unique(total_LAI$Site))
set.seed(1)
folds <- replicate(n_rep, setNames(sample(rep_len(seq_len(k), length(sites))), sites), simplify = FALSE)
message(nrow(total_LAI), " records at ", length(sites), " sites; ", k, " folds of ",
        paste(range(table(folds[[1]])), collapse = "-"), " sites, ", n_rep, " random splits")

#############################################################################
# 1. helpers (as in 5a)
#############################################################################
pull_p <- function(m) {
  f <- summary(m)$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = FALSE); attributes(p) <- NULL
  if (p < 0.001) '<0.001' else signif(p, 2)
}
vi_forms <- function(v, ref = "LAI_LX_total") c(
  linear       = paste(ref, "~", v),
  logx         = paste(ref, "~ log(", v, "+1)"),
  rootx        = paste(ref, "~ sqrt(", v, "+1)"),
  mixed_linear = paste(ref, "~", v, "+ MVG"),
  mixed_log    = paste(ref, "~ log(", v, "+1) + MVG"),
  mixed_root   = paste(ref, "~ sqrt(", v, "+1) + MVG"))
index_forms <- c("linear", "logx", "rootx")
best_rank   <- function(rmse) order(signif(rmse, 3))
# a model's name in tables and figures: "NDVI (NBART), linear", "SL2P 10 m (calibrated)"
model_label <- function(stem, call) ifelse(grepl("^leaf_", stem),
  sprintf("%s (%s)", sub("^leaf_(\\d+)m$", "SL2P \\1 m", stem), ifelse(call == "identity", "uncalibrated", "calibrated")),
  sprintf("%s (%s), %s", sub("^raw", "", stem), ifelse(grepl("^raw", stem), "Sen2Cor", "NBART"), call))

num <- function(x, dp = 2) sub("^-", "−", formatC(x, format = "f", digits = dp))
sgn <- function(x) { x <- round(x, 2) + 0   # + 0 turns -0 into 0, so no "−0.00"
  ifelse(is.na(x), "", paste0(ifelse(x < 0, "", "+"), num(x))) }
eqn <- function(m, x, form = "linear") {
  b <- unname(coef(m))
  t <- switch(sub("^mixed_", "", form), logx = , log = sprintf("ln(%s + 1)", x), rootx = , root = sprintf("√(%s + 1)", x), x)
  if (grepl("^mixed_", form)) return(sprintf("LAI = %s × %s + MVG intercept", num(b[2]), t))
  sprintf("LAI = %s %s %s × %s", num(b[1]), if (b[2] < 0) "−" else "+", num(abs(b[2])), t)
}
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
fmt_stats <- function(s, model = model_label(s$Model.stem, s$Model.call)) data.frame(
  Model = model, "Regression model" = s$eq, "n (training)" = s$n_train, "n (validation)" = s$n_val,
  "r² (training)" = ifelse(is.na(s$r2_train), "–", num(s$r2_train)),
  p = ifelse(is.na(s$p), "–", sub("<", "< ", s$p)), "r² (validation)" = num(s$r2_val),
  Bias = sgn(s$bias), RMSE = num(s$RMSE), NRMSE = num(s$NRMSE), check.names = FALSE)

#############################################################################
# 2. out-of-fold predictions and their statistics
#############################################################################
# one column per repeat; fit(train, test) returns the test predictions
oof <- function(d, fit) sapply(seq_len(n_rep), function(r) {
  f <- folds[[r]][d$Site]; p <- rep(NA_real_, nrow(d))
  for (j in seq_len(k)) { te <- f == j; if (any(te)) p[te] <- fit(d[!te, ], d[te, ]) }
  p
})
lm_fit <- function(form, mixed) function(tr, te) {
  m  <- lm(form, data = tr)
  ok <- if (mixed) te$MVG %in% tr$MVG else rep(TRUE, nrow(te))
  p  <- rep(NA_real_, nrow(te)); if (any(ok)) p[ok] <- predict(m, newdata = te[ok, ]); p
}
cv_score <- function(P, y) rowMeans(apply(P, 2, function(p) {
  ok <- is.finite(p); e <- p[ok] - y[ok]
  c(n = sum(ok), r2 = cor(p[ok], y[ok])^2, bias = mean(e), RMSE = sqrt(mean(e^2)), NRMSE = sqrt(mean(e^2)) / mean(y[ok]))
}))

# kind = "index" fits the linear/log/root (+ MVG) forms; "lai" scores the product
# as it comes and calibrated by a linear fit. keep = TRUE also returns each
# record's mean out-of-fold prediction.
model_table <- function(stem, kind = c("index", "lai"), composite = focus_composite, keep = FALSE) {
  kind <- match.arg(kind)
  v <- paste0(stem, composite)
  if (!v %in% names(total_LAI)) return(NULL)
  d <- total_LAI[is.finite(total_LAI[[v]]) & is.finite(total_LAI$LAI_LX_total), ]
  if (nrow(d) < 10) return(NULL)
  x_lab <- if (kind == "lai") "SL2P" else sub("^raw", "", stem)
  row <- function(call, P, mod) {
    s <- cv_score(P, d$LAI_LX_total)
    out <- data.frame(Model = paste0(stem, "_", call), Model.stem = stem, Model.call = call,
                      composite = sub("^_", "", composite),
                      eq = if (is.null(mod)) paste("LAI =", x_lab) else eqn(mod, x_lab, call),
                      n_train = nrow(d), n_val = round(s[["n"]]),
                      r2_train = if (is.null(mod)) NA else summary(mod)$r.squared,
                      p = if (is.null(mod)) NA else as.character(pull_p(mod)),
                      r2_val = s[["r2"]], bias = s[["bias"]], RMSE = s[["RMSE"]], NRMSE = s[["NRMSE"]],
                      stringsAsFactors = FALSE)
    if (keep) attr(out, "pred") <- data.frame(d, .pred = rowMeans(P, na.rm = TRUE))
    out
  }
  if (kind == "lai") {
    f <- as.formula(paste("LAI_LX_total ~", v))
    return(list(identity   = row("identity", matrix(d[[v]], nrow(d), n_rep), NULL),
                calibrated = row("calibrated", oof(d, lm_fit(f, FALSE)), lm(f, data = d))))
  }
  forms <- vi_forms(v)
  setNames(lapply(names(forms), function(nm) {
    f <- as.formula(forms[[nm]])
    m <- try(lm(f, data = d), silent = TRUE); if (inherits(m, "try-error")) return(NULL)
    P <- try(oof(d, lm_fit(f, grepl("^mixed", nm))), silent = TRUE); if (inherits(P, "try-error")) return(NULL)
    row(nm, P, m)
  }), names(forms))
}
bind <- function(x) do.call(rbind, Filter(Negate(is.null), x))

#############################################################################
# 3. what models are available
#############################################################################
nm <- names(total_LAI)
all_stems  <- unique(sub("_mean$", "", grep("_mean$", nm, value = TRUE)))
raw_stems  <- grep("^raw",  all_stems, value = TRUE)
leaf_stems <- grep("^leaf.*m$", all_stems, value = TRUE)   # not the _err uncertainty layers
vi_stems   <- grep("^(raw|leaf|prosail)", all_stems, value = TRUE, invert = TRUE)

#############################################################################
# 4. summary tables, one per composite, in plots/kfold/model_summary_stats
#############################################################################
t0 <- Sys.time()
all_stats <- do.call(rbind, lapply(composites, function(comp) rbind(
  bind(lapply(c(vi_stems, raw_stems), function(s) bind(model_table(s, "index", comp)))),
  bind(lapply(leaf_stems, function(s) bind(model_table(s, "lai", comp)))))))
message(nrow(all_stats), " models cross-validated in ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

best_of <- function(s) do.call(rbind, lapply(unique(s$Model.stem), function(st) {
  a <- s[s$Model.stem == st, ]
  if (all(a$Model.call %in% c("identity", "calibrated"))) return(a)
  a <- a[a$Model.call %in% index_forms, ]; a[best_rank(a$RMSE)[1], ] }))
for (comp in unique(all_stats$composite)) {
  s <- all_stats[all_stats$composite == comp, ]
  save_table(fmt_stats(s),          file.path("model_summary_stats", paste0("model_summary_stats_kfold_", comp)))
  save_table(fmt_stats(best_of(s)), file.path("model_summary_stats", paste0("model_summary_stats_kfold_best_", comp)))
}
save_table(cbind(Composite = all_stats$composite, fmt_stats(all_stats)),
           file.path("model_summary_stats", "model_summary_stats_kfold_all_composites"))

#############################################################################
# 5. manuscript models: key statistics, Table 2 and the scatter panel
#############################################################################
# the six models of 5a's scatter panel, each index in its best form by
# cross-validated RMSE, over the closest-in-time composite
panel_index <- c("NDVI", "EVI")
closest <- all_stats[all_stats$composite == sub("^_", "", focus_composite), ]
best_form <- function(stem) { a <- closest[closest$Model.stem == stem & closest$Model.call %in% index_forms, ]
  a$Model.call[best_rank(a$RMSE)[1]] }
panel <- list(
  "NDVI (Sen2Cor)"       = model_table("rawNDVI", "index", keep = TRUE)[[best_form("rawNDVI")]],
  "EVI (Sen2Cor)"        = model_table("rawEVI",  "index", keep = TRUE)[[best_form("rawEVI")]],
  "SL2P 10 m (uncalibrated)" = model_table("leaf_10m", "lai", keep = TRUE)$identity,
  "NDVI (NBART)"         = model_table("NDVI", "index", keep = TRUE)[[best_form("NDVI")]],
  "EVI (NBART)"          = model_table("EVI",  "index", keep = TRUE)[[best_form("EVI")]],
  "SL2P 10 m (calibrated)"   = model_table("leaf_10m", "lai", keep = TRUE)$calibrated)

key_stats <- fmt_stats(bind(panel), model = names(panel))
save_table(key_stats, "validation_key_stats_kfold")
print(key_stats, row.names = FALSE)

# bias by reference LAI class, from each record's mean out-of-fold prediction
lai_bins <- seq(0, ceiling(2 * max(total_LAI$LAI_LX_total)) / 2, by = 0.5)
lai_cls  <- function(d) cut(d$LAI_LX_total, lai_bins, include.lowest = TRUE)
cls_n    <- table(lai_cls(attr(panel[[1]], "pred")))
bias_cls <- t(sapply(panel, function(s) { d <- attr(s, "pred"); sgn(tapply(d$.pred - d$LAI_LX_total, lai_cls(d), mean)) }))
colnames(bias_cls) <- sprintf("%s%s – %s (n = %d)", c("", rep("> ", length(cls_n) - 1)),
                              head(lai_bins, -1), lai_bins[-1], cls_n)
table2 <- cbind(key_stats[, c("Model", "Regression model", "NRMSE", "p")], bias_cls, row.names = NULL)
save_table(table2[order(grepl("SL2P", table2$Model)), ], "Table2_validation_models_kfold",   # SL2P rows last
           top = "LAI bias", top_cols = 4 + seq_along(cls_n))

# reference against mean out-of-fold prediction on one axis range, annotated
# with the cross-validated statistics (mean over the repeats)
lims <- range(unlist(lapply(panel, function(s) { d <- attr(s, "pred"); c(d$.pred, d$LAI_LX_total) })), na.rm = TRUE)
step <- diff(lims) * 0.07
p <- Map(function(s, label) {
  d   <- attr(s, "pred")
  ann <- c(sprintf("n = %d, %d sites", s$n_val, length(unique(d$Site))), sprintf("bias = %+.2f", round(s$bias, 2) + 0),
           sprintf("RMSE = %.2f", s$RMSE), sprintf("italic(r)^2 == '%.2f'", s$r2_val))   # quoted keeps a trailing zero
  form <- if (s$Model.call %in% c("identity", "calibrated")) "" else sprintf(", %s", s$Model.call)
  g <- ggplot(d, aes(x = LAI_LX_total, y = .pred)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = 'grey40') +
    geom_smooth(method = "lm", formula = y ~ x, colour = 'black', alpha = 0.4, linewidth = 0.5) +
    geom_point(aes(colour = MVG), size = 1.4) + scale_colour_manual(values = mvg_cols)
  for (j in seq_along(ann))
    g <- g + annotate("text", x = lims[1], y = lims[2] - (j - 1) * step, label = ann[j],
                      parse = (j == 4), hjust = 0, vjust = 1, size = txt.size / 3.2)
  g + coord_equal(xlim = lims, ylim = lims) +
    labs(x = expression(reference~LAI~(m^2~m^-2)), y = expression(modelled~LAI~(m^2~m^-2)),
         title = paste0(label, form)) +
    theme_bw() +
    theme(plot.title = element_text(size = txt.size, face = "plain"), legend.position = 'bottom',
          legend.title = element_blank(), legend.text = element_text(size = txt.size), panel.grid = element_blank(),
          axis.text = element_text(size = txt.size, colour = 'black'), axis.title = element_text(size = txt.size)) +
    guides(colour = guide_legend(ncol = 3))
}, panel, names(panel))
ggsave(file.path(plot_path, "model_vs_reference_scatter_kfold.tif"),
       Reduce(`+`, p) + plot_layout(ncol = 3, guides = 'collect') & theme(legend.position = 'bottom'),
       width = 21.3, height = 18, units = "cm", dpi = 300)
message("done: ", plot_path)
