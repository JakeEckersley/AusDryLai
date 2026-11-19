################################################################################
#### 5c - pixelwise mixed effects models across LAI sources                 ####
################################################################################
# Growth length and LAI change of growth phases against the rain of the lagged
# window, split into its n_events largest events (RE1 ... REn) and the rest
# (RE>=n+1), with senescence length (SL) and landform, for each source in mod_sources:
#   y ~ ns(RE1) + ... + ns(RE>=n+1) + SL + Landform + (each term):Landform + (1 | Site)
# R2 is partitioned individually (Shapley, Lai et al. 2022) and uniquely
# (leave one term out), with parametric bootstrap intervals as in partR2. The
# rain figures refit LAI change on every phase, growth and decline, and that
# model is partitioned the same way.
# Everything is written to plots/one_event, two_events or three_events. Stages
# are saved as they finish and skipped when their file exists; delete the
# mem_*.rds files there to refit.

library(ggplot2)
library(patchwork)
library(lme4)
library(splines)
library(ggeffects)
library(RColorBrewer)
library(future)
library(future.apply)

options(future.globals.maxSize = 5 * 1024^3)

end_path <- 'outputs'
n_events <- 2           # largest events kept apart (1 to 3, as stored by 5b); the rest of the rain is one term
fig_path <- file.path('plots', c('one_event', 'two_events', 'three_events')[n_events])
dir.create(fig_path, showWarnings = FALSE)

vi_best     <- readRDS(file.path(end_path, 'Optimum_closest_VI_choice.rds'))$stem
mod_sources <- c(vi_best, 'SL2P 10 m calibrated')
src_short   <- setNames(c(vi_best, 'SL2Pcal'), mod_sources)
src_lab     <- function(x, f = '%s NBART') sub(paste0('^', vi_best), sprintf(f, vi_best), x)   # in figures; src_short keys the saved bootstraps

lag_days  <- 14         # rainfall lag; 5b stores 0 to 56 d in 7 d steps
booties   <- 50
n_workers <- 20
txt.size  <- 15

rain_terms <- c(paste0('RE', seq_len(n_events)), paste0('REge', n_events + 1))
rain_cols  <- setNames(paste0(c('LRE', if (n_events > 1) paste0('RE', 2:n_events), paste0('MRE', if (n_events > 1) n_events)),
                              '_', lag_days, '_Dlag'), rain_terms)
rain_labs  <- setNames(c(sprintf('RE[%d]', seq_len(n_events)), sprintf('RE[phantom() >= %d]', n_events + 1)), rain_terms)
x_terms    <- c(rain_terms, 'SL')
resps      <- c('LAI_change', 'phase_length')

pred_file <- file.path(fig_path, 'mem_partial_effects_by_source.rds')
boot_file <- file.path(fig_path, 'mem_results_by_source.rds')

phase_all <- readRDS(file.path(end_path, 'LAI_phase_stats_by_source.rds'))
stopifnot(all(mod_sources %in% phase_all$source))
lf_levels <- sort(setdiff(unique(phase_all$Landform[phase_all$source %in% mod_sources]), c(NA, 'Bare channel')))
get_dat <- function(s, growth) {
  d <- phase_all[phase_all$source == s & !is.na(phase_all$Landform) & phase_all$Landform != 'Bare channel', ]
  if (growth) d <- d[which(d$turn == 1), ]
  na.omit(data.frame(Site = d$Site, Landform = factor(d$Landform, levels = lf_levels), SL = d$length_prev_phase,
                     LAI_change = d$LAI_change, phase_length = as.numeric(d$phase_length),
                     setNames(d[rain_cols], rain_terms)))
}
dat     <- lapply(setNames(mod_sources, mod_sources), get_dat, growth = TRUE)
dat_all <- lapply(setNames(mod_sources, mod_sources), get_dat, growth = FALSE)
rm(phase_all); invisible(gc())

# Pearson and Spearman correlation between the rain terms over the growth phases
print(do.call(rbind, lapply(mod_sources, function(s) {
  i <- which(upper.tri(diag(length(rain_terms))), arr.ind = TRUE)
  data.frame(source = src_short[[s]], n = nrow(dat[[s]]), term_1 = rain_terms[i[, 1]], term_2 = rain_terms[i[, 2]],
             pearson  = round(cor(dat[[s]][rain_terms])[i], 2),
             spearman = round(cor(dat[[s]][rain_terms], method = 'spearman')[i], 2))
})), row.names = FALSE)

# rain splines have one knot, at the median of the non-zero values
term_code <- function(d) {
  m <- c(setNames(sprintf('ns(%s, knots = %.6g)', rain_terms,
                          sapply(rain_terms, function(v) median(d[[v]][d[[v]] > 0]))), rain_terms), SL = 'SL')
  c(m, Landform = 'Landform', setNames(paste0(m, ':Landform'), paste0(names(m), ':LF')))
}
mem_form <- function(y, tc, set = names(tc)) paste(y, '~', paste(tc[set], collapse = ' + '), '+ (1 | Site)')

########################################################################################################################
## full models and partial effects
########################################################################################################################
# LAI change against rain uses every phase; against SL, and growth length, growth phases
pe_models <- list(LAI_change_all = list(y = 'LAI_change',   d = dat_all, x = rain_terms),
                  LAI_change     = list(y = 'LAI_change',   d = dat,     x = 'SL'),
                  phase_length   = list(y = 'phase_length', d = dat,     x = x_terms))
if (!file.exists(pred_file)) {
  pred <- list()
  for (s in mod_sources) for (k in names(pe_models)) {
    d <- pe_models[[k]]$d[[s]]
    m <- lmer(as.formula(mem_form(pe_models[[k]]$y, term_code(d))), data = d)
    for (x in pe_models[[k]]$x) {
      q  <- tapply(d[[x]], d$Landform, quantile, 0.99)        # curves stop at their landform's 99th percentile
      at <- paste(signif(seq(min(d[[x]]), max(q), length.out = 100), 6), collapse = ',')
      pr <- as.data.frame(predict_response(m, terms = c(sprintf('%s [%s]', x, at), 'Landform'), type = 'fixed'))
      pred[[s]][[k]][[x]] <- pr[pr$x <= q[as.character(pr$group)], c('x', 'predicted', 'conf.low', 'conf.high', 'group')]
    }
    message(src_short[[s]], ' ', k, ' fitted')
  }
  saveRDS(pred, pred_file)
  rm(m); invisible(gc())
}
pred <- readRDS(pred_file)

########################################################################################################################
## R2 of every submodel, with parametric bootstraps
########################################################################################################################
# players: Landform, each term and each term's landform interaction; an
# interaction only enters a model holding its term and Landform
items <- c('Landform', x_terms, paste0(x_terms, ':LF'))
bit   <- 2^(seq_along(items) - 1)
in_m  <- function(m) items[bitwAnd(m, bit) > 0]
valid <- Filter(function(m) { u <- in_m(m); i <- sub(':LF$', '', grep(':LF$', u, value = TRUE))
                              all(c(i, if (length(i)) 'Landform') %in% u) }, 0:(2^length(items) - 1))

# every model of pe_models: growth-phase LAI change and growth length, and LAI change over every phase
boot <- if (file.exists(boot_file)) readRDS(boot_file) else list()
for (s in mod_sources) for (k in names(pe_models)) {
  p <- paste(src_short[[s]], k)
  if (!is.null(boot[[p]])) next
  y  <- pe_models[[k]]$y; d <- pe_models[[k]]$d[[s]]; tc <- term_code(d)
  m  <- lmer(as.formula(mem_form(y, tc)), data = d)
  mu <- as.vector(getME(m, 'X') %*% fixef(m)); su <- attr(VarCorr(m)$Site, 'stddev'); se <- sigma(m)
  site  <- as.integer(factor(d$Site)); rm(m)
  forms <- setNames(sapply(valid[-1], function(k) mem_form(y, tc, in_m(k))), valid[-1])
  t0 <- Sys.time()
  plan(multisession, workers = n_workers)
  res <- future_lapply(forms[order(-nchar(forms))], function(f) {
    r2 <- function(m) { vf <- var(as.vector(getME(m, 'X') %*% fixef(m))); vr <- VarCorr(m)$Site[1]
                        c(vf, vf + vr) / (vf + vr + sigma(m)^2) }
    w <- character(0)
    m <- withCallingHandlers(lmer(as.formula(f), data = d),
                             warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart('muffleWarning') })
    out <- matrix(r2(m), 1)
    for (b in seq_len(booties)) {
      set.seed(b)        # the same simulated responses for every submodel
      out <- rbind(out, r2(suppressWarnings(refit(m, newresp = mu + rnorm(max(site), 0, su)[site] + rnorm(length(mu), 0, se)))))
    }
    list(r2 = out, warning = paste(unique(w), collapse = '; '))
  }, future.seed = TRUE, future.chunk.size = 1, future.packages = c('lme4', 'splines'))
  plan(sequential)
  res <- res[names(forms)]
  boot[[p]] <- list(source = src_short[[s]], response = k, R2_M = cbind(`0` = 0, sapply(res, function(r) r$r2[, 1])),
                    R2_C = res[[as.character(max(valid))]]$r2[, 2],
                    warning = sapply(res, `[[`, 'warning'), booties = booties)
  saveRDS(boot, boot_file)
  message(p, ': ', length(res), ' models, ', round(difftime(Sys.time(), t0, units = 'mins')), ' min')
}

########################################################################################################################
## variance partitions
########################################################################################################################
# Shapley weights: valid orders reaching each model (f_up) and completing it (g_dn)
pc   <- sapply(valid, function(m) length(in_m(m)))
f_up <- g_dn <- setNames(numeric(length(valid)), valid)
for (m in valid[order(pc)])  f_up[[as.character(m)]] <- if (m == 0) 1 else {
  q <- m - bit[bitwAnd(m, bit) > 0];  sum(f_up[as.character(q[q %in% valid])]) }
for (m in valid[order(-pc)]) g_dn[[as.character(m)]] <- if (m == max(valid)) 1 else {
  q <- m + bit[bitwAnd(m, bit) == 0]; sum(g_dn[as.character(q[q %in% valid])]) }
W <- matrix(0, length(items), length(valid), dimnames = list(items, valid))
for (i in seq_along(items)) for (S in valid[bitwAnd(valid, bit[i]) == 0 & (valid + bit[i]) %in% valid]) {
  w <- f_up[[as.character(S)]] * g_dn[[as.character(S + bit[i])]] / f_up[[as.character(max(valid))]]
  W[i, as.character(S + bit[i])] <- W[i, as.character(S + bit[i])] + w
  W[i, as.character(S)]          <- W[i, as.character(S)] - w
}

part <- function(v, cond) {
  R    <- function(u) v[[as.character(sum(bit[items %in% u]))]]
  less <- function(u) setdiff(items, u)
  lf   <- function(x) paste0(x, ':LF')
  full <- R(items)
  phi  <- setNames(as.vector(W %*% v), items)
  u_solo <- c(sapply(x_terms, function(x) R(less(lf(x))) - R(less(c(x, lf(x))))),
              R(less(lf(x_terms))) - R(less(c(lf(x_terms), 'Landform'))),
              R(less(lf(rain_terms))) - R(less(c(rain_terms, lf(rain_terms)))))
  u_lf   <- c(sapply(x_terms, function(x) full - R(less(lf(x)))), full - R(less(lf(x_terms))), full - R(less(lf(rain_terms))))
  i_solo <- c(phi[x_terms], phi[['Landform']], sum(phi[rain_terms]))
  i_lf   <- c(phi[lf(x_terms)], sum(phi[lf(x_terms)]), sum(phi[lf(rain_terms)]))
  terms  <- c(x_terms, 'Landform', 'Rainfall terms')
  rbind(data.frame(partition = 'model', term = c('Full model', 'Fixed effects', 'Site effect'),
                   solo = c(cond, full, cond - full), total = NA),
        data.frame(partition = 'individual', term = terms, solo = i_solo, total = i_solo + i_lf),
        data.frame(partition = 'unique', term = terms, solo = u_solo, total = u_solo + u_lf))
}

vp <- do.call(rbind, lapply(boot, function(b) {
  reps <- lapply(seq_len(nrow(b$R2_M)), function(i) part(b$R2_M[i, ], b$R2_C[i]))
  ci   <- function(col, q) apply(sapply(reps[-1], `[[`, col), 1, quantile, q, na.rm = TRUE)
  cbind(source = b$source, response = b$response, reps[[1]],
        solo_lower = ci('solo', 0.025), solo_upper = ci('solo', 0.975),
        total_lower = ci('total', 0.025), total_upper = ci('total', 0.975))
}))
write.csv(vp, file.path(fig_path, 'mem_varpart_outputs.csv'), row.names = FALSE, na = '')
fits_warn <- unlist(lapply(boot, `[[`, 'warning'))
message(sum(grepl('converge', fits_warn)), ' of ', length(fits_warn), ' submodel fits gave a convergence warning')

########################################################################################################################
## side by side comparison across sources
########################################################################################################################
cols    <- c(setNames(c('#1b6ca8', '#c1443c'), src_short),                               # on its own
             setNames(c('#8fbde0', '#eba49e'), paste(src_short, 'with landform interaction')))
lab_map <- c(`Full model` = 'Full~model', `Fixed effects` = 'Fixed~effects', `Site effect` = 'Site~effect',
             rain_labs, SL = 'SL', Landform = 'Landform', `Rainfall terms` = 'Rainfall~terms')
ttl     <- c(LAI_change = 'italic(Delta * LAI[t[1] %->% t[2]]) ~ "(growth phases)"', phase_length = 'italic(GL[t[1] %->% t[2]])',
             LAI_change_all = 'italic(Delta * LAI[t[1] %->% t[2]]) ~ "(all phases)"')
th <- theme_bw() +
  theme(panel.grid = element_blank(), axis.title.y = element_blank(), axis.title.x = element_text(size = 15.6),
        axis.text = element_text(size = 15.6, colour = 'black'), plot.title = element_text(size = 15.6),
        legend.text = element_text(size = 15.6), legend.key.spacing.x = unit(12, 'pt'), legend.position = 'bottom')

# each estimate drawn over its bootstrap interval: dark circle = the term on its own,
# light triangle = with its landform interaction; within a row, top to bottom: NDVI,
# NDVI with interaction, SL2Pcal, SL2Pcal with interaction
r2_panel <- function(y, pt) {
  d  <- vp[vp$response == y & vp$partition %in% c('model', pt), ]
  lv <- unique(d$term)
  e  <- rbind(data.frame(d[c('term', 'source')], est = d$solo,  lo = d$solo_lower,  hi = d$solo_upper,  lf = FALSE),
              data.frame(d[c('term', 'source')], est = d$total, lo = d$total_lower, hi = d$total_upper, lf = TRUE))
  e  <- e[is.finite(e$est), ]
  e$key <- ifelse(e$lf, paste(e$source, 'with landform interaction'), e$source)
  e$slot <- 2 * match(e$source, src_short) - 1 + e$lf
  k <- ave(e$slot, e$term, FUN = rank); n <- ave(e$slot, e$term, FUN = length)
  e$x <- length(lv) + 1 - match(e$term, lv) + ((n + 1) / 2 - k) * 0.22
  ggplot(e, aes(x = x, colour = key, shape = key)) + th +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.5, na.rm = TRUE, show.legend = FALSE) +
    geom_point(aes(y = est), size = 2.86) +
    coord_flip(ylim = c(-0.05, 0.6)) +
    scale_x_continuous(breaks = rev(seq_along(lv)), labels = parse(text = lab_map[lv]), expand = expansion(add = 0.5)) +
    scale_colour_manual(values = cols, breaks = names(cols)[c(1, 3, 2, 4)], labels = function(b) src_lab(b, '%s (NBART)'),
                        name = NULL) +
    scale_shape_manual(values = setNames(c(16, 16, 17, 17), names(cols)), breaks = names(cols)[c(1, 3, 2, 4)],
                       labels = function(b) src_lab(b, '%s (NBART)'), name = NULL) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3.85))) +
    labs(y = bquote(.(paste(pt, 'variance')) ~ (italic(R)^2)), title = parse(text = ttl[[y]]))
}
# growth-phase models side by side, and LAI change over every phase on its own
r2_fig <- function(parts, file, h) {
  p <- unlist(lapply(parts, function(pt) lapply(resps, r2_panel, pt = pt)), recursive = FALSE)
  for (i in seq(2, length(p), 2)) p[[i]] <- p[[i]] + theme(axis.text.y = element_blank())
  ggsave(wrap_plots(p, ncol = 2, guides = 'collect') & theme(legend.position = 'bottom'),
         filename = file.path(fig_path, file), width = 30, height = h, units = 'cm', dpi = 300)
  ggsave(wrap_plots(lapply(parts, r2_panel, y = 'LAI_change_all'), ncol = 1, guides = 'collect') & theme(legend.position = 'bottom') &
           guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3.85))),   # one key per line to fit 18 cm
         filename = file.path(fig_path, sub('(\\.tiff?)$', ' (all phases)\\1', file)), width = 18, height = h, units = 'cm', dpi = 300)
}
r2_fig(c('individual', 'unique'), 'CMP R2 overlay by source.tiff', (32 + 4 * n_events) * 2 / 3)
r2_fig('unique', 'CMP R2 unique by source.tif', (32 + 4 * n_events) / 3 + 3)   # the main figure

# kernel density (contours at log steps from 0.1% to 100% of each panel's peak, so one
# legend serves all) and the partial effect per landform, other terms at their means and the site effect at
# zero. LAI change takes a signed log axis (linear within about 0.05 LAI of zero), growth
# length a log axis
y_scale <- list(LAI_change   = list(transform = scales::pseudo_log_trans(sigma = 0.05, base = 10),
                                    breaks = c(-4, -2, -1, -0.5, -0.1, 0, 0.1, 0.5, 1, 2, 4)),
                phase_length = list(transform = 'log10', breaks = c(3, 10, 30, 100, 300), limits = c(NA, 350)))
dens_panel <- function(d, xvar, yvar, y_lab, pr, zero = FALSE) {
  g <- ggplot(d, aes(x = .data[[xvar]], y = .data[[yvar]])) + theme_bw() +
    geom_point(shape = 21, size = 3, stroke = 1, colour = 'grey80', alpha = 0.1, show.legend = FALSE) +
    stat_density_2d(geom = 'polygon', contour = TRUE, contour_var = 'ndensity', breaks = 10^seq(-3, 0, by = 0.25),
                    aes(fill = after_stat(level)), colour = NA) +
    scale_fill_gradientn(name = 'Kernel\ndensity', colors = alpha(brewer.pal(9, 'YlOrRd')[3:9], alpha = .5), transform = 'log10',
                         limits = c(1e-3, 1), breaks = 10^(-3:0), labels = scales::label_log(), oob = scales::squish) +
    scale_colour_discrete(name = 'Landform', limits = lf_levels) +
    scale_x_continuous(name = x_labs[[xvar]]) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          axis.title = element_text(size = 1.1 * txt.size), axis.text = element_text(size = 1.1 * (txt.size - 2)),
          legend.text = element_text(size = 1.1 * txt.size), legend.title = element_text(size = 1.1 * txt.size))
  if (zero) g <- g + geom_hline(yintercept = 0, linetype = 'dashed')
  g <- g + geom_ribbon(data = pr, inherit.aes = FALSE, alpha = 0.18, fill = 'grey25',
                       aes(x = x, ymin = conf.low, ymax = conf.high, group = group)) +
    geom_line(data = pr, inherit.aes = FALSE, linewidth = 1, aes(x = x, y = predicted, colour = group)) +
    coord_cartesian(xlim = c(min(d[[xvar]]), quantile(d[[xvar]], 0.99)))      # x cut at the 99th percentile
  g + do.call(scale_y_continuous, c(list(name = y_lab), y_scale[[yvar]]))
}
x_labs <- c(lapply(rain_labs, function(l) bquote(italic(.(str2lang(l))) ~ '(mm)')),
            list(SL = bquote(italic(SL) ~ '(days)')))
leg_th <- theme(legend.position = 'bottom', legend.title.position = 'top', legend.key.height = unit(0.6, 'cm'),
                legend.spacing.x = unit(1.5, 'cm'), legend.title = element_text(size = 1.1 * txt.size, margin = margin(b = 8)))
leg_gd <- guides(colour = guide_legend(nrow = 2),
                 fill = guide_colourbar(title = 'Kernel density', theme = theme(legend.key.width = unit(6, 'cm'))))
figs <- list()
for (x in x_terms) {
  k <- if (x == 'SL') 'LAI_change' else 'LAI_change_all'
  p <- list()
  for (s in mod_sources)
    p[[paste('lai', s)]] <- dens_panel(pe_models[[k]]$d[[s]], x, 'LAI_change',
                                       bquote(italic(Delta * LAI[t[1] %->% t[2]]) ~ .(sprintf('(%s)', src_lab(src_short[[s]])))),
                                       pred[[s]][[k]][[x]], zero = x != 'SL')
  for (s in mod_sources)
    p[[paste('gl', s)]] <- dens_panel(dat[[s]], x, 'phase_length',
                                      bquote(italic(GL[t[1] %->% t[2]]) ~ .(sprintf('(days, %s)', src_lab(src_short[[s]])))),
                                      pred[[s]]$phase_length[[x]])
  figs[[x]] <- p
  fig <- wrap_plots(p, nrow = 2, guides = 'collect') & leg_th & leg_gd
  ggsave(fig, filename = file.path(fig_path, paste0('CMP LAIchange', if (x != 'SL') ' (all phases)',
                                                    ' and growthlength vs ', x, '.tiff')),
         width = 33, height = 28, units = 'cm')
}
# the largest event and the remaining rain together
ggsave(wrap_plots(c(figs[[1]], figs[[n_events + 1]]), ncol = 2, guides = 'collect') & leg_th & leg_gd,
       filename = file.path(fig_path, sprintf('CMP LAIchange (all phases) and growthlength vs %s and %s.tiff',
                                              rain_terms[1], rain_terms[n_events + 1])),
       width = 33, height = 46, units = 'cm')

# rain at which the all-phase LAI change curve first rises through zero, by linear
# interpolation; the 95% band edges bracket it (NA: no crossing below the 99th percentile)
zero_cross <- function(x, y) {
  o <- order(x); x <- x[o]; y <- y[o]
  i <- which(y[-1] >= 0 & y[-length(y)] < 0)[1]
  if (is.na(i)) return(NA_real_)
  x[i] + (0 - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
}
crossing <- do.call(rbind, lapply(mod_sources, function(s) do.call(rbind, lapply(rain_terms, function(x) {
  pr <- pred[[s]]$LAI_change_all[[x]]
  do.call(rbind, lapply(levels(pr$group), function(l) {
    q <- pr[pr$group == l, ]
    data.frame(source = src_short[[s]], term = x, Landform = l, rain_mm = zero_cross(q$x, q$predicted),
               upper_limit_crosses = zero_cross(q$x, q$conf.high), lower_limit_crosses = zero_cross(q$x, q$conf.low))
  }))
}))))
write.csv(crossing, file.path(fig_path, 'RE_zero_crossing_by_landform.csv'), row.names = FALSE, na = '')

message('done: ', paste(names(boot), collapse = ', '))
