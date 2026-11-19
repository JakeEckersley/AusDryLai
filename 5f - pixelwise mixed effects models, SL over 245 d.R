################################################################################
#### 5f - pixelwise mixed effects models, SL over 245 d                     ####
################################################################################
# Script 5c on the growth phases that follow a senescence phase longer than
# sl_min days. SL belongs to growth phases only, so there is no all-phase model:
# the rain figures and their curves come from growth-phase LAI change.
#   y ~ ns(RE1) + ... + ns(RE>=n+1) + SL + Landform + (each term):Landform + (1 | Site)
# R2 is partitioned individually (Shapley, Lai et al. 2022) and uniquely
# (leave one term out), with parametric bootstrap intervals as in partR2.
# Everything is written to plots/<n>_events/SL_subset (SL_subset_no_zero_floor with
# drop_zero). Stages are saved as they finish and skipped when their file exists;
# delete the mem_*.rds files there to refit. Once both runs exist, the last section
# compares them.

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
sl_min   <- 245         # days; growth phases after a longer senescence phase only
drop_zero <- FALSE      # TRUE: leave out growth phases whose trough (LAI at t1) is on the zero floor
fig_path <- file.path('plots', c('one_event', 'two_events', 'three_events')[n_events],
                      if (drop_zero) 'SL_subset_no_zero_floor' else 'SL_subset')
dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

vi_best     <- readRDS(file.path(end_path, 'Optimum_closest_VI_choice.rds'))$stem
mod_sources <- c(vi_best, 'SL2P 10 m calibrated')
src_short   <- setNames(c(vi_best, 'SL2Pcal'), mod_sources)

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
get_dat <- function(s, growth, no_zero = drop_zero) {
  d <- phase_all[phase_all$source == s & !is.na(phase_all$Landform) & phase_all$Landform != 'Bare channel', ]
  if (growth) d <- d[which(d$turn == 1 & (!no_zero | d$LAI_prev > 0)), ]
  na.omit(data.frame(Site = d$Site, Landform = factor(d$Landform, levels = lf_levels), SL = d$length_prev_phase,
                     LAI_change = d$LAI_change, phase_length = as.numeric(d$phase_length),
                     setNames(d[rain_cols], rain_terms)))
}
dat <- lapply(lapply(setNames(mod_sources, mod_sources), get_dat, growth = TRUE), function(d) droplevels(d[d$SL > sl_min, ]))
rm(phase_all); invisible(gc())
for (s in mod_sources) message(src_short[[s]], ': ', nrow(dat[[s]]), ' growth phases with SL over ', sl_min, ' d at ',
                               length(unique(dat[[s]]$Site)), ' sites')

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
pe_models <- list(LAI_change   = list(y = 'LAI_change',   d = dat, x = x_terms),
                  phase_length = list(y = 'phase_length', d = dat, x = x_terms))
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

# every model of pe_models: growth-phase LAI change and growth length
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
ttl     <- c(LAI_change = 'italic(Delta * LAI[t[1] %->% t[2]]) ~ "(growth phases)"', phase_length = 'italic(GL[t[1] %->% t[2]])')
th <- theme_bw() +
  theme(panel.grid = element_blank(), axis.title.y = element_blank(), axis.title.x = element_text(size = 15.6),
        axis.text = element_text(size = 15.6, colour = 'black'), plot.title = element_text(size = 15.6),
        legend.text = element_text(size = 15.6), legend.key.spacing.x = unit(12, 'pt'), legend.position = 'bottom')

# each estimate drawn over its bootstrap interval: dark circle = the term on its own,
# light triangle = with its landform interaction; within a row, top to bottom: NDVI,
# NDVI with interaction, SL2Pcal, SL2Pcal with interaction
# the R2 axis runs to 0.6, or further when an estimate or interval does
r2_panel <- function(y, pt, v = vp, xmax = max(0.6, ceiling(10 * max(v$solo_upper, v$total_upper, na.rm = TRUE)) / 10)) {
  d  <- v[v$response == y & v$partition %in% c('model', pt), ]
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
    coord_flip(ylim = c(-0.05, xmax)) +
    scale_x_continuous(breaks = rev(seq_along(lv)), labels = parse(text = lab_map[lv]), expand = expansion(add = 0.5)) +
    scale_colour_manual(values = cols, breaks = names(cols)[c(1, 3, 2, 4)], name = NULL) +
    scale_shape_manual(values = setNames(c(16, 16, 17, 17), names(cols)), breaks = names(cols)[c(1, 3, 2, 4)], name = NULL) +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3.85))) +
    labs(y = bquote(.(paste(pt, 'variance')) ~ (italic(R)^2)), title = parse(text = ttl[[y]]))
}
# growth-phase models side by side
r2_fig <- function(parts, file, h) {
  p <- unlist(lapply(parts, function(pt) lapply(resps, r2_panel, pt = pt)), recursive = FALSE)
  for (i in seq(2, length(p), 2)) p[[i]] <- p[[i]] + theme(axis.text.y = element_blank())
  ggsave(wrap_plots(p, ncol = 2, guides = 'collect') & theme(legend.position = 'bottom'),
         filename = file.path(fig_path, file), width = 30, height = h, units = 'cm', dpi = 300)
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
          axis.title = element_text(size = txt.size), axis.text = element_text(size = txt.size - 2),
          legend.text = element_text(size = txt.size), legend.title = element_text(size = txt.size))
  if (zero) g <- g + geom_hline(yintercept = 0, linetype = 'dashed')
  g <- g + geom_ribbon(data = pr, inherit.aes = FALSE, alpha = 0.18, fill = 'grey25',
                       aes(x = x, ymin = conf.low, ymax = conf.high, group = group)) +
    geom_line(data = pr, inherit.aes = FALSE, linewidth = 1, aes(x = x, y = predicted, colour = group)) +
    coord_cartesian(xlim = c(min(d[[xvar]]), quantile(d[[xvar]], 0.99)))      # x cut at the 99th percentile
  g + do.call(scale_y_continuous, c(list(name = y_lab), y_scale[[yvar]]))
}
x_labs <- c(lapply(rain_labs, function(l) bquote(italic(.(str2lang(l))) ~ '(mm)')),
            list(SL = bquote(italic(SL) ~ '(days)')))
for (x in x_terms) {
  p <- list()
  for (s in mod_sources)
    p[[paste('lai', s)]] <- dens_panel(dat[[s]], x, 'LAI_change',
                                       bquote(italic(Delta * LAI[t[1] %->% t[2]]) ~ .(sprintf('(%s)', src_short[[s]]))),
                                       pred[[s]]$LAI_change[[x]])
  for (s in mod_sources)
    p[[paste('gl', s)]] <- dens_panel(dat[[s]], x, 'phase_length',
                                      bquote(italic(GL[t[1] %->% t[2]]) ~ .(sprintf('(days, %s)', src_short[[s]]))),
                                      pred[[s]]$phase_length[[x]])
  fig <- wrap_plots(p, nrow = 2, guides = 'collect') &
    theme(legend.position = 'right', legend.key.height = unit(0.9, 'cm'), legend.spacing.y = unit(0.7, 'cm'),
          legend.title = element_text(size = txt.size, margin = margin(b = 8)))
  ggsave(fig, filename = file.path(fig_path, paste0('CMP LAIchange and growthlength vs ', x, '.tiff')),
         width = 33, height = 28, units = 'cm')
}

########################################################################################################################
## sites and seasons with SL over sl_min days
########################################################################################################################
# growth phases after SL > sl_min per source, site, landform and wet season (July to June) of t1, the end of the
# senescence phase: the share of the landform's pixels at that site they cover and their median SL, t0 and t1
ph <- readRDS(file.path(end_path, 'LAI_phase_stats_by_source.rds'))
ph <- ph[which(ph$source %in% mod_sources & ph$turn == 1 & !is.na(ph$Landform) & ph$Landform != 'Bare channel' &
               !is.na(ph$length_prev_phase) & (!drop_zero | ph$LAI_prev > 0)), c('source', 'Site', 'Landform', 'pixel', 'Date', 'phase_length', 'length_prev_phase')]
ph$t1 <- ph$Date - ph$phase_length
ph$t0 <- ph$t1 - ph$length_prev_phase
ph$season <- as.integer(format(ph$t1 - 181, '%Y'))
n_pix <- tapply(ph$pixel, paste(ph$source, ph$Site, ph$Landform), function(p) length(unique(p)))

long <- ph[ph$length_prev_phase > sl_min, ]
by <- long[c('source', 'Site', 'Landform', 'season')]
sl_tab <- aggregate(list(phases = long$pixel), by, length)
sl_tab$pixels <- aggregate(long$pixel, by, function(p) length(unique(p)))$x
sl_tab$landform_pixels <- as.vector(n_pix[paste(sl_tab$source, sl_tab$Site, sl_tab$Landform)])
sl_tab$pct_pixels <- round(100 * sl_tab$pixels / sl_tab$landform_pixels, 1)
sl_tab$SL_median <- aggregate(long$length_prev_phase, by, median)$x
sl_tab$t0_median <- as.Date(aggregate(as.numeric(long$t0), by, median)$x)
sl_tab$t1_median <- as.Date(aggregate(as.numeric(long$t1), by, median)$x)
sl_tab$source <- src_short[sl_tab$source]
sl_tab$season <- sprintf('%d-%02d', sl_tab$season, (sl_tab$season + 1) %% 100)
sl_tab <- sl_tab[order(sl_tab$source, -sl_tab$pct_pixels), ]
for (s in src_short) {
  write.csv(sl_tab[sl_tab$source == s, -1], file.path(fig_path, sprintf('SL_over_%d_d_by_site_landform_and_season_%s.csv', sl_min, s)),
            row.names = FALSE)
  print(sl_tab[sl_tab$source == s & sl_tab$pct_pixels >= 20, -1], row.names = FALSE)
}
for (p in c(50, 20)) for (s in src_short) {
  hi <- sl_tab[sl_tab$pct_pixels >= p & sl_tab$source == s, ]
  message(s, ': SL over ', sl_min, ' d in >= ', p, '% of a landform\'s pixels at ', nrow(unique(hi[c('Site', 'Landform')])),
          ' site landforms (', length(unique(hi$Site)), ' sites, ', sum(hi$phases), ' phases) in ', length(unique(hi$season)),
          ' wet seasons')
}

########################################################################################################################
## all growth phases against those without zero-floor troughs
########################################################################################################################
# needs a run with each drop_zero setting: SL kernel density with the landform curves, one figure per source, and the
# individual R2 of both sources, one row per set of growth phases, on a shared R2 axis
sc_path <- setNames(file.path(dirname(fig_path), c('SL_subset', 'SL_subset_no_zero_floor')), c('all', 'no_zero'))
sc_lab  <- c(all = 'All growth phases', no_zero = 'Zero-floor troughs removed')
if (all(file.exists(file.path(sc_path, 'mem_varpart_outputs.csv')))) {
  phase_all <- readRDS(file.path(end_path, 'LAI_phase_stats_by_source.rds'))
  sc_dat <- lapply(c(all = FALSE, no_zero = TRUE), function(z) lapply(setNames(mod_sources, mod_sources), function(s) {
    d <- get_dat(s, TRUE, z)
    droplevels(d[d$SL > sl_min, ])
  }))
  rm(phase_all); invisible(gc())
  sc_pred <- lapply(sc_path, function(p) readRDS(file.path(p, 'mem_partial_effects_by_source.rds')))

  for (s in mod_sources) {
    y_lab <- list(LAI_change   = bquote(italic(Delta * LAI[t[1] %->% t[2]]) ~ .(sprintf('(%s)', src_short[[s]]))),
                  phase_length = bquote(italic(GL[t[1] %->% t[2]]) ~ .(sprintf('(days, %s)', src_short[[s]]))))
    p <- list()
    for (k in resps) for (z in names(sc_path)) {
      g <- dens_panel(sc_dat[[z]][[s]], 'SL', k, y_lab[[k]], sc_pred[[z]][[s]][[k]]$SL) +
        expand_limits(x = sl_min - 50, y = if (k == 'LAI_change') c(-0.1, 4) else 1)   # density grid past the data edges, so contours close
      p[[paste(k, z)]] <- if (k == 'LAI_change') g + labs(title = sc_lab[[z]]) else g
    }
    ggsave(wrap_plots(p, nrow = 2, guides = 'collect') &
             theme(legend.position = 'right', legend.key.height = unit(0.9, 'cm'), legend.spacing.y = unit(0.7, 'cm'),
                   legend.title = element_text(size = txt.size, margin = margin(b = 8)), plot.title = element_text(size = txt.size)),
           filename = file.path(sc_path[['all']], sprintf('CMP LAIchange and growthlength vs SL, all vs no zero floor, %s.tiff', src_short[[s]])),
           width = 33, height = 28, units = 'cm')
  }

  sc_vp <- lapply(sc_path, function(p) read.csv(file.path(p, 'mem_varpart_outputs.csv')))
  xmax  <- max(0.6, ceiling(10 * max(sapply(sc_vp, function(v) max(v$solo_upper, v$total_upper, na.rm = TRUE)))) / 10)
  p <- list()
  for (z in names(sc_path)) for (k in resps) {
    g <- r2_panel(k, 'individual', sc_vp[[z]], xmax)
    g <- if (k == resps[1]) g + labs(x = sc_lab[[z]]) + theme(axis.title.y = element_text(size = 15.6, angle = 90)) else
      g + theme(axis.text.y = element_blank())
    p[[paste(z, k)]] <- if (z == 'all') g else g + labs(title = NULL)
  }
  ggsave(wrap_plots(p, ncol = 2, guides = 'collect') & theme(legend.position = 'bottom'),
         filename = file.path(sc_path[['all']], 'CMP R2 individual, all vs no zero floor.tiff'),
         width = 30, height = (32 + 4 * n_events) * 2 / 3, units = 'cm', dpi = 300)
}

message('done: ', paste(names(boot), collapse = ', '))
