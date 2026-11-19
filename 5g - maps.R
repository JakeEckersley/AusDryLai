################################################################################
#### 5g - maps                                                              ####
################################################################################
# Every map, over Esri World Imagery through maptiles (no API key; the credit
# line is printed for the captions):
#   1. Fig. 2: two example sites, KC (unconfined riparian) and HDWW (confined
#      riparian), with the NVIS major vegetation groups, landform and the PCA
#      target points
#   2. vegetation composition of the four test sites, and the area of every
#      class at every site
#   3. supplementary: the Fig. 6 test sites with the Fortescue catchments
#   4. modelled LAI at the test sites on the two campaign dates

library(ggplot2)
library(patchwork)
library(terra)
suppressMessages({library(sf); library(maptiles); library(tidyterra)})

plot_path  <- "plots"
txt.size   <- 10
test_sites <- c('MCK', 'CW', 'YU', 'DR')
MVG_shp    <- "data/Vegetation NVIS classification 15-3-23.shp"
rip_shp    <- "data/riparian vegetation boundary edited 2025.shp"
pca_shp    <- "data/LAI_target_locations_for_mapping.shp"
catch_shp  <- "data/Fortescue catchments.shp"
S2_path    <- "D:/PilbaraLaiProject/FliteredTimeseries"   # filtered VI stacks
leaf_path  <- "D:/PilbaraLaiProject/FilteredLeaf"         # QC-filtered SL2P stacks
end_path   <- "outputs"
min_cover  <- 0.95                                       # share of a site's pixels an LAI map image must cover
lai_fs     <- 1.2 * txt.size                             # LAI map text size
tile_cache <- file.path(tempdir(), "maptiles")
dir.create(tile_cache, showWarnings = FALSE, recursive = TRUE)

# the MVG polygons, and one colour per class across the study (R/mvg_palette.R, shared with 5a)
source("R/mvg_palette.R")
veg <- vect(MVG_shp)
veg$MVG_Nam[veg$MVG_Nam == "Naturally bare - sand, rock, claypan, mudflat"] <- "Naturally bare"
mvg_lv   <- sort(unique(veg$MVG_Nam))
stopifnot(all(mvg_lv %in% names(mvg_cols)))

# Esri tiles stop at zoom 17 here: at 18 they come back as a blank "map data not
# yet available" tile
esri <- function(x, zoom = 17) {
  t <- try(get_tiles(x, provider = "Esri.WorldImagery", crop = TRUE, zoom = zoom,
                     cachedir = tile_cache, verbose = FALSE), silent = TRUE)
  if (inherits(t, "try-error")) { message("basemap tiles could not be fetched - drawn without imagery"); NULL } else t
}

#############################################################################
# 1. Fig. 2: two example sites
#############################################################################
# vegetation classes (translucent), landform outlines and the PCA targets. The
# landform layer labels HDWW's channel 'Unconfined riparian'; the figure shows it
# as 'Confined riparian', as the published version did.
sf_use_s2(FALSE)
veg_sf <- st_make_valid(st_transform(st_as_sf(veg), 4326))
rip    <- st_make_valid(st_transform(st_read(rip_shp, quiet = TRUE), 4326))
rip    <- suppressWarnings(st_make_valid(st_intersection(rip, st_union(veg_sf))))
rip    <- rip[rip$Landform != "Bare channel", ]
# the NVIS polygons leave slivers between them, which the clip turns into holes
# drawn as stray rings inside the outlines; keep each polygon's outer ring only
st_geometry(rip) <- st_sfc(lapply(st_geometry(rip), function(g) {
  if (inherits(g, "POLYGON")) st_polygon(g[1])
  else if (inherits(g, "MULTIPOLYGON")) st_multipolygon(lapply(g, function(p) p[1]))
  else g }), crs = st_crs(rip))
rip$Landform[rip$Site == "HDWW" & rip$Landform == "Unconfined riparian"] <- "Confined riparian"
pca    <- st_transform(st_read(pca_shp, quiet = TRUE), 4326)
lf_lt  <- c(`Confined riparian` = "solid", Floodplain = "dotted", Hillslope = "dotdash", `Unconfined riparian` = "dashed")

site_panel <- function(s) {
  w   <- veg_sf[veg_sf$Site == s, ]
  bb  <- st_bbox(w); pad <- 0.15 * max(diff(bb[c(1, 3)]), diff(bb[c(2, 4)]))
  bb[1:2] <- bb[1:2] - pad; bb[3:4] <- bb[3:4] + pad
  t   <- esri(st_as_sfc(bb))
  g   <- ggplot()
  if (!is.null(t)) g <- g + geom_spatraster_rgb(data = t, maxcell = 1e7)
  # show.legend = TRUE draws every key, so both panels share one legend
  g + geom_sf(data = w, aes(fill = MVG_Nam), colour = NA, alpha = 0.45, show.legend = TRUE) +
    geom_sf(data = rip[rip$Site == s, ], aes(linetype = Landform), colour = "black", fill = NA, linewidth = 1,
            show.legend = TRUE) +
    geom_sf(data = st_crop(pca, bb), colour = "yellow2", size = 0.5, show.legend = FALSE) +
    coord_sf(xlim = bb[c(1, 3)], ylim = bb[c(2, 4)], expand = FALSE) +
    scale_fill_manual(values = mvg_cols, limits = mvg_lv[mvg_lv %in% veg_sf$MVG_Nam[veg_sf$Site %in% c("KC", "HDWW")]]) +
    scale_linetype_manual(values = lf_lt, limits = names(lf_lt)) +
    labs(x = "Longitude (degrees)", y = "Latitude (degrees)") +
    theme_bw() +
    theme(axis.text = element_text(colour = "black", size = 15), text = element_text(size = 15),
          axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
}
fig2 <- site_panel("KC") + site_panel("HDWW") + theme(axis.title.y = element_blank()) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 15),
        legend.background = element_rect(colour = "black", fill = "white", linetype = "solid")) &
  guides(fill = guide_legend(nrow = 4, override.aes = list(linewidth = 0.5)),
         linetype = guide_legend(nrow = 4, override.aes = list(linewidth = 0.75)))
ggsave(file.path(plot_path, "site_examples_map.tif"), fig2, width = 30, height = 20, units = "cm", dpi = 600)

#############################################################################
# 2. vegetation composition of the test sites
#############################################################################
# one feature per class per site
site_polys <- lapply(test_sites, function(s) {
  w <- veg[veg$Site == s, ]
  if (!nrow(w)) return(NULL)
  a <- terra::aggregate(w, by = "MVG_Nam")
  a <- st_transform(st_make_valid(st_as_sf(a)), 3857)
  names(a)[names(a) == "MVG_Nam"] <- "MVG"
  a$MVG <- factor(a$MVG, levels = mvg_lv)
  a
})
names(site_polys) <- test_sites

site_ha <- vapply(test_sites, function(s)
  sum(terra::expanse(veg[veg$Site == s, ], unit = "ha")), numeric(1))

# only the classes these four sites hold, so the legend has no empty keys
comp_lv <- sort(unique(unlist(lapply(site_polys, function(w)
  if (is.null(w)) NULL else as.character(w$MVG)))))
site_polys <- lapply(site_polys, function(w) {
  if (!is.null(w)) w$MVG <- factor(as.character(w$MVG), levels = comp_lv)
  w
})

site_tiles <- lapply(site_polys, function(w) {
  if (is.null(w)) return(NULL)
  esri(st_bbox(st_buffer(st_as_sfc(st_bbox(w)), 0.06 * max(diff(st_bbox(w)[c(1, 3)]),
                                                            diff(st_bbox(w)[c(2, 4)])))))
})

# translucent fill over the imagery
comp_panel <- function(s, y_right = FALSE) {
  w <- site_polys[[s]]; if (is.null(w)) return(NULL)
  ha <- site_ha[[s]]
  p <- ggplot()
  if (!is.null(site_tiles[[s]]))
    p <- p + tidyterra::geom_spatraster_rgb(data = site_tiles[[s]], maxcell = 4e6)
  p <- p +
    geom_sf(data = w, aes(fill = MVG), colour = "white", linewidth = 0.3,
            alpha = 0.55) +
    scale_fill_manual(values = mvg_cols, limits = comp_lv, guide = "none")

  p + coord_sf(expand = FALSE, label_graticule = if (y_right) "SE" else "SW") +
    scale_x_continuous(n.breaks = 3) + scale_y_continuous(n.breaks = 4) +
    labs(x = NULL, y = NULL, title = sprintf("%s  (%.1f ha)", s, ha)) +
    theme_bw() +
    theme(plot.title = element_text(size = txt.size, face = "plain"),
          axis.text = element_text(size = txt.size - 2, colour = "black"),
          axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank(),
          plot.margin = margin(3, if (y_right) 3 else 1.5, 3,
                               if (y_right) 1.5 else 3, "mm"))
}

# build the legend
legend_grob <- function() {
  d <- do.call(rbind, lapply(site_polys, function(w) w[, "MVG"]))
  k <- ggplot(d) +
    geom_sf(aes(fill = MVG), alpha = 0.55, colour = "white") +
    scale_fill_manual(values = mvg_cols, limits = comp_lv,
                      guide = guide_legend(nrow = 3)) +
    theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = txt.size))
  g <- ggplot_gtable(ggplot_build(k))
  g$grobs[[which(grepl("guide-box-bottom", g$layout$name))[1]]]
}

maps <- Filter(Negate(is.null), lapply(seq_along(test_sites), function(i)
  comp_panel(test_sites[i], y_right = i %% 2 == 0)))   # ncol = 2, so evens sit right
if (length(maps)) {
  plt <- patchwork::wrap_plots(maps, ncol = 2) /
    patchwork::wrap_elements(full = legend_grob()) +
    patchwork::plot_layout(heights = c(20, 3))
  ggsave(file.path(plot_path, "site_vegetation_composition.tif"),
         plt, width = 26, height = 26, units = "cm", dpi = 300)
  message("composition map written")
}

# area of every class at every site (from polygons)
veg_area <- local({
  a <- expanse(veg, unit = "ha")
  d <- aggregate(list(area_ha = a),
                 by = list(Site = veg$Site, MVG = veg$MVG_Nam), FUN = sum)
  d$area_ha <- round(d$area_ha, 3)
  d[order(d$Site, -d$area_ha), ]
})
veg_area$site_total_ha <- round(
  ave(veg_area$area_ha, veg_area$Site, FUN = sum), 3)
veg_area$pct_of_site <- round(100 * veg_area$area_ha / veg_area$site_total_ha, 1)
write.csv(veg_area, file.path(plot_path, "site_vegetation_area.csv"), row.names = FALSE)

message("vegetation area: ", nrow(veg_area), " site x class combinations, ",
        length(unique(veg_area$Site)), " sites, ",
        length(unique(veg_area$MVG)), " classes")

#############################################################################
# 3. supplementary map of the Fig. 6 test sites
#############################################################################
# a site's point is the centre of its NVIS polygons; the frame is square around
# the catchments
catch    <- st_transform(st_read(catch_shp, quiet = TRUE), 4326)
site_pts <- local({
  v <- project(veg, "EPSG:4326")
  d <- do.call(rbind, lapply(test_sites, function(s) { e <- as.vector(ext(v[v$Site == s, ]))
    data.frame(Site = s, x = mean(e[1:2]), y = mean(e[3:4])) }))
  st_as_sf(d, coords = c("x", "y"), crs = 4326)
})
# labels to the right of each point; MCK and YU are ~15 km apart, so one goes below and one above
lab_xy <- data.frame(Site = site_pts$Site, x = st_coordinates(site_pts)[, 1] + 0.35,
                     y = st_coordinates(site_pts)[, 2] + c(MCK = -0.12, YU = 0.12)[site_pts$Site])
lab_xy$y[is.na(lab_xy$y)] <- st_coordinates(site_pts)[is.na(lab_xy$y), 2]
bb    <- st_bbox(catch); ctr <- c(mean(bb[c(1, 3)]), mean(bb[c(2, 4)]))
half  <- max(diff(bb[c(1, 3)]), diff(bb[c(2, 4)])) / 2 + 0.4
frame <- st_bbox(c(xmin = ctr[1] - half, ymin = ctr[2] - half, xmax = ctr[1] + half, ymax = ctr[2] + half), crs = st_crs(4326))
tiles <- esri(st_as_sfc(frame), zoom = 9)
site_map <- ggplot()
if (!is.null(tiles)) site_map <- site_map + geom_spatraster_rgb(data = tiles, maxcell = 1e7)
site_map <- site_map +
  geom_sf(data = catch, colour = "black", fill = NA) +
  geom_sf(data = site_pts, colour = "red", size = 1.5) +
  geom_text(data = lab_xy, aes(x = x, y = y, label = Site), size = 4, colour = "black") +
  coord_sf(xlim = frame[c(1, 3)], ylim = frame[c(2, 4)], expand = FALSE) +
  labs(x = "Longitude (degrees)", y = "Latitude (degrees)") +
  theme_bw() +
  theme(axis.text = element_text(size = 15, colour = "black"), text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(plot_path, "SuppTimeseriesMap.tif"), site_map, width = 15, height = 15, units = "cm", dpi = 600)

#############################################################################
# 4. modelled LAI at the test sites
#############################################################################
# one row per test site: the study-area outline over imagery, then NDVI and calibrated
# SL2P LAI (as in script 5b) over the whole image, on the image nearest each campaign
# that covers min_cover of the site in both. MCK was sampled in October 2021 only; its
# second image is matched to the other sites' median March 2022 date. Blank pixels have no
# LAI on that date (QC of scripts 3b / 3c, or outside the product). Drawn on the image
# grid (UTM).
paired <- readRDS(file.path(end_path, "paired_ref_and_sat_dataset.rds"))
camp <- do.call(rbind, lapply(test_sites, function(s) sort(unique(as.Date(paired$Date[paired$Site == s])))[1:2]))
camp[is.na(camp[, 2]), 2] <- median(camp[, 2], na.rm = TRUE)
rownames(camp) <- test_sites

paired <- paired[!(paired$n_records_understorey < 8 | paired$n_records_canopy < 8 | paired$n_pixels < 8), ]
vi_choice <- readRDS(file.path(end_path, "Optimum_closest_VI_choice.rds"))
vi_mod <- readRDS(file.path(end_path, "Optimum_closest_VI_lin_reg_LAI_mod.rds"))[[vi_choice$key]]
leaf_cal <- lm(LAI_LX_total ~ leaf_10m_closest, data = paired[!paired$Site %in% test_sites, ])

to_lai <- function(r, mod, var) {
  values(r) <- pmax(predict(mod, setNames(data.frame(values(r)[, 1]), var)), 0)
  r
}

lai_maps <- lapply(setNames(test_sites, test_sites), function(s) {
  vi <- rast(file.path(S2_path, s, paste0(vi_choice$stem, "-ts.tif")))
  sl <- rast(file.path(leaf_path, s, "leaf_10m-ts.tif"))
  v <- project(veg[veg$Site == s, ], crs(vi))
  m <- rasterize(v, vi[[1]])
  cover <- function(r) setNames(global(mask(!is.na(r), m), "mean", na.rm = TRUE)[, 1], substr(names(r), 1, 10))
  cv <- cover(vi); cs <- cover(sl)
  ok <- intersect(names(cv)[cv >= min_cover], names(cs)[cs >= min_cover])
  maps <- lapply(1:2, function(j) {
    d <- ok[which.min(abs(as.numeric(as.Date(ok)) - camp[s, j]))]
    list(date = as.Date(d), r = c(to_lai(vi[[match(d, names(cv))]], vi_mod, paste0(vi_choice$stem, "_closest")),
                                  to_lai(sl[[match(d, names(cs))]], leaf_cal, "leaf_10m_closest")))
  })
  # outer ring of the vegetation polygons, closing the slivers between them
  u <- st_cast(st_union(st_buffer(st_as_sf(v), 1)), "POLYGON")
  u <- st_union(st_sfc(lapply(u, function(p) st_polygon(p[1])), crs = st_crs(u)))
  list(maps = maps, outline = u, tiles = esri(vi[[1]]))
})
lai_max <- quantile(unlist(lapply(lai_maps, function(x) lapply(x$maps, function(y) values(y$r)))), 0.99, na.rm = TRUE)

# graticule ticks on every panel, at least three per axis; latitude labelled on the outline panel only
map_panel <- function(s, layers, title, lab = NULL, labelled = FALSE) {
  x <- lai_maps[[s]]
  e <- as.vector(ext(x$maps[[1]]$r))
  ll <- crds(project(vect(cbind(e[c(1, 2, 1)], e[c(3, 3, 4)]), crs = crs(x$maps[[1]]$r)), "EPSG:4326"))
  grat <- function(r) {
    for (n in 3:8) { b <- scales::breaks_extended(n)(r); b <- b[b > min(r) & b < max(r)]; if (length(b) >= 3) break }
    b
  }
  ggplot() + layers +
    coord_sf(crs = st_crs(x$outline), xlim = e[1:2], ylim = e[3:4], expand = FALSE, label_graticule = "SW") +
    scale_x_continuous(breaks = grat(ll[1:2, 1])) + scale_y_continuous(breaks = grat(ll[c(1, 3), 2])) +
    labs(x = NULL, y = lab, title = title) +
    theme_bw() +
    theme(plot.title = element_text(size = lai_fs, face = "plain"), axis.title.y = element_text(size = lai_fs),
          axis.text = element_text(size = lai_fs - 2.4, colour = "black"), axis.text.x = element_text(angle = 30, hjust = 1),
          axis.text.y = if (!labelled) element_blank(), panel.grid = element_blank())
}

# tan to dark green on a log scale; LAI under 0.05 (incl. zero) takes the lowest colour
lai_layers <- function(r) list(
  geom_spatraster(data = r, maxcell = 1e6),
  scale_fill_gradientn(name = expression(LAI~(m^2~m^-2)), colours = c("#f3e9d2", "#dccb8f", "#a9b85a", "#5e9a3e", "#22693a", "#0b3d2c"),
                       transform = "log10",
                       limits = c(0.05, lai_max), breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2),
                       labels = scales::label_number(drop0trailing = TRUE), oob = scales::squish, na.value = NA,
                       guide = guide_colourbar(theme = theme(legend.key.width = unit(8, "cm")))))

lai_lab <- rep(c("NDVI (NBART)", "SL2P calibrated"), each = 2)
jk <- list(c(1, 1), c(2, 1), c(1, 2), c(2, 2))   # campaign and product of each LAI panel
panels <- unlist(lapply(seq_along(test_sites), function(i) {
  s <- test_sites[i]
  x <- lai_maps[[s]]
  top <- map_panel(s, list(if (!is.null(x$tiles)) geom_spatraster_rgb(data = x$tiles, maxcell = 4e6),
                           geom_sf(data = x$outline, fill = NA, colour = "white", linewidth = 0.5)),
                   if (i == 1) "Study area", s, labelled = TRUE)
  c(list(top), lapply(1:4, function(k) {
    y <- x$maps[[jk[[k]][1]]]
    d <- format(y$date, "%d %b %Y")
    map_panel(s, lai_layers(y$r[[jk[[k]][2]]]), if (i == 1) paste0(lai_lab[k], "\n", d) else d)
  }))
}), recursive = FALSE)

asp <- vapply(lai_maps, function(x) { e <- as.vector(ext(x$maps[[1]]$r)); diff(e[1:2]) / diff(e[3:4]) }, numeric(1))
lai_fig <- wrap_plots(panels, nrow = 4, byrow = TRUE, heights = 1 / asp, guides = "collect") &
  theme(legend.position = "bottom", legend.title.position = "top", legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = lai_fs), legend.title = element_text(size = lai_fs))
ggsave(file.path(plot_path, "test_site_LAI_maps.tif"), lai_fig, width = 30, height = 30.7, units = "cm", dpi = 600,
       compression = "lzw")   # height fills the 30 cm width with panels
message("LAI map images: ", paste(sprintf("%s %s", rep(test_sites, each = 2),
                                          format(do.call(c, lapply(lai_maps, function(x) c(x$maps[[1]]$date, x$maps[[2]]$date))))),
                                  collapse = ", "))

message("imagery credit for the captions: ", get_credit("Esri.WorldImagery"))
