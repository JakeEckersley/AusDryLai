library(terra)
dir.create('data')

# moving files
# NVIS_map <- "C:/Users/eck025/OneDrive - UWA/My Documents/PhD 2020/Field sites and data/GIS/Vegetation NVIS classification 15-3-23.shp"
# vect_data <- vect(NVIS_map)
# writeVector(vect_data,"data/Vegetation NVIS classification 15-3-23.shp",overwrite=TRUE)

NVIS_map <- "data/Vegetation NVIS classification 15-3-23.shp"
vect_data <- vect(NVIS_map)

sites <- unique(vect_data$Site)
bbox_list <- vector("list", length(sites))

dir.create('data/SiteBoundsForSatDownloads',showWarnings = F)
for (i in seq_along(sites)) {
  site_sub <- vect_data[vect_data$Site == sites[i], ]
  site_sub <- buffer(site_sub,150)
  site_ext <- ext(site_sub)
  bbox_poly <- as.polygons(site_ext)
  bbox_poly$Site <- sites[i]
  bbox_list[[i]] <- bbox_poly
  bbox_poly <- set.crs(bbox_poly,vect_data)
  
  writeVector(bbox_poly,filename = paste0('data/SiteBoundsForSatDownloads/',sites[i],'.shp'),overwrite=T)
}

site_bboxes <- do.call(rbind, bbox_list)
site_bboxes <- set.crs(site_bboxes,vect_data)

pilbara_bbox <- buffer(site_bboxes,9000)
pilbara_bbox <- ext(pilbara_bbox)
pilbara_bbox <- as.polygons(pilbara_bbox)
pilbara_bbox <- set.crs(pilbara_bbox,vect_data)
writeVector(pilbara_bbox,filename = paste0('data/PilbaraBbox.shp'),overwrite=T)

##########################
# Kimberley included
##########################
nw_sites <- vect("H:/ecology/08_Phase2_GBINCB/04_GIS&Modelling/FieldworkGISLayersMaster/Sites/QuadratLocations.shp")
nw_bbox <- ext(nw_sites)
nw_bbox <- as.polygons(nw_bbox)
nw_bbox <- set.crs(nw_bbox,vect_data)
# writeVector(nw_bbox,filename = paste0('data/MultiProjBbox.shp'),overwrite=T)

Proj_bboxes <- rbind(pilbara_bbox,nw_bbox)
Proj_bboxes <- as.polygons(ext(Proj_bboxes))
Proj_bboxes <- set.crs(Proj_bboxes,vect_data)
writeVector(Proj_bboxes,filename = paste0('data/MultiProjBbox.shp'),overwrite=T)
