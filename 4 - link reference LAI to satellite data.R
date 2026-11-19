################################################################################
#### Link Sat LAI to veg maps. Get PAI and then convert to LAI from lin reg ####
################################################################################
library(ggplot2)
library(patchwork)
library(terra)
library(future.apply)
library(parallelly)

end_path <- "outputs"

R_LAI <- readRDS(paste0(end_path,"/reference_LAI.rds"))
FV2200_LAI <- read.delim(paste0(end_path,"/FV2200_results_5_rings_27-08-2024.txt"), sep="\t",fill=TRUE,header=T)

# Filtered Satelite images
S2_path = "D:/PilbaraLaiProject/FliteredTimeseries"
S2_site_list <- list.dirs(S2_path,recursive = F,full.names = T)
S2_imgs <- lapply(S2_site_list,function(x){Sys.glob(file.path(x, "*.tif"))})
names(S2_imgs) <- basename(S2_site_list)

# Raw satelite images.
# Sen2Cor L2A reflectance from Earth Engine (script 2g), screened there as the LEAF
# toolbox screens its inputs (QA60, SCL classes 4 and 5; masked pixels are NA). Nothing
# else is applied - no spike, wet-soil or coverage gate.
raw_sat_root <- "D:/PilbaraLaiProject/Sen2Cor_LEAFfilter_stacks"
raw_stack_path <- "D:/PilbaraLaiProject/RawTimeseries_Sen2Cor"
raw_bands <- c("NDVI", "EVI", "MSAVI2", "NIRv", "SR", "DVI", "EVI2", "MSR", "kNDVI", "GCVI", "NDGI")
raw_prefix <- "raw"

# One observation per date in the raw VI arm. Sentinel-2 delivers each overpass
# once per MGRS tile, so a site near a boundary carries the same pass twice.
# Most pixels covered wins, lowest view angle breaks ties.
one_pass_a_day <- TRUE


# Modelled LAI
# rtm_path <- "D:/PilbaraLaiProject/FilteredProsail"        # PROSAIL inversions, filtered by script 3f
leaf_ts_path <- "D:/PilbaraLaiProject/FilteredLeaf"       # SL2P retrievals, filtered by script 3c
leaf_band <- "estimateLAI"

# Vegetation classes over all sites
MVG_shp <- "data/Vegetation NVIS classification 15-3-23.shp"
understorey_PAI_to_LAI <- read.csv("data/PAI to LAI understorey ratios visual assement data.csv",skip=2)
trap_env_data <- read.csv("C:/Users/eck025/OneDrive - UWA/My Documents/PhD 2020/Field sites and data/Litterfall trap env data as at 18-03-2024.csv")
site_env_data <- read.csv("C:/Users/eck025/OneDrive - UWA/My Documents/PhD 2020/Field sites and data/Litterfall analysis/Site data.csv") #### read in enviro data

# Sampling window for summary statistics
avg_window <- 15

#############################################################################
# undo FV2200 naming abbreviation
#############################################################################

sample_period_list <- unique(R_LAI$Sample_period)
site_list <- unique(R_LAI$Site)
MVG_list <- unique(R_LAI$MVG)
# add MVG
FV2200_LAI$MVG <- NA
for (i in 1:length(MVG_list)){
  idx <- which(grepl(paste0(gsub(" ","_",MVG_list[i],fixed=T),"_"),FV2200_LAI$LAI_File))
  FV2200_LAI$MVG[idx] <- MVG_list[i]
}
# add Site
FV2200_LAI$Site <- NA
for (i in 1:length(site_list)){
  idx <- which(grepl(paste0(gsub(" ","_",site_list[i],fixed=T),"_"),FV2200_LAI$LAI_File))
  FV2200_LAI$Site[idx] <- site_list[i]
}
# add sample period
FV2200_LAI$Sample_period <- NA
for (i in 1:length(sample_period_list)){
  idx <- which(grepl(paste0(gsub(" ","_",sample_period_list[i],fixed=T),"_"),FV2200_LAI$LAI_File))
  FV2200_LAI$Sample_period[idx] <- sample_period_list[i]
}

# fix the errors from where MVG names were too big
FV2200_LAI$Strata <- substr(FV2200_LAI$LAI_File,(nchar(FV2200_LAI$LAI_File)+1)-1,nchar(FV2200_LAI$LAI_File))
for (i in 1:length(FV2200_LAI$Strata)){
  if(grepl('Naturally_bare',FV2200_LAI$LAI_File[i])){
    FV2200_LAI$MVG[i] <- MVG_list[4]
    if(FV2200_LAI$LAI_File[i]==FV2200_LAI$LAI_File[i-1]&FV2200_LAI$Strata[i-1]=='G'){
      FV2200_LAI$Strata[i] <- 'C'
    }else{
      FV2200_LAI$Strata[i] <- 'G'
    }
  }
  if(grepl('Inland_Aquatic',FV2200_LAI$LAI_File[i])){
    FV2200_LAI$MVG[i] <- MVG_list[10]
    if(FV2200_LAI$LAI_File[i]==FV2200_LAI$LAI_File[i-1]&FV2200_LAI$Strata[i-1]=='G'){
      FV2200_LAI$Strata[i] <- 'C'
    }else{
      FV2200_LAI$Strata[i] <- 'G'
    }
  }
}
# get n records
n_records <- strsplit(FV2200_LAI$Records," ",fixed=T)
for (i in 1:length(n_records))n_records[[i]] <- as.numeric(gsub("B","",n_records[[i]][2]))
FV2200_LAI$n_records <- do.call(rbind,n_records)
#create dates
FV2200_LAI$Date <- as.Date(as.POSIXct(FV2200_LAI$Date, format="%Y%m%d %H:%M:%S"))

#rename LAI column
colnames(FV2200_LAI)[which(colnames(FV2200_LAI)%in%c('LAI','Le'))] <- c('PAI_LX_LAI2200','PAIe_LAI2200')

#############################################################################
# get LAI from overstory PAI using Hemi-cam model
#############################################################################
# note - this results in values below 0 in some cases!!
# read in model and LAI datafiles from R
LAI_all <- readRDS("C:/Users/eck025/OneDrive - UWA/My Documents/PhD 2020/LAI_validation_analysis/Hemi_cam_PAI_and_LAI_data.rds")
PAI2200toLAI_lm <- readRDS("C:/Users/eck025/OneDrive - UWA/My Documents/PhD 2020/LAI_validation_analysis/PAI2200toLAI_lm.rds")
PAI2200toLAI_lm_byGen <- readRDS("C:/Users/eck025/OneDrive - UWA/My Documents/PhD 2020/LAI_validation_analysis/PAI2200toLAI_lm.rds")



#alternative groupings
# old.names <- as.character(unique(canopy_LAI$MVG))
# new.names <- c('Acacia dominated','Eucalupt dominated','Mixed grasslands','Mixed shrublands','Melaleuca dominated','Acacia dominated','Mixed grasslands','Acacia dominated','Eucalupt dominated')
# dom.struct <- c('Forest and woodlands','Forest and woodlands','Grasslands','Shrublands','Forest and woodlands','Shrublands','Grasslands','Forest and woodlands','Forest and woodlands')
# 
# canopy_LAI$MVG.combined <- NA
# canopy_LAI$veg.strata <- NA
# for(i in 1:nrow(canopy_LAI)){
#   idx <- which(old.names==canopy_LAI$MVG[i])
#   canopy_LAI$MVG.combined[i] <- new.names[idx]
#   canopy_LAI$veg.strata[i] <- dom.struct[idx]
# }


#split out overstorey data
canopy_LAI <- FV2200_LAI[FV2200_LAI$Strata=="C",]

#apply model to estimate LAI
canopy_LAI$LAI_LX_canopy <- predict(PAI2200toLAI_lm,newdata = canopy_LAI)
canopy_LAI$LAI_LX_canopy[canopy_LAI$LAI_LX_canopy<0] <- 0
colnames(canopy_LAI)[which(colnames(canopy_LAI)=='n_records')] <- 'n_records_canopy'
colnames(canopy_LAI)[which(colnames(canopy_LAI)=='PAI_LX_LAI2200')] <- 'PAI_LX_LAI2200_canopy'
colnames(canopy_LAI)[which(colnames(canopy_LAI)=='PAIe_LAI2200')] <- 'PAIe_LAI2200_canopy'

#############################################################################
# calculate understorey PAI and convert to LAI using %greenness dataset
#############################################################################
understorey_PAI_to_LAI$Date <- as.Date(understorey_PAI_to_LAI$Date,format = "%d/%m/%Y")
understorey_PAI_to_LAI$SamplePeriod <- ifelse(understorey_PAI_to_LAI$Date<'2022-01-01'&understorey_PAI_to_LAI$Date>'2021-09-01','October21',NA)
understorey_PAI_to_LAI$SamplePeriod <- ifelse(understorey_PAI_to_LAI$Date>'2022-01-01'&understorey_PAI_to_LAI$Date<'2022-04-01','March22',understorey_PAI_to_LAI$SamplePeriod)

index2 <- which(colnames(trap_env_data)%in%c('Trap.code','MVG.Name'))
trap_env_data <- trap_env_data[,index2]
Site.cluster <- site_env_data[,which(colnames(site_env_data)%in%c("Site","Site.cluster"))]

understorey_PAI_to_LAI <- merge(understorey_PAI_to_LAI,trap_env_data,by.x='Code',by.y='Trap.code')
understorey_PAI_to_LAI <- merge(understorey_PAI_to_LAI,Site.cluster)


# find the unique over/understorey combinations
colnames(understorey_PAI_to_LAI)[5] <- 'prop'
v <- vect(MVG_shp)
plot(v[v$Site=='CS',])
combs <- strsplit(unique(paste0(v$Site,'_',v$MVG_Nam,'_',v$Undrstr)),'_')
combs <- do.call(rbind,
  lapply(combs,function(x){
    df <- as.data.frame(t(x))
    colnames(df) <- c('Site','MVG.Name','Understorey')
    return(df)
    })
)

combs <- combs[combs$Site%in%FV2200_LAI$Site,]
combs <- combs[order(combs$Site),]
combs$Understorey[combs$Understorey=='tussock grassland'] <- 'tussock'
combs$Understorey[combs$Understorey=='hummock grass'] <- 'hummock'
combs$Understorey[combs$Understorey=='sedgelands'] <- 'sedge'

# mapped area of each understorey class, to weight the conversions where one MVG is mapped with two classes
und_area <- data.frame(Site = v$Site, MVG.Name = v$MVG_Nam, Understorey = v$Undrstr, area_ha = expanse(v, unit = "ha"))
und_area$Understorey[und_area$Understorey=='tussock grassland'] <- 'tussock'
und_area$Understorey[und_area$Understorey=='hummock grass'] <- 'hummock'
und_area$Understorey[und_area$Understorey=='sedgelands'] <- 'sedge'
und_area <- aggregate(area_ha ~ Site + MVG.Name + Understorey, data = und_area, FUN = sum)
combs <- merge(combs, und_area, by = c('Site','MVG.Name','Understorey'), all.x = TRUE)

# add the sampling period
SampPrd <- strsplit(unique(paste0(FV2200_LAI$Date,'_',FV2200_LAI$Site)),'_')
SampPrd <- do.call(rbind,
                   lapply(SampPrd,function(x){
                     df <- as.data.frame(t(x))
                     colnames(df) <- c('Date','Site')
                     return(df)
                     })
                   )
SampPrd$SamplePeriod <- ifelse(SampPrd$Date<'2022-01-01'&SampPrd$Date>'2021-09-01','October21',NA)
SampPrd$SamplePeriod <- ifelse(SampPrd$Date>'2022-01-01'&SampPrd$Date<'2022-04-01','March22',SampPrd$SamplePeriod)
combs <- merge(combs,SampPrd,by='Site')

# backup - use nearby sites
combs <- merge(combs,Site.cluster)

# empty cols for data
combs$conversion <- NA
combs$sd <- NA

# find the best conversion
for(i in 1:nrow(combs)){
  site <- combs$Site[i]
  mvg <- combs$MVG.Name[i]
  us <- combs$Understorey[i]
  sp <- combs$SamplePeriod[i]
  clust <- combs$Site.cluster[i]
  
  # 1. optimum - same site, mvg, & understorey
  idx <- which(understorey_PAI_to_LAI$Site == site &
                 understorey_PAI_to_LAI$Understorey == us &
                 understorey_PAI_to_LAI$SamplePeriod == sp &
                 understorey_PAI_to_LAI$MVG.Name == mvg)
  
  if(length(idx)>1){
    combs$conversion[i] <- mean(understorey_PAI_to_LAI$prop[idx],na.omit=T)
    combs$sd[i] <- sd(understorey_PAI_to_LAI$prop[idx])
    next
  }
  
  # 2. same site & understorey
  idx <- which(understorey_PAI_to_LAI$Site == site &
                 understorey_PAI_to_LAI$Understorey == us &
                 understorey_PAI_to_LAI$SamplePeriod == sp)
  
  if(length(idx)>1){
    combs$conversion[i] <- mean(understorey_PAI_to_LAI$prop[idx],na.omit=T)
    combs$sd[i] <- sd(understorey_PAI_to_LAI$prop[idx])
    next
  }
  
  # 3. same site
  idx <- which(understorey_PAI_to_LAI$Site == site &
                 understorey_PAI_to_LAI$SamplePeriod == sp)
  
  if(length(idx)>1){
    combs$conversion[i] <- mean(understorey_PAI_to_LAI$prop[idx],na.omit=T)
    combs$sd[i] <- sd(understorey_PAI_to_LAI$prop[idx])
    next
  }
  
  # 4. sites within the cluster (mvg + us)
  idx <- which(understorey_PAI_to_LAI$Site.cluster == clust &
                 understorey_PAI_to_LAI$Understorey == us &
                 understorey_PAI_to_LAI$SamplePeriod == sp &
                 understorey_PAI_to_LAI$MVG.Name == mvg)
  
  if(length(idx)>1){
    combs$conversion[i] <- mean(understorey_PAI_to_LAI$prop[idx],na.omit=T)
    combs$sd[i] <- sd(understorey_PAI_to_LAI$prop[idx])
    next
  }
  
  # 5. sites within the cluster (us)
  idx <- which(understorey_PAI_to_LAI$Site.cluster == clust &
                 understorey_PAI_to_LAI$Understorey == us &
                 understorey_PAI_to_LAI$SamplePeriod == sp)
  
  if(length(idx)>1){
    combs$conversion[i] <- mean(understorey_PAI_to_LAI$prop[idx],na.omit=T)
    combs$sd[i] <- sd(understorey_PAI_to_LAI$prop[idx])
    next
  }
  
  # 6. sites outside the cluster (us)
  idx <- which(understorey_PAI_to_LAI$Understorey == us &
                 understorey_PAI_to_LAI$SamplePeriod == sp)
  
  if(length(idx)>1){
    combs$conversion[i] <- mean(understorey_PAI_to_LAI$prop[idx],na.omit=T)
    combs$sd[i] <- sd(understorey_PAI_to_LAI$prop[idx])
    next
  }
}



# apply the conversion
understorey_LAI <- FV2200_LAI[!FV2200_LAI$Strata=="C",]

# get understorey PAI values
for (i in 1:nrow(understorey_LAI)){
  idx <- which(canopy_LAI$Date==understorey_LAI$Date[i]&
                 canopy_LAI$Site==understorey_LAI$Site[i]&
                 canopy_LAI$MVG==understorey_LAI$MVG[i])
  understorey_LAI$PAI_understorey[i] <- understorey_LAI$PAI_LX_LAI2200[i]-canopy_LAI$PAI_LX_LAI2200[idx]
  if(understorey_LAI$PAI_understorey[i]<0){
    understorey_LAI$PAI_understorey[i] <- 0
  }
  if(is.na(understorey_LAI$PAI_understorey[i])){
    print(paste0('missing canopy record from i = ',i))
  }
}

combs$conversion[is.na(combs$conversion)] <- mean(combs$conversion,na.rm=T)

colnames(combs)[colnames(combs)=="MVG.Name"] <- "MVG"

understorey_LAI <- merge(understorey_LAI[colnames(understorey_LAI)!="Sample_period"],combs,by=c('Date','Site','MVG'))

understorey_LAI$LAI_LX_understorey <- understorey_LAI$PAI_understorey * understorey_LAI$conversion

hist(understorey_LAI$LAI_LX_understorey)
understorey_LAI[understorey_LAI$LAI_LX_understorey>1,]
colnames(understorey_LAI)[colnames(understorey_LAI)%in%c('n_records',"conversion","sd")] <- c('n_records_understorey',"conversion_us","sd_us")

# One reference value per Site x MVG x Date. Where an MVG is mapped with two understorey classes (hummock and
# tussock at CS, shrub and tussock at HD) the one LAI-2200C measurement has a conversion for each class and the
# merge above made a row for each. They are not two visits, so the understorey LAI is the mean of the
# conversions weighted by the mapped area of each class.
understorey_LAI$area_ha[is.na(understorey_LAI$area_ha)] <- 1       # single-class MVGs: the weight cancels
understorey_LAI <- do.call(rbind, lapply(split(understorey_LAI, paste(understorey_LAI$Site, understorey_LAI$MVG, understorey_LAI$Date)), function(x) {
  y <- x[1,]
  y$LAI_LX_understorey <- weighted.mean(x$LAI_LX_understorey, x$area_ha, na.rm = TRUE)
  y$conversion_us      <- weighted.mean(x$conversion_us, x$area_ha, na.rm = TRUE)
  y$sd_us              <- weighted.mean(x$sd_us, x$area_ha, na.rm = TRUE)
  y$Understorey        <- paste(sort(unique(x$Understorey)), collapse = '+')
  y
}))
rownames(understorey_LAI) <- NULL



#############################################################################
# combine to get total values
#############################################################################
total_LAI <- list()
for (i in 1:nrow(understorey_LAI)){
  idx <- which(canopy_LAI$Date==understorey_LAI$Date[i]&
               canopy_LAI$Site==understorey_LAI$Site[i]&
               canopy_LAI$MVG==understorey_LAI$MVG[i])
  
  # understorey_LAI holds one row per stratum and date (the area-weighted conversions above);
  # Understorey names the class or classes that were averaged
  total_LAI[[i]] <- as.data.frame(cbind(understorey_LAI[i,which(colnames(understorey_LAI)%in%c('Site','Date','MVG','Understorey','n_records_understorey','PAI_understorey','LAI_LX_understorey'))],
                   canopy_LAI[idx,which(colnames(canopy_LAI)%in%c('n_records_canopy','PAI_understorey','LAI_LX_canopy'))]))
}
total_LAI <- do.call(rbind,total_LAI)
total_LAI$LAI_LX_total <- total_LAI$LAI_LX_canopy+total_LAI$LAI_LX_understorey


#############################################################################
#############################################################################
# combine with satelite data
#############################################################################
#############################################################################
site_date_str <- unique(paste0(total_LAI$Site,"_",total_LAI$Date))
site_date_list <- list()
for(i in 1:length(site_date_str)){
  site_date_list[[i]] <- as.data.frame(t(strsplit(site_date_str[[i]],'_',fixed=T)[[1]]))
  colnames(site_date_list[[i]]) <- c('Site','Date')
}

get_sum_stats <- function(x,# = site_date_list[[1]],
                          sat_img_list = S2_imgs,
                          Veg_class,
                          days_buffer=avg_window,
                          na_thresh = 0.2, # fraction of NA or non-valid pixels across entire image that rejects closest in time
                          return_rast = F){
  
  #Set project CRS
  proj_crs<-"EPSG:4326"
  
  samp_date <- x[,which(colnames(x)=='Date')]
  site <- x[,which(colnames(x)=='Site')]
  
  if(is.character(samp_date)==T) samp_date <- as.Date(samp_date)
  
  # get the rast images matching the site code
  site_sat_img_list <- sat_img_list[which(site==names(sat_img_list))][[1]]
  
  # limit shape to the right site
  Veg_class <- vect(Veg_class)
  Site_MVG_shape <- Veg_class[which(Veg_class[['Site']]==site),]
  
  # loop through each layer
  sum_rast <- list()
  for(layer_idx in 1:length(site_sat_img_list)){ # loop through each individual band
    # read in as rast and fix projection
    S2_MS_img <- rast(site_sat_img_list[layer_idx])
    S2_MS_img <- terra::project(S2_MS_img,proj_crs)
    
    band <- unlist(strsplit(basename(site_sat_img_list[layer_idx]),'-'))[1]
    
    # get the dates of each capture
    date_list <- as.Date(substr(names(S2_MS_img), 1, 10))
    
    # date range
    min_date <- samp_date - days_buffer
    max_date <- samp_date + days_buffer
    date_range_idx <- which(date_list > min_date & date_list < max_date)
    closest_date_idx <- order(as.numeric(date_list - samp_date) ^2)
    
    S2_DateRange <- S2_MS_img[[date_range_idx]]
    
    #check how much of the image is NA and remove them if >20% is NA
    retain_idx <- list()
    for(xx in 1:length(date_range_idx)){
      n_NA <- global(S2_DateRange[[xx]], fun="isNA")
      n_pixels <- dim(S2_DateRange[[xx]])[1]*dim(S2_DateRange[[xx]])[2]
      if(n_NA<n_pixels * na_thresh){
        retain_idx[xx] <- xx
      }else{
        retain_idx[xx] <- NULL
      }
    }
    retain_idx <- unlist(retain_idx)
    
    if (length(retain_idx)) S2_DateRange <- S2_DateRange[[retain_idx]]
    
    # read in the closest in time image and remove it if there's too many NAs
    # walk the candidates in order of closeness, take the first usable one, making sure it's within 7 days of sampling
    im <- 1L
    ok <- FALSE
    while (im <= length(closest_date_idx)) {
      cand <- S2_MS_img[[closest_date_idx[im]]]
      n_NA <- global(cand, fun = "isNA")
      n_pixels <- dim(cand)[1] * dim(cand)[2]
      if (n_NA <= n_pixels * na_thresh &&
          abs(date_list[closest_date_idx[im]] - samp_date) <= 7) {
        S2_MS_img_closest <- cand
        ok <- TRUE
        break
      }
      im <- im + 1L
    }
    if (!ok) {
      S2_MS_img_closest <- S2_MS_img[[closest_date_idx[1]]]
      values(S2_MS_img_closest) <- NA
      message('error: ', site, ' failed NA flag')
    }
    
    sum_rast[[layer_idx]] <- rast(list(mean(S2_DateRange,na.rm=T),
                                       median(S2_DateRange,na.rm=T),
                                       stdev(S2_DateRange,na.rm=T),
                                       max(S2_DateRange,na.rm=T),
                                       min(S2_DateRange,na.rm=T),
                                       S2_MS_img_closest))
    
    labels <- c('mean','median','temp_sd','max','min','closest')
    labels <- paste0(band,"_",labels)
    names(sum_rast[[layer_idx]]) <- labels
  }
  
  # Use one of the Sentinel-2 images as a reference grid
  ref_grid <- sum_rast[[1]][[1]]

  # downsamples the 20 m to identical 10 m cells for common reference.
  for (k in seq_along(sum_rast)) {
    if (!terra::compareGeom(sum_rast[[k]], ref_grid, stopOnError = FALSE,
                            messages = FALSE)) {
      keep_names <- names(sum_rast[[k]])
      sum_rast[[k]] <- terra::resample(sum_rast[[k]], ref_grid, method = "near")
      names(sum_rast[[k]]) <- keep_names
    }
  }
  
  # add the # of passes for each sample
  n_passes <- length(retain_idx)
  
  # Rasterize the shapefile to match the S2 grid, assign the maximum overlapping vegetation class ID to each pixel
  # could also calculate a weighting of each - might be better but this seems defensible
  sum_rast[[layer_idx+1]] <- terra::rasterize(Site_MVG_shape, ref_grid, field = "MVG_Nam", fun = "max", background = NA)
  sum_rast <- rast(sum_rast)
  sum_df <- as.data.frame(sum_rast,xy=T)
  sum_df <- sum_df[!is.na(sum_df[['MVG_Nam']]),]
  
  # take the mean of each veg class
  # how many pixels each class holds when nothing is masked, for the gate below
  class_total <- table(as.data.frame(sum_rast[['MVG_Nam']], na.rm = TRUE)[[1]])

  MVG_list <- unique(Veg_class$MVG_Nam)
  sum_stat <- list()
  for(xx in 1:length(MVG_list)){
    index <- which(sum_df$MVG_Nam==MVG_list[xx])
    if(length(index)>0){
      sum_df_MVG <- sum_df[index,-c(which(colnames(sum_df)%in%c('x','y','MVG_Nam')))]
      sum_stat[[xx]] <- cbind(as.data.frame(t(colMeans(sum_df_MVG))))
      sd <- cbind(as.data.frame(t(apply(sum_df_MVG[,which(grepl('mean',colnames(sum_df_MVG)))],2,sd)))) # basically the same concept as colmeans but applies sd to the 2nd dimension (columns)
      colnames(sd) <- gsub('_mean','_sd_spat',colnames(sd))
      n_pixels <- nrow(sum_df_MVG)
      
      # Fraction of the pixels from the veg class subset actually retained on this date. as.data.frame
      # Dates surviving on a small, spatially biased remnant
      # of the polygon was the main source of spikes in the model timeseries.
      cls_tot <- as.numeric(class_total[MVG_list[xx]])
      frac_valid <- if (is.na(cls_tot) || cls_tot == 0) NA_real_ else n_pixels / cls_tot
      sd$class_pixels <- cls_tot
      sd$frac_valid   <- round(frac_valid, 4)
      sd <- cbind(sd,n_pixels)
      sd$n_passes <- n_passes
      sum_stat[[xx]] <- cbind(sum_stat[[xx]],sd)
      sum_stat[[xx]]$MVG_Name <- MVG_list[xx]
      
      # above the gate the row is kept but emptied, so the merge still records
      # that the date was visited and why it carries no value
      # na_thresh is the tolerated fraction of invalid pixels, as at the gates above
      # so the row is emptied once more than that share was masked.
      if (!is.na(frac_valid) && (1 - frac_valid) > na_thresh) {
        v <- grep('_mean$|_median$|_temp_sd$|_max$|_min$|_closest$|_sd_spat$',
                  colnames(sum_stat[[xx]]), value = TRUE)
        sum_stat[[xx]][, v] <- NA_real_
      }
    }else{
      sum_stat[[xx]] <- NULL
    }
  }
  sum_stat <- do.call(rbind,sum_stat)
  
  sum_stat$Site <- site
  sum_stat$Date <- samp_date
  sum_df$Site <- site
  sum_df$Date <- samp_date
  
  output <- list(sum_stat,sum_df)
  names(output) <- c('mean_summary_statistics','raw_pixel_statistics')
  if(return_rast==T){
    return(sum_rast)
  }else{
    return(output)
  }
}


##############################################################################################################
#open multisession
n_workers <- min(8L, parallelly::availableCores(omit = 2),length(site_list))
plan(multisession, workers = n_workers)

  output_list <- future_lapply(site_date_list,
                               get_sum_stats,
                               Veg_class = MVG_shp,
                               sat_img_list = S2_imgs,
                               days_buffer = avg_window,
                               future.seed = TRUE,
                               future.packages = c("terra"))

#close multisession
plan(sequential)
##############################################################################################################
# Save file
saveRDS(output_list,file = paste0(end_path,"/VI_summary_stats_for_validation.rds"))

output_list <- readRDS(paste0(end_path,"/VI_summary_stats_for_validation.rds"))
# Extract each list of DFs
Sat_timeser <- do.call(rbind,lapply(output_list,FUN=function(x)x=x[[1]]))# extract the summary stats for each veg type
#temp_sum_stats <- lapply(output_list,FUN=function(x)x=x[[2]]) # this is just the summary stats
colnames(Sat_timeser)[which(colnames(Sat_timeser)=='MVG_Name')] <- 'MVG'
total_LAI <- merge(total_LAI,Sat_timeser,by=c('Site','MVG','Date'),all.x=T)

saveRDS(total_LAI,paste0(end_path,"/paired_ref_and_sat_dataset.rds"))




################################################################################
#### Model-predicted LAI: LEAF-Toolbox (script 2f) and PROSAIL (script 3e)   ####
################################################################################
# SL2P products come from script 3c and PROSAIL from 3f, both already QC filtered

##############################################################################
# pull the band out of the raw VIs to match filtered timeseries.
##############################################################################
source("R/stack_io.R")
sat_root  <- raw_sat_root
stack_idx <- stack_index(sat_root)

# One layer per date, chosen on the pixels that survive the masks (NA in the
# stacks). NULL per site means take every layer.
raw_sel <- NULL
if (one_pass_a_day && !is.null(stack_idx)) {
  sites_i <- sort(unique(stack_idx$Site))
  raw_sel <- setNames(lapply(sites_i, function(s) {
    b4 <- try(site_band_series(s, "B4", stack_idx), silent = TRUE)
    m  <- site_meta(s, stack_idx)
    if (inherits(b4, "try-error") || is.null(b4) || is.null(m)) return(NULL)
    if (nlyr(b4) != nrow(m)) return(NULL)
    cov <- 1 - colMeans(is.na(values(b4)))
    sort(unname(vapply(split(seq_along(cov), m$Date), function(j)
      j[order(-cov[j], m$MEAN_INCIDENCE_ZENITH_ANGLE_B8A[j])][1], integer(1))))
  }), sites_i)
}

if (!is.null(stack_idx)) {
  for (s in sort(unique(stack_idx$Site))) {
    outs <- file.path(raw_stack_path, s, paste0(raw_prefix, raw_bands, "-ts.tif"))
    if (all(file.exists(outs))) next
    dir.create(file.path(raw_stack_path, s), showWarnings = FALSE, recursive = TRUE)
    n_written <- 0L
    for (k in seq_along(raw_bands)) {
      if (file.exists(outs[k])) next
      r <- try(site_vi_series(s, raw_bands[k], stack_idx, lyrs = raw_sel[[s]]),
               silent = TRUE)
      if (inherits(r, "try-error") || is.null(r)) {
        message("  no ", raw_bands[k], " for ", s); next
      }
      writeRaster(r, outs[k], overwrite = TRUE)
      n_written <- n_written + 1L
    }
    if (n_written) message("re-stacked raw ", s, ": ", n_written, " bands")
  }
}
raw_site_dirs <- list.dirs(raw_stack_path, recursive = FALSE, full.names = TRUE)
raw_S2_imgs <- lapply(raw_site_dirs, function(d) Sys.glob(file.path(d, "*-ts.tif")))
names(raw_S2_imgs) <- basename(raw_site_dirs)
raw_S2_imgs <- raw_S2_imgs[lengths(raw_S2_imgs) > 0]


##############################################################################
# compile raw (unfiltered) VIs
##############################################################################
raw_site_date_list <- site_date_list[
  vapply(site_date_list, function(x) x[, "Site"] %in% names(raw_S2_imgs), logical(1))]

if (length(raw_site_date_list)) {
  plan(multisession, workers = n_workers)
  raw_output_list <- future_lapply(raw_site_date_list,
                                   get_sum_stats,
                                   Veg_class = MVG_shp,
                                   sat_img_list = raw_S2_imgs,
                                   days_buffer = avg_window,
                                   future.seed = TRUE,
                                   future.packages = c("terra"))
  plan(sequential)
  saveRDS(raw_output_list, file = paste0(end_path, "/raw_VI_summary_stats_for_validation.rds"))
} else {
  raw_output_list <- NULL
  message("no raw stacks found - skipping the raw VI summaries")
}

# ##############################################################################
# # summarise custom RadTransModel products
# ##############################################################################
# rtm_sites <- sort(unique(basename(list.dirs(rtm_path,     recursive = FALSE, full.names = TRUE))))
# RTM_imgs <- lapply(rtm_sites, function(s) {
#   f <- Sys.glob(file.path(rtm_path,     s, "*-ts.tif"))
#   f[!endsWith(f, "_qc-ts.tif")]        # QC stacks are read alongside, not here
# })
# names(RTM_imgs) <- rtm_sites
# RTM_imgs <- RTM_imgs[lengths(RTM_imgs) > 0]
# message("model products per site: ",
#         paste(unique(unlist(lapply(RTM_imgs, function(f)
#           sub("-ts.tif", "", basename(f), fixed = TRUE)))), collapse = ", "))
# 
# # only the sites that actually have reference data
# rtm_site_date_list <- site_date_list[
#   vapply(site_date_list, function(x) x[, "Site"] %in% names(RTM_imgs), logical(1))]
# 
# n_workers <- min(8L, parallelly::availableCores(omit = 2), length(rtm_site_date_list))
# plan(multisession, workers = n_workers)
# rtm_output_list <- future_lapply(rtm_site_date_list,
#                                  get_sum_stats,
#                                  Veg_class = MVG_shp,
#                                  sat_img_list = RTM_imgs,
#                                  days_buffer = avg_window,
#                                  future.packages = c("terra"))
# plan(sequential)
# 
# saveRDS(rtm_output_list, file = paste0(end_path, "/RTM_summary_stats_for_validation.rds"))
# 
##############################################################################
# summarise LEAF products
##############################################################################
leaf_sites <- sort(unique(basename(list.dirs(leaf_ts_path, recursive = FALSE, full.names = TRUE))))
leaf_imgs <- lapply(leaf_sites, function(s) {
  f <- Sys.glob(file.path(leaf_ts_path, s, "*-ts.tif"))
  f[!endsWith(f, "_qc-ts.tif")]        # QC stacks are read alongside, not here
})
names(leaf_imgs) <- leaf_sites
leaf_imgs <- leaf_imgs[lengths(leaf_imgs) > 0]
message("model products per site: ",
        paste(unique(unlist(lapply(leaf_imgs, function(f)
          sub("-ts.tif", "", basename(f), fixed = TRUE)))), collapse = ", "))

# only the sites that actually have reference data
leaf_site_date_list <- site_date_list[
  vapply(site_date_list, function(x) x[, "Site"] %in% names(leaf_imgs), logical(1))]

n_workers <- min(8L, parallelly::availableCores(omit = 2), length(leaf_site_date_list))
plan(multisession, workers = n_workers)
leaf_output_list <- future_lapply(leaf_site_date_list,
                                 get_sum_stats,
                                 Veg_class = MVG_shp,
                                 sat_img_list = leaf_imgs,
                                 days_buffer = avg_window,
                                 future.seed = TRUE,
                                 future.packages = c("terra"))
plan(sequential)

saveRDS(leaf_output_list, file = paste0(end_path, "/LEAF_summary_stats_for_validation.rds"))

##############################################################################
# merge everything: reference LAI + VI summary stats + model-predicted LAI
##############################################################################
# PROSAIL is not carried into the analysis; see scripts 3e and 3f.
# RTM_timeser <- do.call(rbind, lapply(rtm_output_list, function(x) x[[1]]))
# colnames(RTM_timeser)[which(colnames(RTM_timeser) == 'MVG_Name')] <- 'MVG'
# 
# # n_pixels / n_passes already came through with the VI merge
# drop <- intersect(c("n_pixels", "n_passes", "class_pixels", "frac_valid"), colnames(RTM_timeser))
# RTM_timeser <- RTM_timeser[, setdiff(colnames(RTM_timeser), drop)]
# 
# total_LAI <- merge(total_LAI, RTM_timeser, by = c('Site', 'MVG', 'Date'), all.x = TRUE)

# leaf outputs ------------
leaf_timeser <- do.call(rbind, lapply(leaf_output_list, function(x) x[[1]]))
colnames(leaf_timeser)[which(colnames(leaf_timeser) == 'MVG_Name')] <- 'MVG'

drop <- intersect(c("n_pixels", "n_passes", "class_pixels", "frac_valid"), colnames(leaf_timeser))
leaf_timeser <- leaf_timeser[, setdiff(colnames(leaf_timeser), drop)]

total_LAI <- merge(total_LAI, leaf_timeser, by = c('Site', 'MVG', 'Date'), all.x = TRUE)

# raw VIs -----------------
if (!is.null(raw_output_list)) {
  raw_timeser <- do.call(rbind, lapply(raw_output_list, function(x) x[[1]]))
  colnames(raw_timeser)[which(colnames(raw_timeser) == 'MVG_Name')] <- 'MVG'
  drop <- intersect(c("n_pixels", "n_passes", "class_pixels", "frac_valid"), colnames(raw_timeser))
  raw_timeser <- raw_timeser[, setdiff(colnames(raw_timeser), drop)]
  total_LAI <- merge(total_LAI, raw_timeser, by = c('Site', 'MVG', 'Date'), all.x = TRUE)
}

saveRDS(total_LAI, paste0(end_path, "/paired_ref_and_sat_dataset.rds"))
message("merged dataset: ", nrow(total_LAI), " rows, ", ncol(total_LAI), " columns")
message("  model LAI columns: ",
        paste(grep("prosail|leaf", colnames(total_LAI), value = TRUE), collapse = ", "))
message("  raw VI columns   : ",
        length(grep(paste0("^", raw_prefix), colnames(total_LAI))))
