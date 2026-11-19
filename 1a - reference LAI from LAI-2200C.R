############################################################
############################################################
############################################################
## This script is to split transects into vegetation classes before being pushed through FV2200
## First go through and clean text files, labelling all scatter readings with a K
## Update: 27/08/24 added code to spit out the combined a/b text file for each veg class to run in FV2200
## Also calculates total LAI

library(data.table)
library(sf)
#library(terra)
#library(dplyr)

#############################
# Functions
#############################
#### 1. function to implement LAI equations
#### edit this to change the minimum based on a pre-defined value....
LAI_totals<-function(ABpairs,Strata,clip=T,exclud_ring=5,sat_LAI=F){
  
  W<-c(0.041,0.131,0.201,0.29,0.337) # weighting factors Eq. 10.15
  
  # recalculate weighting if ring is excluded
  if(is.numeric(exclud_ring)&!exclud_ring==0){
    W<-W[-exclud_ring]/sum(W[-exclud_ring])
  }
  
  angles<-c(7,23,38,53,68) #
  angles_rad<-angles * ( pi / 180.0 )
  
  b<-ABpairs[,5:9] # extract b readings
  
  
  if(is.numeric(sat_LAI)){
    max_trans<-list()
    for (rngs in 1:length(angles_rad)){
      max_trans[[rngs]] <- exp(-sat_LAI / (2 * cos(angles_rad[rngs])))
    }
    max_trans<-do.call(rbind,max_trans)
  }

  dodgy_reading<-b==0
  if(length(which(dodgy_reading))>0){
    dodgy_reading<-list()
    for(cols in 1:length(W)){
      idx<-which(b[,cols]==0)
      b[which(b[,cols]==0),]<-max_trans[cols]
      if(length(idx)>0){
        dodgy_reading[[cols]]<-data.frame(dataset=unique(ABpairs[,22]),strata=Strata,column=cols,row=as.character(ABpairs[idx,3]))
      }
    }
    dodgy_reading<-do.call(rbind,dodgy_reading)
    print(paste0(nrow(dodgy_reading),' saturated readings in ',unique(ABpairs$Site)," ",unique(ABpairs$MVG.Name)," ",Strata))
    print(dodgy_reading)
  }

  idx<-which(colnames(ABpairs) %in% c("A_Ring 1","A_Ring 2","A_Ring 3","A_Ring 4","A_Ring 5" ))
  a<-ABpairs[,idx]
  a[a==0]<-NA
  

  ############
  #Eq 10.14
  LAI_pairs<-as.data.frame(matrix(ncol=length(W))) #blank end df
  transmittance<-b/a
  
  
  if(clip==T){
    clipped_readings<-table(transmittance[,1:length(W)]>1)['TRUE']/sum(table(transmittance>1)['FALSE']+table(transmittance>1)['TRUE'])
    transmittance[transmittance>1]<-1
  }else{
    transmittance[transmittance>1]<-NA
    removed_readings<-nrow(transmittance)-nrow(na.omit(transmittance))
    
    #remove anything with an NA
    transmittance<-na.omit(transmittance)
  }
  


  # adjust trans as per saturated LAI
  if(is.numeric(sat_LAI)){
    for (i in 1:nrow(transmittance)){
      for (j in 1:ncol(transmittance)){
        if(transmittance[i,j]<max_trans){
          transmittance[i,j]<-max_trans[j]
        }
      }
    }
  }
  
  
  for (ring in 1:length(W)){
    for (i in 1:nrow(transmittance)){
      LAI_pairs[i,ring]<--log(transmittance[i,ring])/(1/cos(angles_rad[ring]))
    }
  }
  ln_K<-colSums(LAI_pairs)*(1/nrow(LAI_pairs))
  rng_tot<-ln_K*W
  
  LAI<-sum(rng_tot)*2
  ############
  
  ########################
  #Eq. 10-18
  ring_KW<-as.data.frame(matrix(ncol = ncol(LAI_pairs),nrow=nrow(LAI_pairs)))
  for(ring in 1:length(W)){
    for (i in 1:nrow(LAI_pairs)){
      ring_KW[i,ring]<-LAI_pairs[i,ring]*W[ring]
    }
  }
  pointLAI<-rowSums(ring_KW)*2
  LAIse<-sqrt(sum((pointLAI-LAI)^2)/length(pointLAI))
  ########################
  
  ################ avg trans per ring eq. 10-11
  avgtrans<-as.data.frame(t(colMeans(transmittance)))
  for (ring in 1:length(W)){
    colnames(avgtrans)[ring]<-paste0("AvgTrans_",ring)
  }
  ########################
  
  ################ Gap frac per ring eq. 10-12
  gaps<-as.data.frame(matrix(ncol=5,nrow=1))
  for (ring in 1:length(W)){
    colnames(gaps)[ring]<-paste0('Gap_',ring)
    gap_frac<-vector()
    for (i in 1:nrow(transmittance)){
      gap_frac[i]<-log(transmittance[i,ring])
      gaps[1,ring]<-exp(sum(gap_frac)*(1/length(gap_frac)))
    }
  }
  ########################
  
  ######################## clumping factor per ring eq. 10-16
  ACFS<-as.data.frame(log(avgtrans)/log(gaps))
  for (ring in 1:length(W)){
    colnames(ACFS)[ring]<-paste0("ACFS_",ring)
  }
  ########################
  
  ######################## total ACF eq.10.17
  LAIe<-2*sum(-(log(colMeans(transmittance[1:length(W)]))/(1/cos(angles_rad[1:length(W)])))*W)
  ACF<-LAIe/LAI
  ########################
  
  #estimate the footprint
  footprint_est<-tan(max(angles_rad[1:length(W)]))*mean(ABpairs$Height)
  
  LAI_results<-data.frame(LAI=LAI,LAIse=LAIse,LAIe=LAIe,ACF=ACF,N_samples=nrow(transmittance),footprint_m=footprint_est) # this might not work.... check
  
  results<-list(LAI_results,pointLAI,dodgy_reading)
  
  return(results)
}

Adj_A_from_match<-function(matchfile){
  a_indx<-seq(from=1,to=nrow(matchfile),by=2)
  b_indx<-seq(from=2,to=nrow(matchfile),by=2)
  A_adj<-colMeans(matchfile[b_indx,5:9]/matchfile[a_indx,5:9])
}

merge_AB<-function(B_readings,A_readings,A_adjust_vector,just_readings){
  if(is.character(A_readings$DateTime)){
    A_readings$DateTime<-as.POSIXct(A_readings$DateTime, format="%Y%m%d %H:%M:%S")
  }
  if(is.character(B_readings$DateTime)){
    B_readings$DateTime<-as.POSIXct(B_readings$DateTime, format="%Y%m%d %H:%M:%S")
  }
  if(just_readings==T){
    end_names<-c(colnames(B_readings),paste0("A_",colnames(A_readings)[5:9])) 
    ABpairs<-as.data.frame(matrix(ncol=ncol(B_readings)+5,nrow=0))    #create a blank matrix to add A readings to
    for (rw in 1:nrow(B_readings)){     #loop through each row
      idx <- which.min(abs(as.numeric(B_readings$DateTime[rw]) - as.numeric(A_readings$DateTime))) #find closest in time
      reading_pair<-cbind(B_readings[rw,],A_readings[idx,5:9]*A_adjust_vector) #bind with B_readings
      colnames(reading_pair)<-colnames(ABpairs) #remove colnames
      ABpairs<-rbind(ABpairs,reading_pair) #add to list
    }
    colnames(ABpairs)<-end_names
  }else{
    end_names<-c(colnames(B_readings),colnames(A_readings))
    ABpairs<-as.data.frame(matrix(ncol=ncol(B_readings)+ncol(A_readings),nrow=0))    #create a blank matrix to add A readings to
    for (rw in 1:nrow(B_readings)){     #loop through each row
      idx <- which.min(abs(as.numeric(B_readings$DateTime[rw]) - as.numeric(A_readings$DateTime))) #find closest in time
      reading_pair<-cbind(B_readings[rw,],A_readings[idx,]) #bind with B_readings
      colnames(reading_pair)<-colnames(ABpairs) #remove colnames
      ABpairs<-rbind(ABpairs,reading_pair) #add to list
      ABpairs[,(ncol(B_readings)+5):(ncol(B_readings)+9)]<-ABpairs[,(ncol(B_readings)+5):(ncol(B_readings)+9)]*A_adjust_vector
    }
    colnames(ABpairs)<-end_names         
  }
  return(ABpairs)
}

Align_AB_FV2200<-function(B_readings,A_readings,A_adjust_vector,K_readings,G_readings,text_head){
  
  # read in scatters
  Ks<-K_readings[,1:10]
  colnames(Ks)<-namevector
  
  # read in As - adjust from match file
  A_readings[,10]<-NA
  As<-A_readings[,1:10]
  As[,5:9]<-round(As[,5:9]*A_adjust_vector,2)
  colnames(As)<-namevector
  
  # read in Bs
  Bs<-B_readings[,1:10]
  Bs[,10]<-NA
  colnames(Bs)<-namevector
  
  Gs<-G_readings
  colnames(Gs)<-namevector
  
  for (colmn in 1:ncol(text_head)){
    text_head[,colmn]<-as.character(text_head[,colmn])
    Ks[,colmn]<-as.character(Ks[,colmn])
  }
  
  if(nrow(Ks)>0){
    Ks[,10]<-NA
    master<-rbind(text_head,Ks)
  }else{
    master<-text_head
  }
  
  for (rw in 1:nrow(Bs)){     #loop through each row
    idx <- which.min(abs(as.numeric(Bs$V3[rw]) - as.numeric(Gs$V3))) #find closest in time
    reading_pair<-rbind(Bs[rw,],Gs[idx,]) #bind with B_readings
    
    idx <- which.min(abs(as.numeric(Bs$V3[rw]) - as.numeric(As$V3))) #find closest in time
    reading_pair<-rbind(reading_pair,As[idx,]) #bind with B_readings
    
    for (colmn in 1:ncol(reading_pair)){
      reading_pair[,colmn]<-as.character(reading_pair[,colmn])
    }
    reading_pair$V3<-gsub("-","",reading_pair$V3,fixed=T)
    master<-rbind(master,reading_pair) #add to list
  }
  
  master
  
  
  return(master)
}

#####################################
### Read in files
#####################################
file_locations<-c("C:/Users/eck025/OneDrive - UWA/My Documents/PhD 2020/Field sites and data/LAI2200C/Results corrected K labelled/October 21/transects","C:/Users/eck025/OneDrive - UWA/My Documents/PhD 2020/Field sites and data/LAI2200C/Results corrected K labelled/March 22/transects")
transects<-Sys.glob(file.path(file_locations, "*.TXT"))
write_gps=FALSE

A_file_locations<-c("C:/Users/eck025/OneDrive - UWA/My Documents/PhD 2020/Field sites and data/LAI2200C/October 2021 results corrected/Transect A records","C:/Users/eck025/OneDrive - UWA/My Documents/PhD 2020/Field sites and data/LAI2200C/March 2022/Transect A readings")
A_transects<-Sys.glob(file.path(A_file_locations, "*.TXT"))

#create an object for adjusting A imports
matchfile<-read.delim("C:/Users/eck025/OneDrive - UWA/My Documents/PhD 2020/Field sites and data/LAI2200C/Match files 22/JE-45MTB matched.TXT", sep="\t",fill=TRUE,header=F,skip=36)
A_adj<-Adj_A_from_match(matchfile)



##############################################################################################
###### classification files ######
##############################################################################################
# read in veg composition map

projcrs <- st_crs("EPSG:4326") # set crs for points (use this rather than epsg strings)

vegmap<-st_read("data/Vegetation NVIS classification 15-3-23.shp")
colnames(vegmap)[c(3,5:(ncol(vegmap)-1))]<-c('ID','Orig_grouping','Understorey','Class_1','Height','proj_crown_cover%','MVG.Name','Key','MVSG.Name')
vegmap$MVG.Name<-as.factor(vegmap$MVG.Name)
vegmap$MVSG.Name<-as.factor(vegmap$MVSG.Name)
vegmap$Understorey<-as.factor(vegmap$Understorey)
vegmap[vegmap$MVG.Name=='Acacia Open Woodlands',]






# check out each class
table(vegmap$Understorey)
levels(vegmap$MVSG.Name)
table(vegmap$MVSG.Name)
levels(vegmap$MVG.Name)
table(vegmap$MVG.Name)


end_path<-"outputs"
dir.create(end_path,showWarnings = F)
LAI_results_destination<-paste0(end_path,"/reference_LAI.rds")


#########################################################
#########################################################
### Main loop
#########################################################
#########################################################

#create a bin for combined results
combined_results_mast<-list()
for (i in 1:length(transects)){
  # find the site and the directory for each file in the list
  site<-strsplit(tools::file_path_sans_ext(basename(transects[i])),"-")[[1]][1]
  wd<-dirname(transects[i])
  sample_period<-gsub(" ","",basename(dirname(wd)))
  
  ### This chunk is used to split A and B readings (or ground and canopy readings) and to remove extra GPS readings
  data <- read.delim(transects[i], sep="\t",fill=TRUE,header=F,skip=35)
  A_readings <- read.delim(A_transects[i], sep="\t",fill=TRUE,header=F,skip=35)
  
  ### read in the data head. it's a fuck for the text file this was from trial and error
  datahead <- read.delim(transects[i], sep="\t",fill=TRUE,header=F,nrows=21)
  datahead[5,2] #gives the site code
  namevector <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10")
  datahead[,namevector[3:10]] <- NA
  datahead2 <- read.delim(transects[i], sep="\t",fill=F,header=F,skip=21,nrows=8)
  datahead2[,namevector[7:10]] <- NA
  datahead3 <- read.delim(transects[i], sep="\t",fill=F,header=F,skip=30,nrows=1)
  datahead3[,namevector[2:10]] <- NA
  datahead4 <- read.delim(transects[i], sep="\t",fill=F,header=F,skip=31,nrows=1)
  datahead4[,namevector[9:10]] <- NA
  spacing1<-c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
  datahead5 <- read.delim(transects[i], sep="\t",fill=F,header=F,skip=34,nrows=1)
  datahead5[,namevector[2:10]] <- NA
  txthead<-rbind(datahead,datahead2,spacing1,datahead3,datahead4,spacing1,spacing1,datahead5)
  
  
  # name the columns as per the fv2200 files
  colnames(data)<-c("ABGK","Observation","DateTime","Type","Ring 1","Ring 2","Ring 3","Ring 4","Ring 5","DateTimeGPS")
  colnames(A_readings)<-c("ABGK","Observation","DateTime","Type","Ring 1","Ring 2","Ring 3","Ring 4","Ring 5")
  
  #change the data type
  data$Type<-as.factor(data$Type)
  data$ABGK<-as.factor(data$ABGK)
  
  
  Canopy<-data[!data$ABGK == 'B',]
  if(data$ABGK[1]=='K'){
    Understorey<-data[!data$ABGK == 'A',]
  } else {
    Understorey<-data[!data$ABGK == 'A',][-1,]
  }
  # Change As to Bs in the canopy file so it can be read into FV2200
  levels(Canopy$ABGK) <- c(levels(Canopy$ABGK), "B") 
  Canopy$ABGK[Canopy$ABGK == 'A'] <- 'B'
  
  ### remove GPS row if the previous row was also a gps point (deletes above or below gps data)
  Canopy<-rbind(Canopy[1,],Canopy[!Canopy$ABGK == data.table::shift(Canopy$ABGK, 1L, type="lag"),][-1,])
  Understorey<-rbind(Understorey[1,],Understorey[!Understorey$ABGK == data.table::shift(Understorey$ABGK, 1L, type="lag"),][-1,])
  
  ### change the Ks to As for FV2200 software
  Canopy$ABGK[Canopy$ABGK == 'K'] <- 'A'
  levels(Understorey$ABGK) <- c(levels(Understorey$ABGK), "A")
  Understorey$ABGK[Understorey$ABGK == 'K'] <- 'A'
  
  
  #to create a standalone gps file
  if (write_gps==TRUE){
    CanopyGPS<-Canopy[!Canopy$ABGK == 'A',]
    CanopyGPS<-CanopyGPS[!CanopyGPS$ABGK == 'B',]
    UnderstoreyGPS<-Understorey[!Understorey$ABGK == 'A',]
    UnderstoreyGPS<-UnderstoreyGPS[!UnderstoreyGPS$ABGK == 'B',]
    colnames(CanopyGPS)<-c("Type","Observation","DateTimeLocal","Type2","Lat","Long","Elevation","N satelites","Accuracy (m)","DateTimeGPS")
    colnames(UnderstoreyGPS)<-c("Type","Observation","DateTimeLocal","Type2","Lat","Long","Elevation","N satelites","Accuracy (m)","DateTimeGPS")
    
    write.csv(UnderstoreyGPS, file = paste0(wd,"/",site,"-G-transect-GPS-10-21.csv"),
              row.names = F,na="",quote=F)
    write.csv(CanopyGPS, file = paste0(wd,"/",site,"-C-transect-GPS-10-21.csv"),
              row.names = F,na="",quote=F)
  }
  
  
  
  
  ## 1. isolate gps from other data
  C1<-Canopy[Canopy$ABGK == 'G',]
  colnames(C1)<-c("Type","Observation","DateTimeLocal","Type2","Lat","Long","Elevation","N satelites","Accuracy (m)","DateTimeGPS")
  C2<-Canopy[!Canopy$ABGK == 'G',]
  G1<-Understorey[Understorey$ABGK == 'G',]
  colnames(G1)<-c("Type","Observation","DateTimeLocal","Type2","Lat","Long","Elevation","N satelites","Accuracy (m)","DateTimeGPS")
  G2<-Understorey[!Understorey$ABGK == 'G',]
  
  
  #2. recombine into one row
  C<-cbind(C2,C1)
  G<-cbind(G2,G1)
  C<-C[!C$Elevation==0,]
  G<-G[!G$Elevation==0,]
  
  
  #3. turn into spatial data
  ##### below is the correct process when using the sf package (sp is now deprecated)
  ##### attached is a good CRS summary https://geocompr.robinlovelace.net/reproj-geo-data.html
  
  C_sf<-st_as_sf(x = C,coords = c("Long", "Lat"),crs = projcrs) # convert points into spatial dataframe
  G_sf<-st_as_sf(x = G,coords = c("Long", "Lat"),crs = projcrs) # convert points into spatial dataframe
  
  
  ### 5. add structural info to points
  sf::sf_use_s2(FALSE) #### this tells sf not to use the s2 spherical coord system
                       #### ok for this application but problematic for doing geometric analyses
                       #### such as buffers and distances. in these situations reproject the data
                       #### to a planar coordinate ref system (currently WGS84)
                       #### see this ref https://geocompr.robinlovelace.net/reproj-geo-data.html
  CC <- data.frame(st_join(C_sf, vegmap))
  GG <- data.frame(st_join(G_sf, vegmap))
  
  CC$DateTime<-as.POSIXct(CC$DateTime, format="%Y%m%d %H:%M:%S")
  GG$DateTime<-as.POSIXct(GG$DateTime, format="%Y%m%d %H:%M:%S")
  A_readings$DateTime<-as.POSIXct(A_readings$DateTime, format="%Y%m%d %H:%M:%S")
  C1$DateTimeLocal<-as.POSIXct(C1$DateTimeLocal, format="%Y%m%d %H:%M:%S")
  G1$DateTimeLocal<-as.POSIXct(G1$DateTimeLocal, format="%Y%m%d %H:%M:%S")
  
  
  #6. create seperate points for scatter and readings. Could be done much earlier.... if (data$AGBK[c(1,3,5,7)]=='K'){Scatter_readings<-data[c(1,3,5,7),]}
  CCA<-CC[CC$ABGK=='A',]
  CCB<-CC[CC$ABGK=='B',]
  GGA<-GG[GG$ABGK=='A',]
  GGB<-GG[GG$ABGK=='B',]
  
  ### 7. Align a readings by closest in time
  Canopy_ABpairs<-merge_AB(CCB,A_readings,A_adj,just_readings=T)
  Ground_ABpairs<-merge_AB(GGB,A_readings,A_adj,just_readings=T)
  

  #8. isolate each vegetation class
  # remove non-existent levels
  C_MVG_list<-levels(droplevels(Canopy_ABpairs$MVG.Name))
  C_MVG_Num_list<-levels(droplevels(as.factor(Canopy_ABpairs$MVG)))
  # create a list of the ones remaining
  canopy_list <- list()
  for (p in 1:length(C_MVG_list)){
    canopy_list[[p]]<-Canopy_ABpairs[Canopy_ABpairs$MVG.Name == C_MVG_list[p],]
  }
  
  
  G_MVG_list<-levels(droplevels(Ground_ABpairs$MVG.Name))
  G_MVG_Num_list<-levels(droplevels(as.factor(Ground_ABpairs$MVG)))
  ground_list <- list()
  for (p in 1:length(G_MVG_list)){
    ground_list[[p]]<-Ground_ABpairs[Ground_ABpairs$MVG.Name == G_MVG_list[p],]
  }
  
  #### 9. apply results
  results_list <- list()
  Strata<-"Canopy"
  combined_results<-list()
  
  for (p in 1:length(canopy_list)){ # canopy loop start
    results_list[[p]]<-LAI_totals(canopy_list[[p]],Strata)
    idx<-which.max(as.data.frame(as.matrix(table(canopy_list[[p]]$Understorey)))[,1])
    Dom_under<-levels(canopy_list[[p]]$Understorey)[idx]
    combined_results[[p]]<-cbind(site,sample_period,Strata,C_MVG_list[[p]],C_MVG_Num_list[[p]],Dom_under,results_list[[p]][[1]])
    }
  combined_results<-do.call(rbind,combined_results)
  colnames(combined_results)<-c('Site','Sample_period','Strata','MVG','MVG_num','Dom_under',colnames(results_list[[p]][[1]]))
  
  
  #### as above for ground
  results_list <- list()
  Strata<-"Ground"
  combined_results.2<-list()
  for (p in 1:length(ground_list)){ # canopy loop start
    results_list[[p]]<-LAI_totals(ground_list[[p]],Strata)
    idx<-which.max(as.data.frame(as.matrix(table(ground_list[[p]]$Understorey)))[,1])
    Dom_under<-levels(ground_list[[p]]$Understorey)[idx]
    combined_results.2[[p]]<-cbind(site,sample_period,Strata,G_MVG_list[[p]],G_MVG_Num_list[[p]],Dom_under,results_list[[p]][[1]])
  }
  combined_results.2<-do.call(rbind,combined_results.2)
  colnames(combined_results.2)<-c('Site','Sample_period','Strata','MVG','MVG_num','Dom_under',colnames(results_list[[p]][[1]]))
  combined_results<-rbind(combined_results,combined_results.2)
  combined_results_mast[[i]]<-combined_results
  
  
  
  
  ############
  #############################################
  #############################################
  ### create files for FV2200
  #############################################
  #############################################
  Canopy_FV2200<-list()
  Ground_FV2200<-list()
  #8. isolate each vegetation class
  for (p in 1:length(C_MVG_list)){
    C<-Align_AB_FV2200(CCB[which(CCB$MVG.Name==C_MVG_list[p]),],A_readings,A_adj,CCA,C1,txthead)
    C[1,2]<-paste0(site,"_MVG-",gsub(' ','_',C_MVG_list[p],fixed=T),"_",sample_period,"_C")
    write.table(C, file = paste0(end_path,'/For_analysis_in_FV2200/',site,"_",C_MVG_list[p],"_",sample_period,"_C.txt"), sep = "\t",quote=F,
                row.names = F,col.names = F,na="",append=F)
  }
  for (p in 1:length(G_MVG_list)){
    G<-Align_AB_FV2200(GGB[which(GGB$MVG.Name==G_MVG_list[p]),],A_readings,A_adj,GGA,G1,txthead)
    G[1,2]<-paste0(site,"_MVG-",gsub(' ','_',G_MVG_list[p],fixed=T),"_",sample_period,"_G")
    write.table(G, file = paste0(end_path,'/For_analysis_in_FV2200/',site,"_",G_MVG_list[p],"_",sample_period,"_G.txt"), sep = "\t",quote=F,
                row.names = F,col.names = F,na="",append=F)
  }
} #main loop (files) end
combined_results<-do.call(rbind,combined_results_mast)
saveRDS(combined_results,LAI_results_destination)



















# Assuming 'data' is your dataset with columns 'LAI', 'MVG', and 'Replicates'
# Load necessary libraries
library(ggplot2)
library(ggpubr)

# Create a ggboxplot
y_position <- aggregate(combined_results[combined_results$Strata == 'Canopy',]$LAI, 
                        by = list(combined_results[combined_results$Strata == 'Canopy',]$MVG), FUN = max)
colnames(y_position)[1]<-'MVG'
combined_results$Sample_period<-as.factor(combined_results$Sample_period)
ggplot(combined_results[combined_results$Strata=='Canopy',], aes(x = MVG, y = ACF)) +
  geom_boxplot() +
  theme_minimal() +
  labs(y = "Canopy clumping",x="") +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1,size = 12)) + 
  scale_x_discrete(labels = function(x) gsub("[ ]", "\n", x))

mean(combined_results[combined_results$Strata=='Canopy',]$ACF)
mean(combined_results[combined_results$Strata=='Ground',]$ACF)


#### extra nonsense.
levels(vegmap$MVSG.Name)
vegmap[vegmap$MVSG.Name=='Melaleuca open forests and woodlands',] #add sedgelands to understorey list; remove CS from here...
vegmap[vegmap$MVSG.Name=='Melaleuca shrublands and open shrublands',] #shrub US
vegmap[vegmap$MVSG.Name=='Other Acacia forests and woodlands',] #Tussock US
vegmap[vegmap$MVSG.Name=='Other shrublands',] #Mostly Tussock US - maybe split...

vegmap$area<-st_area(vegmap)
aggregate(as.numeric(vegmap$area)/10000, list(vegmap$MVSG.Name), FUN=sum) #MVSG area in ha  




