# to determine number of cores on computer, this one is 4
parallel::detectCores()

##########################################
#### USE SDM PACKAGE TO CREATE MODELS ####
##########################################

rm(list = ls())
gc()

# devtools::install_github('babaknaimi/sdm')
library(sdm)
library(raster)
library(maptools)
library(rgdal)
set.seed(3)

# 100 m for MACRO and SPONGE
setwd("D:/working/Shapefiles/100m/Tables")
GT_100m <- read.csv("01_GT_10mExtract_100m.csv")
GT_100m$bpi25 <- as.numeric(GT_100m$bpi25)
GT_100m$bpi3 <- as.numeric(GT_100m$bpi3)

South_100m <- subset(GT_100m, Area == "Osprey")
coordinates(South_100m)=~Long+Lat
proj4string(South_100m)<- CRS("+proj=longlat +datum=WGS84")
South_100m<-spTransform(South_100m,CRS("+proj=longlat"))

North_100m <- subset(GT_100m, Area != "Osprey")
coordinates(North_100m)=~Long+Lat
proj4string(North_100m)<- CRS("+proj=longlat +datum=WGS84")
North_100m<-spTransform(North_100m,CRS("+proj=longlat"))

class(North_100m) # it is a SpatialPointsDataFrame
# plot the data, P and A
plot(North_100m[North_100m$CORAL == 1,],col='blue',pch=16)
points(North_100m[North_100m$CORAL == 0,],col='red',pch=16)

class(South_100m) # it is a SpatialPointsDataFrame
# plot the data, P and A
plot(South_100m[South_100m$CORAL == 1,],col='blue',pch=16)
points(South_100m[South_100m$CORAL == 0,],col='red',pch=16)

# Let's read predictor variables (raster datasets)
lst <- list.files(path= "D:/working/GRD_Mask_Files/10m",pattern='.grd$',full.names = T) # list the name of files in the specified path,
lst 
preds <- stack(lst)
preds


#### SPONGE MODELS ####

# Data preperation
d_10mSP <- sdmData(SPONGE ~ depth + bscatter + slope + bpi25 + bpi3 + curv + hyp5 + vrm3 + rugosity,
                   train = North_100m, test = South_100m, predictors=preds)
d_10mSP

# Model Fitting
m10m_SPONGE <- sdm( ~ depth + bscatter + slope + bpi25 + bpi3 + curv + hyp5 + vrm3 + rugosity,
                    data = d_10mSP,
                    methods = c('brt', 'rpart', 'glm', 'mars', 'maxent', 'rf'),
                    replication = c('cv'),
                    cv.folds = 5,
                    n = 10)
m10m_SPONGE

# to save data and model products
setwd("D:/working/sdm_outputs_update/sdm_files")
write.sdm(d_10mSP,'sdm_course_10m_SPONGE.sdd', overwrite = TRUE)
write.sdm(m10m_SPONGE,'sdm_course_10m_SPONGE.sdm', overwrite = TRUE)

# MODEL EVALUATION
meta_10m <- getModelInfo(m10m_SPONGE)
ev_10m_test <- getEvaluation(m10m_SPONGE, stat = c('AUC', 'COR', 'TSS', 'Sens', 'Spe','threshold'), opt = 2, wtest = 'test.dep')
ev_10m_eval <- getEvaluation(m10m_SPONGE, stat = c('AUC', 'COR', 'TSS', 'Sens', 'Spe','threshold'), opt = 2, wtest = 'test.indep')

# need to extract bulk response curve information
vi_df <- as.data.frame(matrix(ncol = 9, nrow = 300))
colnames(vi_df)<-c("depth","bscatter", "slope","bpi25","bpi3","curv","hyp5", "vrm3", "rugosity")
model <- 1:300
for (i in model) {
  vi <- getVarImp(m10m_SPONGE, id = i)
  if(is.null(vi) == FALSE) importance <- vi@varImportance$AUCtest else importance <- NA
  vi_df[i,] <- importance  
}
vi_df$modelID <- meta_10m$modelID

# combine and export evaluation data
merge <- merge(meta_10m, ev_10m_test, by.x = "modelID", by.y = "modelID", all.x = T)
merge2 <- merge(merge, ev_10m_eval, by.x = "modelID", by.y = "modelID", all.x = T)
merge3 <- merge(merge2, vi_df, by.x = "modelID", by.y = "modelID", all.x = T)
merge3$Scenario <- "10 m"
colnames(merge3) <- gsub(".x", ".test", colnames(merge3))
setwd("D:/working/sdm_outputs_update/evaluation")
write.csv(merge3, "Evaluation_10m_SPONGE.csv", row.names = F)

# PREDICT MODEL
setwd("D:/working/sdm_outputs_update/predictions")
predict(m10m_SPONGE, preds, filename = 'predictions_10m_SPONGE.img', mean = TRUE,
        parallelSettings = list(ncore=6,method='parallel'))

# ENSEMBLES
setwd("D:/working/sdm_outputs_update/predictions")
inf <- getModelInfo(m10m_SPONGE)
ev <- getEvaluation(m10m_SPONGE, stat = c('AUC', 'COR', 'TSS', 'Sens', 'Spe'), opt = 2)
inf_ev <- merge(inf, ev, all.x = T)

SPONGE <- inf$modelID[inf$species == "SPONGE"]

# build using all models
e_all_SPONGE <- ensemble(m10m_SPONGE, preds, filename = 'ensemble_10m_SPONGE_all.img',
                         setting = list(id = SPONGE, method = 'weighted', stat = 'AUC'), overwrite = TRUE)

# keep those with AUC > 0.7
id3 <- inf_ev$modelID[which(inf_ev$AUC >= 0.7 & inf_ev$species == "SPONGE")]
e_AUC07_SPONGE <- ensemble(m10m_SPONGE, preds, filename = 'ensemble_10m_SPONGE_AUC07.img',
                           setting = list(id = id3, method = 'weighted', stat = 'AUC'), overwrite = TRUE)

gc()

setwd("D:/working/sdm_outputs_update/r_data")
# save.image(file = "sdm_10m_resolution.RData")
save.image(file = "sdm_10m_resolution_SPONGE.RData")
