##########################################
#### USE SDM PACKAGE TO CREATE MODELS ####
##########################################

# devtools::install_github('babaknaimi/sdm')
library(sdm)
library(raster)
library(maptools)
library(rgdal)
set.seed(3)

# 25 m for CORAL
setwd("C:/working/Shapefiles/25m/Tables")
GT_25m <- read.csv("02_GT_50mExtract_25m.csv")
GT_25m$bpi25 <- as.numeric(GT_25m$bpi25)
GT_25m$bpi3 <- as.numeric(GT_25m$bpi3)

South_25m <- subset(GT_25m, Area == "Osprey")
coordinates(South_25m)=~Long+Lat
proj4string(South_25m)<- CRS("+proj=longlat +datum=WGS84")
South_25m<-spTransform(South_25m,CRS("+proj=longlat"))

North_25m <- subset(GT_25m, Area != "Osprey")
coordinates(North_25m)=~Long+Lat
proj4string(North_25m)<- CRS("+proj=longlat +datum=WGS84")
North_25m<-spTransform(North_25m,CRS("+proj=longlat"))

class(North_25m) # it is a SpatialPointsDataFrame
# plot the data, P and A
plot(North_25m[North_25m$CORAL == 1,],col='blue',pch=16)
points(North_25m[North_25m$CORAL == 0,],col='red',pch=16)

class(South_25m) # it is a SpatialPointsDataFrame
# plot the data, P and A
plot(South_25m[South_25m$CORAL == 1,],col='blue',pch=16)
points(South_25m[South_25m$CORAL == 0,],col='red',pch=16)


# 100m m for MACRO and SPONGE
setwd("C:/working/Shapefiles/100m/Tables")
GT_100m <- read.csv("02_GT_50mExtract_100m.csv")
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
lst <- list.files(path= "C:/working/GRD_Mask_Files/50m/01_GRD_Files",pattern='.grd$',full.names = T) # list the name of files in the specified path,
lst 
preds <- stack(lst)
preds

#### CORAL MODELS ####

# Data preperation
d_25m <- sdmData(CORAL ~ depth + bscatter + slope + bpi25 + bpi3 + curv + hyp5 + vrm3 + rugosity,
                 train = North_25m, test = South_25m, predictors=preds)
d_25m

# Model Fitting
m50m_CORAL <- sdm( ~ depth + bscatter + slope + bpi25 + bpi3 + curv + hyp5 + vrm3 + rugosity,
                    data = d_25m,
                    methods = c('brt', 'rpart', 'glm', 'mars', 'maxent', 'rf'),
                    replication = c('cv'),
                    cv.folds = 5,
                    n = 10)
m50m_CORAL

# to save data and model products
setwd("C:/working/sdm_outputs_update/sdm_files")
write.sdm(d_25m,'sdm_course_50m_CORAL.sdd', overwrite = TRUE)
write.sdm(m50m_CORAL,'sdm_course_50m_CORAL.sdm', overwrite = TRUE)

# MODEL EVALUATION
meta_50m <- getModelInfo(m50m_CORAL)
ev_50m_test <- getEvaluation(m50m_CORAL, stat = c('AUC', 'COR', 'TSS', 'Sens', 'Spe','threshold'), opt = 2, wtest = 'test.dep')
ev_50m_eval <- getEvaluation(m50m_CORAL, stat = c('AUC', 'COR', 'TSS', 'Sens', 'Spe','threshold'), opt = 2, wtest = 'test.indep')

# need to extract bulk response curve information
vi_df <- as.data.frame(matrix(ncol = 9, nrow = 300))
colnames(vi_df)<-c("depth","bscatter", "slope","bpi25","bpi3","curv","hyp5", "vrm3", "rugosity")
model <- 1:300
for (i in model) {
  vi <- getVarImp(m50m_CORAL, id = i)
  if(is.null(vi) == FALSE) importance <- vi@varImportance$AUCtest else importance <- NA
  vi_df[i,] <- importance  
}
vi_df$modelID <- meta_50m$modelID

# combine and export evaluation data
merge <- merge(meta_50m, ev_50m_test, by.x = "modelID", by.y = "modelID", all.x = T)
merge2 <- merge(merge, ev_50m_eval, by.x = "modelID", by.y = "modelID", all.x = T)
merge3 <- merge(merge2, vi_df, by.x = "modelID", by.y = "modelID", all.x = T)
merge3$Scenario <- "50 m"
colnames(merge3) <- gsub(".x", ".test", colnames(merge3))
setwd("C:/working/sdm_outputs_update/evaluation")
write.csv(merge3, "Evaluation_50m_CORAL.csv", row.names = F)

# PREDICT MODEL
setwd("C:/working/sdm_outputs_update/predictions")
predict(m50m_CORAL, preds, filename = 'predictions_50m_CORAL.img', mean = TRUE)

# ENSEMBLES
setwd("C:/working/sdm_outputs_update/predictions")
inf <- getModelInfo(m50m_CORAL)
ev <- getEvaluation(m50m_CORAL, stat = c('AUC', 'COR', 'TSS', 'Sens', 'Spe'), opt = 2)
inf_ev <- merge(inf, ev, all.x = T)

CORAL <- inf$modelID[inf$species == "CORAL"]

# build using all models
e_all_CORAL <- ensemble(m50m_CORAL, preds, filename = 'ensemble_50m_coral_all.img',
                        setting = list(id = CORAL, method = 'weighted', stat = 'AUC'), overwrite = TRUE)

# keep those with AUC > 0.7, n = 227
id <- inf_ev$modelID[which(inf_ev$AUC >= 0.7 & inf_ev$species == "CORAL")]
e_AUC07_CORAL <- ensemble(m50m_CORAL, preds, filename = 'ensemble_50m_coral_AUC07.img',
                          setting = list(id = id, method = 'weighted', stat = 'AUC'), overwrite = TRUE)

gc()

#### MACRO MODELS ####

# Data preperation
d_100m <- sdmData(MACRO ~ depth + bscatter + slope + bpi25 + bpi3 + curv + hyp5 + vrm3 + rugosity,
                  train = North_100m, test = South_100m, predictors=preds)
d_100m

# Model Fitting
m50m_MACRO <- sdm( ~ depth + bscatter + slope + bpi25 + bpi3 + curv + hyp5 + vrm3 + rugosity,
                    data = d_100m,
                    methods = c('brt', 'rpart', 'glm', 'mars', 'maxent', 'rf'),
                    replication = c('cv'),
                    cv.folds = 5,
                    n = 10)
m50m_MACRO

# to save data and model products
setwd("C:/working/sdm_outputs_update/sdm_files")
write.sdm(d_100m,'sdm_course_50m_MACRO.sdd', overwrite = TRUE)
write.sdm(m50m_MACRO,'sdm_course_50m_MACRO.sdm', overwrite = TRUE)

# MODEL EVALUATION
meta_50m <- getModelInfo(m50m_MACRO)
ev_50m_test <- getEvaluation(m50m_MACRO, stat = c('AUC', 'COR', 'TSS', 'Sens', 'Spe','threshold'), opt = 2, wtest = 'test.dep')
ev_50m_eval <- getEvaluation(m50m_MACRO, stat = c('AUC', 'COR', 'TSS', 'Sens', 'Spe','threshold'), opt = 2, wtest = 'test.indep')

# need to extract bulk response curve information
vi_df <- as.data.frame(matrix(ncol = 9, nrow = 300))
colnames(vi_df)<-c("depth","bscatter", "slope","bpi25","bpi3","curv","hyp5", "vrm3", "rugosity")
model <- 1:300
for (i in model) {
  vi <- getVarImp(m50m_MACRO, id = i)
  if(is.null(vi) == FALSE) importance <- vi@varImportance$AUCtest else importance <- NA
  vi_df[i,] <- importance  
}
vi_df$modelID <- meta_50m$modelID

# combine and export evaluation data
merge <- merge(meta_50m, ev_50m_test, by.x = "modelID", by.y = "modelID", all.x = T)
merge2 <- merge(merge, ev_50m_eval, by.x = "modelID", by.y = "modelID", all.x = T)
merge3 <- merge(merge2, vi_df, by.x = "modelID", by.y = "modelID", all.x = T)
merge3$Scenario <- "50 m"
colnames(merge3) <- gsub(".x", ".test", colnames(merge3))
setwd("C:/working/sdm_outputs_update/evaluation")
write.csv(merge3, "Evaluation_50m_MACRO.csv", row.names = F)

# PREDICT MODEL
setwd("C:/working/sdm_outputs_update/predictions")
predict(m50m_MACRO, preds, filename = 'predictions_50m_MACRO.img', mean = TRUE)

# ENSEMBLES
setwd("C:/working/sdm_outputs_update/predictions")
inf <- getModelInfo(m50m_MACRO)
ev <- getEvaluation(m50m_MACRO, stat = c('AUC', 'COR', 'TSS', 'Sens', 'Spe'), opt = 2)
inf_ev <- merge(inf, ev, all.x = T)

MACRO <- inf$modelID[inf$species == "MACRO"]

# build using all models
e_all_MACRO <- ensemble(m50m_MACRO, preds, filename = 'ensemble_50m_MACRO_all.img',
                        setting = list(id = MACRO, method = 'weighted', stat = 'AUC'), overwrite = TRUE)

# keep those with AUC > 0.5
id2 <- inf_ev$modelID[which(inf_ev$AUC >= 0.5 & inf_ev$species == "MACRO")]
e_AUC07_MACRO <- ensemble(m50m_MACRO, preds, filename = 'ensemble_50m_MACRO_AUC05.img',
                          setting = list(id = id2, method = 'weighted', stat = 'AUC'), overwrite = TRUE)

gc()

#### SPONGE MODELS ####

# Data preperation
d_100mSP <- sdmData(SPONGE ~ depth + bscatter + slope + bpi25 + bpi3 + curv + hyp5 + vrm3 + rugosity,
                    train = North_100m, test = South_100m, predictors=preds)
d_100mSP

# Model Fitting
m50m_SPONGE <- sdm( ~ depth + bscatter + slope + bpi25 + bpi3 + curv + hyp5 + vrm3 + rugosity,
                     data = d_100mSP,
                     methods = c('brt', 'rpart', 'glm', 'mars', 'maxent', 'rf'),
                     replication = c('cv'),
                     cv.folds = 5,
                     n = 10)
m50m_SPONGE

# to save data and model products
setwd("C:/working/sdm_outputs_update/sdm_files")
write.sdm(d_100mSP,'sdm_course_50m_SPONGE.sdd', overwrite = TRUE)
write.sdm(m50m_SPONGE,'sdm_course_50m_SPONGE.sdm', overwrite = TRUE)

# MODEL EVALUATION
meta_50m <- getModelInfo(m50m_SPONGE)
ev_50m_test <- getEvaluation(m50m_SPONGE, stat = c('AUC', 'COR', 'TSS', 'Sens', 'Spe','threshold'), opt = 2, wtest = 'test.dep')
ev_50m_eval <- getEvaluation(m50m_SPONGE, stat = c('AUC', 'COR', 'TSS', 'Sens', 'Spe','threshold'), opt = 2, wtest = 'test.indep')

# need to extract bulk response curve information
vi_df <- as.data.frame(matrix(ncol = 9, nrow = 300))
colnames(vi_df)<-c("depth","bscatter", "slope","bpi25","bpi3","curv","hyp5", "vrm3", "rugosity")
model <- 1:300
for (i in model) {
  vi <- getVarImp(m50m_SPONGE, id = i)
  if(is.null(vi) == FALSE) importance <- vi@varImportance$AUCtest else importance <- NA
  vi_df[i,] <- importance  
}
vi_df$modelID <- meta_50m$modelID

# combine and export evaluation data
merge <- merge(meta_50m, ev_50m_test, by.x = "modelID", by.y = "modelID", all.x = T)
merge2 <- merge(merge, ev_50m_eval, by.x = "modelID", by.y = "modelID", all.x = T)
merge3 <- merge(merge2, vi_df, by.x = "modelID", by.y = "modelID", all.x = T)
merge3$Scenario <- "50 m"
colnames(merge3) <- gsub(".x", ".test", colnames(merge3))
setwd("C:/working/sdm_outputs_update/evaluation")
write.csv(merge3, "Evaluation_50m_SPONGE.csv", row.names = F)

# PREDICT MODEL
setwd("C:/working/sdm_outputs_update/predictions")
predict(m50m_SPONGE, preds, filename = 'predictions_50m_SPONGE.img', mean = TRUE)

# ENSEMBLES
setwd("C:/working/sdm_outputs_update/predictions")
inf <- getModelInfo(m50m_SPONGE)
ev <- getEvaluation(m50m_SPONGE, stat = c('AUC', 'COR', 'TSS', 'Sens', 'Spe'), opt = 2)
inf_ev <- merge(inf, ev, all.x = T)

SPONGE <- inf$modelID[inf$species == "SPONGE"]

# build using all models
e_all_SPONGE <- ensemble(m50m_SPONGE, preds, filename = 'ensemble_50m_SPONGE_all.img',
                         setting = list(id = SPONGE, method = 'weighted', stat = 'AUC'), overwrite = TRUE)

# keep those with AUC > 0.7
id3 <- inf_ev$modelID[which(inf_ev$AUC >= 0.7 & inf_ev$species == "SPONGE")]
e_AUC07_SPONGE <- ensemble(m50m_SPONGE, preds, filename = 'ensemble_50m_SPONGE_AUC07.img',
                           setting = list(id = id3, method = 'weighted', stat = 'AUC'), overwrite = TRUE)

setwd("C:/working/sdm_outputs_update/r_data")
# save.image(file = "sdm_50m_resolution.RData")
save.image(file = "sdm_50m_resolution.RData")
