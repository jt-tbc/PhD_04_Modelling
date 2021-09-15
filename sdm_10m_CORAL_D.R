# to determine number of cores on computer, this one is 4
parallel::detectCores()

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
setwd("D:/working/Shapefiles/25m/Tables")
GT_25m <- read.csv("01_GT_10mExtract_25m.csv")
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


# Let's read predictor variables (raster datasets)
lst <- list.files(path= "D:/working/GRD_Mask_Files/10m",pattern='.grd$',full.names = T) # list the name of files in the specified path,
lst 
preds <- stack(lst)
preds

#### CORAL MODELS ####

# Data preperation
d_25m <- sdmData(CORAL ~ depth + bscatter + slope + bpi25 + bpi3 + curv + hyp5 + vrm3 + rugosity,
                 train = North_25m, test = South_25m, predictors=preds)
d_25m

# Model Fitting
m10m_CORAL <- sdm( ~ depth + bscatter + slope + bpi25 + bpi3 + curv + hyp5 + vrm3 + rugosity,
                    data = d_25m,
                    methods = c('brt', 'rpart', 'glm', 'mars', 'maxent', 'rf'),
                    replication = c('cv'),
                    cv.folds = 5,
                    n = 10)
m10m_CORAL

# to save data and model products
setwd("D:/working/sdm_outputs_update/sdm_files")
write.sdm(d_25m,'sdm_course_10m_CORAL.sdd', overwrite = TRUE)
write.sdm(m10m_CORAL,'sdm_course_10m_CORAL.sdm', overwrite = TRUE)

# MODEL EVALUATION
meta_10m <- getModelInfo(m10m_CORAL)
ev_10m_test <- getEvaluation(m10m_CORAL, stat = c('AUC', 'COR', 'TSS', 'Sens', 'Spe','threshold'), opt = 2, wtest = 'test.dep')
ev_10m_eval <- getEvaluation(m10m_CORAL, stat = c('AUC', 'COR', 'TSS', 'Sens', 'Spe','threshold'), opt = 2, wtest = 'test.indep')

# need to extract bulk response curve information
vi_df <- as.data.frame(matrix(ncol = 9, nrow = 300))
colnames(vi_df)<-c("depth","bscatter", "slope","bpi25","bpi3","curv","hyp5", "vrm3", "rugosity")
model <- 1:300
for (i in model) {
  vi <- getVarImp(m10m_CORAL, id = i)
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
write.csv(merge3, "Evaluation_10m_CORAL.csv", row.names = F)

# PREDICT MODEL
setwd("D:/working/sdm_outputs_update/predictions")
predict(m10m_CORAL, preds, filename = 'predictions_10m_CORAL.img', mean = TRUE,
        parallelSettings = list(ncore=6,method='parallel'))

# ENSEMBLES
setwd("D:/working/sdm_outputs_update/predictions")
inf <- getModelInfo(m10m_CORAL)
ev <- getEvaluation(m10m_CORAL, stat = c('AUC', 'COR', 'TSS', 'Sens', 'Spe'), opt = 2)
inf_ev <- merge(inf, ev, all.x = T)

CORAL <- inf$modelID[inf$species == "CORAL"]

# build using all models
ensemble(m10m_CORAL, preds, filename = 'ensemble_10m_coral_all.img',
                        setting = list(id = CORAL, method = 'weighted', stat = 'AUC'), overwrite = TRUE)

# keep those with AUC > 0.75
id <- inf_ev$modelID[which(inf_ev$AUC >= 0.75 & inf_ev$species == "CORAL")]
e_AUC075_CORAL <- ensemble(m10m_CORAL, preds, filename = 'ensemble_10m_coral_AUC075.img',
                          setting = list(id = id, method = 'weighted', stat = 'AUC'), overwrite = TRUE)


gc()

setwd("D:/working/sdm_outputs_update/r_data")
# save.image(file = "sdm_10m_resolution.RData")
save.image(file = "sdm_10m_resolution_CORAL.RData")



####################
#### NULL MODEL ####
####################

# length of the raster
.d <- c(1,1,1,1,1,0,0,0,0,0,0,0)
# get .p from the data to determine probability of prescence, this is the prevelance
# use this to build a random layer
# you can then compare this surface to the generated layers
# wilcoxin test, non-parametric test, or paired t-test, check which test is appropriate
.p <- length(which(.d == 1)) / length(.d)
.pr <- rep(NA,12)
for (i in 1:12) .pr[i] <- rbinom(c(0,1),1,prob=c(1-.p,.p))
# generate p value about whether significantly different


install_github("danlwarren/ENMTools")
library(ENMTools)

library(devtools)
install_local("D:/Users/TUR262/Downloads/ENMTools-master.zip")
library(ENMTools)

