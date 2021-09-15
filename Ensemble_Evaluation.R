library(sdm)

setwd("C:/Users/TUR262/Documents/Chapters 4 and 5 Update/Resolution_ReAnalysis/Data/Ensembles")

# CORAL ensembles
lst_CORAL <- list.files(path= "C:/Users/TUR262/Documents/Chapters 4 and 5 Update/Resolution_ReAnalysis/Data/Ensembles",
                  pattern='Coral')

number <- c(1:8)
evaluation_CORAL <- as.data.frame(matrix(ncol = 7, nrow = 0))
colnames(evaluation_CORAL)<-c("csv","AUC","COR","TSS","sensitivity","specificity","threshold")

for (i in number) {
  df <- read.csv(lst_CORAL[i], header = T)
  e <- evaluates(df$CORAL, df$predicted)
  evaluation_CORAL[i,1] <- lst_CORAL[i] 
  evaluation_CORAL[i,2] <- e@statistics$AUC
  evaluation_CORAL[i,3] <- e@statistics$COR[1]
  evaluation_CORAL[i,4] <- e@threshold_based$TSS[2]
  evaluation_CORAL[i,5] <- e@threshold_based$sensitivity[2]
  evaluation_CORAL[i,6] <- e@threshold_based$specificity[2]
  evaluation_CORAL[i,7] <- e@threshold_based$threshold[2]
}

# MACRO ensembles
lst_MACRO <- list.files(path= "C:/Users/TUR262/Documents/Chapters 4 and 5 Update/Resolution_ReAnalysis/Data/Ensembles",
                        pattern='Macro')

number <- c(1:8)
evaluation_MACRO <- as.data.frame(matrix(ncol = 7, nrow = 0))
colnames(evaluation_MACRO)<-c("csv","AUC","COR","TSS","sensitivity","specificity","threshold")

for (i in number) {
  df <- read.csv(lst_MACRO[i], header = T)
  e <- evaluates(df$MACRO, df$predicted)
  evaluation_MACRO[i,1] <- lst_MACRO[i] 
  evaluation_MACRO[i,2] <- e@statistics$AUC
  evaluation_MACRO[i,3] <- e@statistics$COR[1]
  evaluation_MACRO[i,4] <- e@threshold_based$TSS[2]
  evaluation_MACRO[i,5] <- e@threshold_based$sensitivity[2]
  evaluation_MACRO[i,6] <- e@threshold_based$specificity[2]
  evaluation_MACRO[i,7] <- e@threshold_based$threshold[2]
}


# SPONGE ensembles
lst_SPONGE <- list.files(path= "C:/Users/TUR262/Documents/Chapters 4 and 5 Update/Resolution_ReAnalysis/Data/Ensembles",
                        pattern='Sponge')

number <- c(1:8)
evaluation_SPONGE <- as.data.frame(matrix(ncol = 7, nrow = 0))
colnames(evaluation_SPONGE)<-c("csv","AUC","COR","TSS","sensitivity","specificity","threshold")

for (i in number) {
  df <- read.csv(lst_SPONGE[i], header = T)
  e <- evaluates(df$SPONGE, df$predicted)
  evaluation_SPONGE[i,1] <- lst_SPONGE[i] 
  evaluation_SPONGE[i,2] <- e@statistics$AUC
  evaluation_SPONGE[i,3] <- e@statistics$COR[1]
  evaluation_SPONGE[i,4] <- e@threshold_based$TSS[2]
  evaluation_SPONGE[i,5] <- e@threshold_based$sensitivity[2]
  evaluation_SPONGE[i,6] <- e@threshold_based$specificity[2]
  evaluation_SPONGE[i,7] <- e@threshold_based$threshold[2]
}


final_evaluation <- rbind(evaluation_CORAL,evaluation_MACRO,evaluation_SPONGE)
setwd("C:/Users/TUR262/Documents/Chapters 4 and 5 Update/Resolution_ReAnalysis/Data")
write.csv(final_evaluation,"01_Ensemble_Evaluation.csv", row.names = F)
