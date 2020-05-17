library(SCENIC)
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$seed <- 123

tsneCompare <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(5,15,50), perpl=c(5,15,50))
tsneCompare <- tsneAUC(scenicOptions, aucType="AUC", nPcs=c(5,15,50), perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/): 
tsneCompare <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
plotTsne_compareSettings(tsneCompare, scenicOptions, showLegend=FALSE)

par(mfcol=c(3,3))
tsneCompare_high <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(tsneCompare_high, scenicOptions, showLegend=FALSE, varName="CellType", cex=.5)

#scenicOptions@settings$defaultTsne$aucType <- "AUC"
#scenicOptions@settings$defaultTsne$dims <- 50
#scenicOptions@settings$defaultTsne$perpl <- 50