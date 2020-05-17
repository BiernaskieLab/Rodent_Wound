setwd("/home/ssinha/SCENIC_Trial/SCENIC_LWD14_fibro_transition")
library(SCENIC)
library(GENIE3)
library(RcisTarget)
library(AUCell)
scenicOptions <- readRDS("int/scenicOptions.Rds")

load("x.Robj")
exprMat <- x@raw.data
exprMat <- as.matrix(exprMat)
genesKept <- loadInt(scenicOptions, "genesKept")
exprMat_filtered <- exprMat[genesKept,]
exprMat_filtered <- log2(exprMat_filtered+1) 

runGenie3(exprMat_filtered, scenicOptions)
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat)
