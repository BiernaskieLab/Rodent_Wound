

library(RColorBrewer)

jpeg(file = "SCENIC_Contour_1.jpeg", width = 20, height = 20, units = "cm", res = 500)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(10, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=8, drawlabels=F)
dev.off()


jpeg(file = "SCENIC_Creb3l1_Gli1_colocalize_1.jpeg", width = 20, height = 20, units = "cm", res = 500)
regulonNames <- c( "Creb3l1","Gli1")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
text(-30,-25, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(-30,-25-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)
dev.off()

jpeg(file = "Hox_genes.jpeg", width = 20, height = 20, units = "cm", res = 500)
regulonNames <- list(red=c("Hoxc5", ""),
                     green=c("Hoxc10"),
                     blue=c( "Hoxd8"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC")
text(-30,-25, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(-30,-25-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)
text(-30,-25-8, attr(cellCol,"blue"), col="blue", cex=.7, pos=4)
dev.off()

jpeg(file = "Hox_genes_trial_1col_hoxc.jpeg", width = 18, height = 18, units = "cm", res = 500)
regulonNames <- list(red=c("Hoxc5","Hoxc10"))
                     #green=c("Hoxa5", "Hoxa9"),
                     #blue=c("Hoxd8"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.01)
text(-30,-25, attr(cellCol,"red"), col="tomato", cex=1, pos=4)
#text(-30,-25-4, attr(cellCol,"green"), col="green", cex=1, pos=4)
#text(-30,-25-8, attr(cellCol,"cornflowerblue"), col="blue", cex=1, pos=4)
dev.off()


jpeg(file = "Hox_genes_trial_1col_hoxa.jpeg", width = 18, height = 18, units = "cm", res = 500)
regulonNames <- list(#red=c("Hoxc5","Hoxc10"),
                     green=c("Hoxa5", "Hoxa9"))
#blue=c("Hoxd8"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.01)
#text(-30,-25, attr(cellCol,"red"), col="tomato", cex=1, pos=4)
text(-30,-25-4, attr(cellCol,"green"), col="green", cex=1, pos=4)
#text(-30,-25-8, attr(cellCol,"cornflowerblue"), col="blue", cex=1, pos=4)
dev.off()

jpeg(file = "Hox_genes_trial_1col_hoxd.jpeg", width = 18, height = 18, units = "cm", res = 500)
regulonNames <- list(blue=c("Hoxd8"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.01)
#text(-30,-25, attr(cellCol,"red"), col="tomato", cex=1, pos=4)
#text(-30,-25-4, attr(cellCol,"green"), col="green", cex=1, pos=4)
text(-30,-25-8, attr(cellCol,"cornflowerblue"), col="blue", cex=1, pos=4)
dev.off()

jpeg(file = "Rarg_regulon.jpeg", width = 18, height = 18, units = "cm", res = 500)
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[(rownames(aucell_regulonAUC))[c("Rarg")],], plots="Binary", alphaOn = 5, borderColor = adjustcolor("lightgray", alpha.f = 0.8))
dev.off()

jpeg(file = "Rarg_regulon.jpeg", width = 18, height = 18, units = "cm", res = 500)
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[(rownames(aucell_regulonAUC))[c("Rarg")],], plots="Binary", alphaOn = 5, borderColor = adjustcolor("lightgray", alpha.f = 0.8))
dev.off()

head(tSNE_scenic$Y)


jpeg(file = "Rar_regulon.jpeg", width = 18, height = 18, units = "cm", res = 500)
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Rarg")],], plots="Binary", alphaOn = 5, borderColor = adjustcolor("lightgray", alpha.f = 0.8))
dev.off()


regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))


pdf("Step4_BinaryRegulonActivity_tSNE_colByAUC.pdf")
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="Binary")
dev.off()




jpeg(file = "3_stage_seperation.jpeg", width = 20, height = 20, units = "cm", res = 500)
regulonNames <- list(red=c("Tead1"),
                     green=c("Nfib"),
                     blue=c( "Mef2c"))
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="Binary")
text(-30,-25, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(-30,-25-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)
text(-30,-25-8, attr(cellCol,"blue"), col="blue", cex=.7, pos=4)
dev.off()





scenicOptions <- readRDS("scenicOptions.Rds")
