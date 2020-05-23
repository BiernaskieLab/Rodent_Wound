
library(Seurat)

###### LWC_D14_fibr - upper lower
load(file = "LWC_D14_fibr.Robj")
save(LWC_D14_fibr_v3, file = "LWC_D14_fibr_v3.Robj")
LWC_D14_fibr_v3 = UpdateSeuratObject(LWC_D14_fibr)

LWC_D14_fibr_v3

View(LWC_D14_fibr_v3@meta.data)
DimPlot(LWC_D14_fibr_v3, reduction = "tsne")


resolutions <- seq(0, 1, 0.1)
LWC_D14_fibr_v3 <- FindNeighbors(LWC_D14_fibr_v3, dims = 1:10)
LWC_D14_fibr_v3 <- FindClusters(LWC_D14_fibr_v3, resolution = resolutions)

colnames(LWC_D14_fibr_v3@meta.data)

pdf("LWC_D14_fibr_v3_res_0_1.pdf", width = 10, height = 10)
DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0", reduction = "tsne")
DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.1", reduction = "tsne")
DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.2", reduction = "tsne")
DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.3", reduction = "tsne")
DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.4", reduction = "tsne")
DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.5", reduction = "tsne")
DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.6", reduction = "tsne")
DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.7", reduction = "tsne")
DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.8", reduction = "tsne")
DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.9", reduction = "tsne")
DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.1", reduction = "tsne")
dev.off()

View(LWC_D14_fibr_v3@meta.data)

install.packages("clustree")
library(clustree)
p12= clustree(LWC_D14_fibr_v3)

library(cowplot)
p1= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0", reduction = "tsne")
p2= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.1", reduction = "tsne")
p3= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.2", reduction = "tsne")
p4= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.3", reduction = "tsne")
p5= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.4", reduction = "tsne")
p6= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.5", reduction = "tsne")
p7= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.6", reduction = "tsne")
p8= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.7", reduction = "tsne")
p9= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.8", reduction = "tsne")
p10= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.9", reduction = "tsne")
p11= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.1", reduction = "tsne")

pdf("LWC_D14_fibr_v3_res_0_1_clustress.pdf", width = 40, height = 30)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
dev.off()


p1= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0", reduction = "tsne", label = T)
p2= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.1", reduction = "tsne", label = T)
p3= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.2", reduction = "tsne", label = T)
p4= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.3", reduction = "tsne", label = T)
p5= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.4", reduction = "tsne", label = T)
p6= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.5", reduction = "tsne", label = T)
p7= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.6", reduction = "tsne", label = T)
p8= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.7", reduction = "tsne", label = T)
p9= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.8", reduction = "tsne", label = T)
p10= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.9", reduction = "tsne", label = T)
p11= DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.1", reduction = "tsne", label = T)
p12= clustree(LWC_D14_fibr_v3)

pdf("LWC_D14_fibr_v3_res_0_1_clustress-LabelT.pdf", width = 40, height = 30)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
dev.off()


pdf("LWC_D14_fibr_v3_res_0_1_clustress-LabelT.pdf", width = 40, height = 30)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
dev.off()





pdf("LWC_D14_fibr_v3_res_0_1_clustress-LabelT-1.pdf", width = 15, height = 15)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11)
dev.off()

jpeg(file = "LWC_D14_fibr_v3_DimPlot_clust0.6.jpeg", width = 15, height = 15, units = "cm", res = 500)
DimPlot(LWC_D14_fibr_v3, group.by = "RNA_snn_res.0.6", reduction = "tsne", label = T)
dev.off()

pdf("LWC_D14_fibr_v3_res_0_1_clustress.pdf", width = 40, height = 30)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
dev.off()



###### LWD14fibr_subset - upper lower
load(file = "LWD14fibr_subset.Robj")
View(LWD14fibr_subset)

LWD14fibr_subset_v3 = UpdateSeuratObject(LWD14fibr_subset)
DimPlot(LWD14fibr_subset_v3, reduction = "tsne")

write.csv(LWD14fibr_subset_v3@meta.data, file = "LWD14fibr_subset_v3@meta.data.csv")
meta_data = read.csv(file = "LWD14fibr_subset_v3@meta.data.csv")
LWD14fibr_subset_v3$sample_ident = meta_data$sample_ident

DimPlot(LWD14fibr_subset_v3, reduction = "tsne", group.by = "sample_ident")

resolutions <- seq(0, 1, 0.1)
LWD14fibr_subset_v3 <- FindNeighbors(LWD14fibr_subset_v3, dims = 1:10)
LWD14fibr_subset_v3 <- FindClusters(LWD14fibr_subset_v3, resolution = resolutions)

p1= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0", reduction = "tsne")
p2= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.1", reduction = "tsne")
p3= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.2", reduction = "tsne")
p4= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.3", reduction = "tsne")
p5= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.4", reduction = "tsne")
p6= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.5", reduction = "tsne")
p7= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.6", reduction = "tsne")
p8= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.7", reduction = "tsne")
p9= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.8", reduction = "tsne")
p10= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.9", reduction = "tsne")
p11= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.1", reduction = "tsne")
p12= clustree(LWD14fibr_subset_v3)

pdf("LWD14fibr_subset_v3_clustress.pdf", width = 40, height = 30)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
dev.off()

pdf("LWD14fibr_subset_v3_clustress-1.pdf", width = 15, height = 15)
plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11)
dev.off()


LWD14fibr_subset_v3@active.ident

p0=DimPlot(LWD14fibr_subset_v3, reduction = "tsne", group.by = "sample_ident")
p1= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0", reduction = "tsne", label = T)
p2= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.1", reduction = "tsne", label = T)
p3= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.2", reduction = "tsne", label = T)
p4= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.3", reduction = "tsne", label = T)
p5= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.4", reduction = "tsne", label = T)
p6= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.5", reduction = "tsne", label = T)
p7= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.6", reduction = "tsne", label = T)
p8= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.7", reduction = "tsne", label = T)
p9= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.8", reduction = "tsne", label = T)
p10= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.9", reduction = "tsne", label = T)
p11= DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.1", reduction = "tsne", label = T)
p12= clustree(LWD14fibr_subset_v3)

pdf("LWD14fibr_subset_v3_clustress_labelT-1.pdf", width = 15, height = 20)
plot_grid(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
dev.off()

pdf("LWD14fibr_subset_v3_clustress_labelT-2.pdf", width = 15, height = 15)
plot_grid(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11)
dev.off()

View(LWD14fibr_subset_v3@meta.data)


jpeg(file = "LWD14fibr_subset_v3_DimPlot_clust0.6.jpeg", width = 15, height = 15, units = "cm", res = 500)
DimPlot(LWD14fibr_subset_v3, group.by = "RNA_snn_res.0.6", reduction = "tsne", label = T)
dev.off()


library(tidyverse)
LWD14fibr_subset_v3_quant = select(LWD14fibr_subset_v3@meta.data, sample_ident, RNA_snn_res.0.6)
LWD14fibr_subset_v3_quant$si_res = paste(LWD14fibr_subset_v3_quant$sample_ident, LWD14fibr_subset_v3_quant$RNA_snn_res.0.6, sep = "_")
LWD14fibr_subset_v3_quant$si_res = matrix(LWD14fibr_subset_v3_quant$si_res)
View(LWD14fibr_subset_v3_quant)
write.csv(LWD14fibr_subset_v3_quant$si_res, file = "LWD14fibr_subset_v3_quant$si_res.csv")

clust_quant <- read.table(pipe("pbpaste"), sep="\t", header = TRUE)
save(clust_quant, file = "clust_quant.csv")
clust_quant <- data.frame(clust_quant)

jpeg(file = "ggplot_clust_quant_LWD14.jpeg", width = 15, height = 15, units = "cm", res = 500)
ggplot(clust_quant, aes(fill=Sample, y=Value, x=Clust_ID)) + 
  geom_bar(position="fill", stat="identity") + scale_x_discrete(limits=c()) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) 
dev.off()

# Stacked + percent
jpeg(file = "ggplot_celltype_quantification.jpeg", width = 15, height = 15, units = "cm", res = 500)
ggplot(data_celltype_quantification, aes(fill=CellType, y=Value, x=Sample)) + 
  geom_bar(position="fill", stat="identity") + scale_x_discrete(limits=c("Cnt", "Los", "Ang", "Ang_Los")) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_fill_manual(values=c("#00bdc3", "#8d8cc4",
                             "#1eace3", "#1eb67b", 
                             "#eb68a6", "#a3a938",
                             "#bd77b2", "#3ab54a",
                             "#f17a6e", "#d99228"))
dev.off()







################################  All Wound
load(file = "all_fibros_4_scenic.Robj")

all_fibros_4_scenic_v3 = UpdateSeuratObject(all_fibros_4_scenic)
DimPlot(all_fibros_4_scenic_v3, reduction = "tsne")
View(all_fibros_4_scenic_v3@meta.data)

all_fibros_4_scenic_v3$sample.ident = all_fibros_4_scenic_v3@active.ident
FeaturePlot(all_fibros_4_scenic_v3, features = "Pdgfra")



load(file = "D14LWC_LWP_SW8_UI.Robj")
D14LWC_LWP_SW8_UI_v3 = UpdateSeuratObject(D14LWC_LWP_SW8_UI)

jpeg(file = "D14LWC_LWP_SW8_UI_v3_clustering_labelT.jpeg", width = 15, height = 15, units = "cm", res = 500)
DimPlot(D14LWC_LWP_SW8_UI_v3, reduction = "tsne", label = T)
dev.off()

View(D14LWC_LWP_SW8_UI_v3@active.ident)
View(D14LWC_LWP_SW8_UI_v3@meta.data)

write.csv(D14LWC_LWP_SW8_UI_v3@meta.data, file = "D14LWC_LWP_SW8_UI_v3@meta.data.csv")


sample_1 <- grep("-1", colnames(D14LWC_LWP_SW8_UI@data), value = T); length(sample_1)
sample_2 <- grep("-2", colnames(D14LWC_LWP_SW8_UI@data), value = T); length(sample_2)
sample_3 <- grep("-3", colnames(D14LWC_LWP_SW8_UI@data), value = T); length(sample_3)
sample_4 <- grep("-4", colnames(D14LWC_LWP_SW8_UI@data), value = T); length(sample_4)
sample_5 <- grep("-5", colnames(D14LWC_LWP_SW8_UI@data), value = T); length(sample_5)
sample_6 <- grep("-6", colnames(D14LWC_LWP_SW8_UI@data), value = T); length(sample_6)
sample_7 <- grep("-7", colnames(D14LWC_LWP_SW8_UI@data), value = T); length(sample_7)
Day18_Day8_Lib1_try1 <- grep("-8", colnames(D14LWC_LWP_SW8_UI@data), value = T); length(Day18_Day8_Lib1_try1)
Day18_Day8_Lib2_try1 <- grep("-9", colnames(D14LWC_LWP_SW8_UI@data), value = T); length(Day18_Day8_Lib2_try1)
Day18_Day8_Lib3_try1 <- grep("-10", colnames(D14LWC_LWP_SW8_UI@data), value = T); length(Day18_Day8_Lib3_try1)
Day18_Day8_Lib4_try1 <- grep("-11", colnames(D14LWC_LWP_SW8_UI@data), value = T); length(Day18_Day8_Lib4_try1)




########### LWP AUC Signatures ########### 
pdf("AUCell_plotTSNE_Thap11.pdf")
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Thap11")],], 
                        plots="Expression",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Thap11")],], 
                        plots="AUC",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)

AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Thap11")],], 
                        plots="Binary",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
dev.off()



pdf("AUCell_plotTSNE_Gtf3a.pdf")
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Gtf3a")],], 
                        plots="Expression",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Gtf3a")],], 
                        plots="AUC",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)

AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Gtf3a")],], 
                        plots="Binary",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
dev.off()


pdf("AUCell_plotTSNE_Rara.pdf")
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Rara")],], 
                        plots="Expression",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Rara")],], 
                        plots="AUC",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)

AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Rara")],], 
                        plots="Binary",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
dev.off()


pdf("AUCell_plotTSNE_Zfp354c.pdf")
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Zfp354c")],], 
                        plots="Expression",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Zfp354c")],], 
                        plots="AUC",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)

AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Zfp354c")],], 
                        plots="Binary",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
dev.off()

pdf("AUCell_plotTSNE_Nfib.pdf")
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Nfib")],], 
                        plots="Expression",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Nfib")],], 
                        plots="AUC",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)

AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Nfib")],], 
                        plots="Binary",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
dev.off()


########### LWC AUC Signatures ########### 

pdf("AUCell_plotTSNE_Cd59a.pdf")
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Cd59a")],], 
                        plots="Expression",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Cd59a")],], 
                        plots="AUC",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)

AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Cd59a")],], 
                        plots="Binary",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
dev.off()

pdf("AUCell_plotTSNE_Nfix.pdf")
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Nfix")],], 
                        plots="Expression",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Nfix")],], 
                        plots="AUC",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)

AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Nfix")],], 
                        plots="Binary",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
dev.off()

pdf("AUCell_plotTSNE_Prrx1.pdf")
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Prrx1")],], 
                        plots="Expression",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Prrx1")],], 
                        plots="AUC",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)

AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Prrx1")],], 
                        plots="Binary",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
dev.off()

pdf("AUCell_plotTSNE_Osr1.pdf")
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Osr1")],], 
                        plots="Expression",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Osr1")],], 
                        plots="AUC",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)

AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Osr1")],], 
                        plots="Binary",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
dev.off()

pdf("AUCell_plotTSNE_Osr1.pdf")
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Osr1")],], 
                        plots="Expression",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Osr1")],], 
                        plots="AUC",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)

AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Osr1")],], 
                        plots="Binary",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
dev.off()

pdf("AUCell_plotTSNE_Osr1.pdf")
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Osr1")],], 
                        plots="Expression",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Osr1")],], 
                        plots="AUC",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)

AUCell::AUCell_plotTSNE(tSNE_scenic$Y, 
                        exprMat, 
                        aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Osr1")],], 
                        plots="Binary",
                        borderColor = adjustcolor("lightgray", alpha.f = 0.8), alphaOff = 0.5, alphaOn = 1, cex = 1.1)
dev.off()







