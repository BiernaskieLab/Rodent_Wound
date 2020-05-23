# Assigning sample #s:
sample_2 <- grep("-2", colnames(macro_recluster@data), value = T); length(sample_2)
sample_4 <- grep("-4", colnames(macro_recluster@data), value = T); length(sample_4)
sample_6 <- grep("-6", colnames(macro_recluster@data), value = T); length(sample_6)
sample_11 <- grep("-11", colnames(macro_recluster@data), value = T); length(sample_11)


########## Macrophase from wound_neg_0.3 and Re-Clustering ######
macros <- names(wound_neg_0.3@ident[wound_neg_0.3@ident %in% c("Macrophages")])
sub.macros <- FilterCells(wound_neg_0.3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = macros)
macro_recluster <- FilterCells(sub.macros, cells.use = c(sample_2, sample_4, sample_6, sample_11), subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
macro_recluster <- FindVariableGenes(macro_recluster, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
length(x = macro_recluster@var.genes)
macro_recluster <- NormalizeData(macro_recluster, normalization.method = "LogNormalize", scale.factor = 10000)
macro_recluster <- ScaleData(macro_recluster, vars.to.regress = c("nUMI", "percent.mito"))
macro_recluster <- RunPCA(macro_recluster, pc.genes = macro_recluster@var.genes, pcs.compute = 50)
PCElbowPlot(macro_recluster, num.pc = 50) ## 22 significant PCs
macro_recluster <- FindClusters(macro_recluster, reduction.type = "pca", dims.use = 1:22, resolution = 1.5, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
macro_recluster <- StashIdent(macro_recluster, save.name = "ClusterNames_macro_recluster")
macro_recluster <- RunTSNE(macro_recluster, dims.use = 1:22, do.fast = T)
TSNEPlot(macro_recluster, do.label = TRUE, pt.size = 0.8)


# Looking at T Cell Clusters by sample
jpeg(file = "Macro_seperated_by_sample_LW.jpeg", width = 15, height = 15, units = "cm", res = 500)
or_tSNE <- as.data.frame(GetCellEmbeddings(macro_recluster_LW, reduction.type = "tsne"))
or_tSNE[sample_2, "Sample"] <- "D14 LWC"
or_tSNE[sample_4, "Sample"] <- "D14 LWP"
or_tSNE[sample_6, "Sample"] <- "D14 SW"
or_tSNE[sample_11, "Sample"] <- "D8 SW"
ggplot(or_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample), size=1.00) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()



# LW ONLY Macrophases
macro_recluster_LW <- FilterCells(macro_recluster, cells.use = c(sample_2, sample_4), subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
macro_recluster_LW <- FindVariableGenes(macro_recluster_LW, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
length(x = macro_recluster_LW@var.genes)
macro_recluster_LW <- NormalizeData(macro_recluster_LW, normalization.method = "LogNormalize", scale.factor = 10000)
macro_recluster_LW <- ScaleData(macro_recluster_LW, vars.to.regress = c("nUMI", "percent.mito"))

macro_recluster_LW <- RunPCA(macro_recluster_LW, pc.genes = macro_recluster_LW@var.genes, pcs.compute = 50)
PCElbowPlot(macro_recluster_LW, num.pc = 50) ## 15 significant PCs
macro_recluster_LW <- FindClusters(macro_recluster_LW, reduction.type = "pca", dims.use = 1:15, resolution = 1.5, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
macro_recluster_LW <- StashIdent(macro_recluster_LW, save.name = "ClusterNames_macro_recluster_LW")
macro_recluster_LW <- RunTSNE(macro_recluster_LW, dims.use = 1:15, do.fast = T)
TSNEPlot(macro_recluster_LW, do.label = TRUE, pt.size = 0.8)

write.csv(macro_recluster_LW@meta.data, "Cluster5_composition.csv")

save(macro_recluster_LW, file = "macro_recluster_LW.Robj")
load(file = "macro_recluster.Robj")
load(file = "macro_recluster_LW.Robj")


jpeg(file = "Macro_seperated_by_sample.jpeg", width = 15, height = 15, units = "cm", res = 500)
or_tSNE <- as.data.frame(GetCellEmbeddings(macro_recluster_LW, reduction.type = "tsne"))
or_tSNE[sample_2, "Sample"] <- "D14 LWC"
or_tSNE[sample_4, "Sample"] <- "D14 LWP"
or_tSNE[sample_6, "Sample"] <- "D14 SW"
or_tSNE[sample_11, "Sample"] <- "D8 SW"
ggplot(or_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample), size=0.80) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()

pdf("macro_LW_reclustering.pdf")
ggplot(or_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample), size=0.80) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
FeaturePlot(macro_recluster_LW, features.plot = c("Cx3cr1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(macro_recluster_LW, features.plot = c("Ccr2"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(macro_recluster_LW, features.plot = c("Cxcl1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(macro_recluster_LW, features.plot = c("Cxcl2"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(macro_recluster_LW, features.plot = c("Cx3cr1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()

FeaturePlot(macro_recluster, features.plot = c("Malat1"), reduction.use = "tsne", cols.use = c("grey", "blue"))


FeaturePlot(wound_neg_, features.plot = c("Lta4h"), reduction.use = "tsne", cols.use = c("grey", "blue"))

FeaturePlot(macro_recluster_LW, features.plot = c("Cxcl12"), reduction.use = "tsne", cols.use = c("grey", "blue"))


regen_macro.markers <- FindMarkers(macro_recluster_LW, ident.1 = c(5), min.pct = 0.20, test.use = "negbinom")
regen_macro.markers <- regen_macro.markers[order(regen_macro.markers$avg_logFC),]
write.csv(regen_macro.markers,file = "regen_macro.markers.csv")


LW_macro_4.markers <- FindMarkers(macro_recluster_LW, ident.1 = c(4), min.pct = 0.20, test.use = "negbinom")
LW_macro_4.markers <- LW_macro_4.markers[order(LW_macro_4.markers$avg_logFC),]
write.csv(LW_macro_4.markers,file = "LW_macro_4.markers.csv")


# Cluster composition of macro_recluster_LW:
##Cluster 5 = Regen cluster - 21 LWO; 94 LWC
##Cluster 4 = Regen cluster - 47 LWO; 78 LWC

##Cluster 2 = Pheriphery cluster - 85 LWO; 56 LWC
##Cluster 11 = Pheriphery cluster -  LWO;  LWC

load(file = "wound_subset.Robj")
sample_1 <- grep("-1", colnames(wound_subset@data), value = T); length(sample_1)
sample_2 <- grep("-2", colnames(wound_subset@data), value = T); length(sample_2)
sample_3 <- grep("-3", colnames(wound_subset@data), value = T); length(sample_3)
sample_4 <- grep("-4", colnames(wound_subset@data), value = T); length(sample_4)
sample_5 <- grep("-5", colnames(wound_subset@data), value = T); length(sample_5)
sample_6 <- grep("-6", colnames(wound_subset@data), value = T); length(sample_6)
sample_7 <- grep("-7", colnames(wound_subset@data), value = T); length(sample_7)
Day18_Day8_Lib1_try1 <- grep("-8", colnames(wound_subset@data), value = T); length(Day18_Day8_Lib1_try1)
Day18_Day8_Lib2_try1 <- grep("-9", colnames(wound_subset@data), value = T); length(Day18_Day8_Lib2_try1)
Day18_Day8_Lib3_try1 <- grep("-10", colnames(wound_subset@data), value = T); length(Day18_Day8_Lib3_try1)
Day18_Day8_Lib4_try1 <- grep("-11", colnames(wound_subset@data), value = T); length(Day18_Day8_Lib4_try1)


jpeg(file = "wound_subset_macrophase_Dpt.jpeg", width = 15, height = 15, units = "cm", res = 500)
FeaturePlot(wound_subset, features.plot = c("Dpt"), reduction.use = "tsne", cols.use = c("yellow", "red"))
dev.off()

epi_mac_fibr <- names(wound_subset@ident[wound_subset@ident %in% c(11,2, 22,15, 21, 6, 10, 1,3,0,13,5,12,20,9,7)])
sub.epi_mac_fibr <- FilterCells(wound_subset, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = epi_mac_fibr)

epi_mac_fibr_SW8 <- FilterCells(sub.epi_mac_fibr, cells.use = c(Day18_Day8_Lib3_try1, Day18_Day8_Lib4_try1), subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
epi_mac_fibr_LW14 <- FilterCells(sub.epi_mac_fibr, cells.use = c(sample_1, sample_2, sample_3, sample_4), subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)


# epi_mac_fibr_SW8 Reclustering:
epi_mac_fibr_SW8 <- FindVariableGenes(epi_mac_fibr_SW8, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
length(x = epi_mac_fibr_SW8@var.genes) #1334
epi_mac_fibr_SW8 <- NormalizeData(epi_mac_fibr_SW8, normalization.method = "LogNormalize", scale.factor = 10000)
epi_mac_fibr_SW8 <- ScaleData(epi_mac_fibr_SW8, vars.to.regress = c("nUMI", "percent.mito"))
epi_mac_fibr_SW8 <- RunPCA(epi_mac_fibr_SW8, pc.genes = epi_mac_fibr_SW8@var.genes, pcs.compute = 50)
jpeg(file = "PCElbowPlot_epi_mac_fibr_SW8.jpeg", width = 15, height = 15, units = "cm", res = 500)
PCElbowPlot(epi_mac_fibr_SW8, num.pc = 50)
dev.off()
epi_mac_fibr_SW8 <- FindClusters(epi_mac_fibr_SW8, reduction.type = "pca", dims.use = 1:22, resolution = 0.6, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
epi_mac_fibr_SW8 <- RunTSNE(epi_mac_fibr_SW8, dims.use = 1:22, do.fast = T)
jpeg(file = "tsne_epi_mac_fibr_SW8.jpeg", width = 15, height = 15, units = "cm", res = 500)
TSNEPlot(epi_mac_fibr_SW8, do.label = TRUE, pt.size = 0.8); dev.off()

# epi_mac_fibr_LW14 Reclustering:
epi_mac_fibr_LW14 <- FindVariableGenes(epi_mac_fibr_LW14, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
length(x = epi_mac_fibr_LW14@var.genes) #1290
epi_mac_fibr_LW14 <- NormalizeData(epi_mac_fibr_LW14, normalization.method = "LogNormalize", scale.factor = 10000)
epi_mac_fibr_LW14 <- ScaleData(epi_mac_fibr_LW14, vars.to.regress = c("nUMI", "percent.mito"))
epi_mac_fibr_LW14 <- RunPCA(epi_mac_fibr_LW14, pc.genes = epi_mac_fibr_LW14@var.genes, pcs.compute = 50)
jpeg(file = "epi_mac_fibr_LW14_epi_mac_fibr_SW8.jpeg", width = 15, height = 15, units = "cm", res = 500)
PCElbowPlot(epi_mac_fibr_LW14, num.pc = 50); dev.off()
epi_mac_fibr_LW14 <- FindClusters(epi_mac_fibr_LW14, reduction.type = "pca", dims.use = 1:25, resolution = 0.6, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
epi_mac_fibr_LW14 <- RunTSNE(epi_mac_fibr_LW14, dims.use = 1:25, do.fast = T)
jpeg(file = "tsne_epi_mac_fibr_LW14.jpeg", width = 15, height = 15, units = "cm", res = 500)
TSNEPlot(epi_mac_fibr_LW14, do.label = TRUE, pt.size = 0.8); dev.off()

save(epi_mac_fibr_LW14_subset, file ="epi_mac_fibr_LW14_subset.Robj")
load(file ="epi_mac_fibr_LW14.Robj")

epi_mac_fibr_LW14_subset = SubsetData(epi_mac_fibr_LW14, subset.raw = TRUE)
epi_mac_fibr_SW8_subset = SubsetData(epi_mac_fibr_SW8, subset.raw = TRUE)




PCElbowPlot(epi_recluster, num.pc = 50) ## 20 significant PCs
epi_recluster <- FindClusters(epi_recluster, reduction.type = "pca", dims.use = 1:20, resolution = 1.5, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
epi_recluster <- StashIdent(epi_recluster, save.name = "ClusterNames_epi_recluster")
epi_recluster <- RunTSNE(epi_recluster, dims.use = 1:20, do.fast = T)
TSNEPlot(epi_recluster, do.label = TRUE, pt.size = 0.8)



# Macrophases - 11, 2, 22
# Epithelial Cells - 15, 21, 6, 10, 1
# Fibroblasts - 3,0,13,5,12,20,9,7



###### Cell Phone Database Analysis on Synergy ######
load(file ="epi_mac_fibr_LW14_subset.Robj")
load(file ="epi_mac_fibr_SW8_subset.Robj")
TSNEPlot(epi_mac_fibr_SW8_subset, do.label = TRUE, pt.size = 0.8)
FeaturePlot(epi_mac_fibr_SW8_subset, features.plot = c("Dpt"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(epi_mac_fibr_SW8_subset, features.plot = c("Cd68"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(epi_mac_fibr_SW8_subset, features.plot = c("Krt10"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(epi_mac_fibr_SW8_subset, features.plot = c("Tnf"), reduction.use = "tsne", cols.use = c("yellow", "red"))

TSNEPlot(epi_mac_fibr_LW14_subset, do.label = TRUE, pt.size = 0.8)
FeaturePlot(epi_mac_fibr_LW14_subset, features.plot = c("Dpt"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(epi_mac_fibr_LW14_subset, features.plot = c("Crabp1"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(epi_mac_fibr_LW14_subset, features.plot = c("Adamts18"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(epi_mac_fibr_LW14_subset, features.plot = c("Pdgfra"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(epi_mac_fibr_LW14_subset, features.plot = c("Bmp3"), reduction.use = "tsne", cols.use = c("yellow", "red"))

FeaturePlot(epi_mac_fibr_LW14_subset, features.plot = c("Cd68"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(epi_mac_fibr_LW14_subset, features.plot = c("Tnf"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(epi_mac_fibr_LW14_subset, features.plot = c("Krt10"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(epi_mac_fibr_LW14_subset, features.plot = c("Wnt3"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(epi_mac_fibr_LW14_subset, features.plot = c("Cd3g"), reduction.use = "tsne", cols.use = c("yellow", "red"))

sample_1 <- grep("-1", colnames(LWD14fibr_subset@data), value = T); length(sample_1)
sample_2 <- grep("-2", colnames(LWD14fibr_subset@data), value = T); length(sample_2)
sample_3 <- grep("-3", colnames(LWD14fibr_subset@data), value = T); length(sample_3)
sample_4 <- grep("-4", colnames(LWD14fibr_subset@data), value = T); length(sample_4)

#jpeg(file = "tsne_LWD14fibr.jpeg", width = 15, height = 15, units = "cm", res = 500)
s_tSNE <- as.data.frame(GetCellEmbeddings(LWD14fibr_subset, reduction.type = "tsne"))
s_tSNE[sample_1, "Sample"] <- "LWC D14"
s_tSNE[sample_2, "Sample"] <- "LWC D14"
s_tSNE[sample_3, "Sample"] <- "LWO D14"
s_tSNE[sample_4, "Sample"] <- "LWO D14"
ggplot(s_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample),size = 0.3) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()



#### 4, 7, 11, 10, 14, 3, 1
LWD14fibr <- names(epi_mac_fibr_LW14_subset@ident[epi_mac_fibr_LW14_subset@ident %in% c(4, 7, 11, 10, 14, 3, 1)])
sub.LWD14fibr <- FilterCells(epi_mac_fibr_LW14_subset, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = LWD14fibr)

LWD14fibr <- FindVariableGenes(sub.LWD14fibr, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
length(x = LWD14fibr@var.genes) #1037
LWD14fibr <- NormalizeData(LWD14fibr, normalization.method = "LogNormalize", scale.factor = 10000)
LWD14fibr <- ScaleData(LWD14fibr, vars.to.regress = c("nUMI", "percent.mito"))
LWD14fibr <- RunPCA(LWD14fibr, pc.genes = LWD14fibr@var.genes, pcs.compute = 50)
jpeg(file = "LWD14fibr.jpeg", width = 15, height = 15, units = "cm", res = 500)
PCElbowPlot(LWD14fibr, num.pc = 50); dev.off()
LWD14fibr <- FindClusters(LWD14fibr, reduction.type = "pca", dims.use = 1:25, resolution = 0.6, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
LWD14fibr <- RunTSNE(LWD14fibr, dims.use = 1:25, do.fast = T)
jpeg(file = "tsne_LWD14fibr.jpeg", width = 15, height = 15, units = "cm", res = 500)
TSNEPlot(LWD14fibr, do.label = TRUE, pt.size = 0.8); dev.off()
save(LWD14fibr_subset, file ="LWD14fibr_subset.Robj")
LWD14fibr_subset = SubsetData(LWD14fibr, subset.raw = TRUE)

############ Sox18 +ve Fibros Charecterization ############
load(file ="LWD14fibr_subset.Robj")
TSNEPlot(LWD14fibr_subset, do.label = TRUE, pt.size = 0.8)
FeaturePlot(LWD14fibr_subset, features.plot = c("Crabp1"), reduction.use = "tsne", cols.use = c("yellow", "red"))

LWD14fibr_subset <- RunPHATE(object = LWD14fibr_subset)
TSNEPlot(LWD14fibr_subset, do.label = F, pt.size = 0.8)
DimPlot(LWD14fibr_subset, reduction.use = 'phate')
FeaturePlot(LWD14fibr_subset, features.plot = c("Fgfr"), reduction.use = "phate", cols.use = c("grey", "blue"))

cc.genes <- readLines(con = "/Users/SarthakSinha/Desktop/Seurat_Wounding/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")
library(stringr)
cc.genes <- str_to_title(cc.genes)
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
LWD14fibr_subset_CC <- CellCycleScoring(LWD14fibr_subset, s.genes = s.genes, 
                                          g2m.genes = g2m.genes, set.ident = TRUE)

TSNEPlot(LWD14fibr_subset_CC, do.label = F, pt.size = 0.8)
DimPlot(LWD14fibr_subset_CC, reduction.use = 'phate')


FeaturePlot(LWD14fibr_subset, features.plot = c("Sox18"), reduction.use = "tsne", cols.use = c("yellow", "red"))

s_tSNE <- as.data.frame(GetCellEmbeddings(LWD14fibr_subset, reduction.type = "tsne"))
s_tSNE[sample_1, "Sample"] <- "LWC D14"
s_tSNE[sample_2, "Sample"] <- "LWC D14"
s_tSNE[sample_3, "Sample"] <- "LWO D14"
s_tSNE[sample_4, "Sample"] <- "LWO D14"


############ Large Wound Center Fibros ONLY Charecterization ############
LWC_D14_fibr <- FilterCells(LWD14fibr_subset, cells.use = c(sample_1, sample_2), subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
LWC_D14_fibr <- FindVariableGenes(LWC_D14_fibr, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
length(x = LWC_D14_fibr@var.genes)
LWC_D14_fibr <- NormalizeData(LWC_D14_fibr, normalization.method = "LogNormalize", scale.factor = 10000)
LWC_D14_fibr <- ScaleData(LWC_D14_fibr, vars.to.regress = c("nUMI", "percent.mito"))
LWC_D14_fibr <- RunPCA(LWC_D14_fibr, pc.genes = LWC_D14_fibr@var.genes, pcs.compute = 50)
PCElbowPlot(LWC_D14_fibr, num.pc = 50) ## 20 significant PCs
LWC_D14_fibr <- FindClusters(LWC_D14_fibr, reduction.type = "pca", dims.use = 1:22, resolution = 0.6, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
LWC_D14_fibr <- RunTSNE(LWC_D14_fibr, dims.use = 1:22, do.fast = T)
TSNEPlot(LWC_D14_fibr, do.label = TRUE, pt.size = 0.8)
FeaturePlot(LWC_D14_fibr, features.plot = c("Rspo3"), reduction.use = "tsne", cols.use = c("yellow", "red"))
LWC_D14_fibr <- RunPHATE(object = LWC_D14_fibr)
TSNEPlot(LWC_D14_fibr, do.label = T, pt.size = 0.8)
DimPlot(LWC_D14_fibr, reduction.use = 'phate')
LWC_D14_fibr_CC <- CellCycleScoring(LWC_D14_fibr, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
TSNEPlot(LWC_D14_fibr, do.label = TRUE, pt.size = 0.8)
load(file = "LWC_D14_fibr.Robj")
FeaturePlot(LWC_D14_fibr, features.plot = c("Fgf20"), reduction.use = "tsne", cols.use = c("yellow", "red"))






save(LWC_D14_fibr, file ="LWC_D14_fibr.Robj")
save(LWC_D14_fibr_CC, file ="LWC_D14_fibr_CC.Robj")
FeaturePlot(LWC_D14_fibr, features.plot = c("Mki67", "Igfbp4", "Crabp1", "Rspo3"), reduction.use = "tsne", cols.use = c("yellow", "red"))
LWC_D14_fibr_subset = SubsetData(LWC_D14_fibr, subset.raw = TRUE)



LWC_D14_fibr_subset <- RunDiffusion(LWC_D14_fibr_subset,genes.use = LWC_D14_fibr_subset@var.genes)
DMPlot(LWC_D14_fibr_subset, do.label = TRUE)
FeaturePlot(LWC_D14_fibr, features.plot = c("Mki67", "Igfbp4", "Crabp1", "Rspo3"), reduction.use = "tsne", cols.use = c("yellow", "red"))

### 4, 2, 3, 7:
fibro_trajec <- names(LWC_D14_fibr_subset@ident[LWC_D14_fibr_subset@ident %in% c(4, 2, 3, 7)])
sub.fibro_trajec <- FilterCells(LWC_D14_fibr_subset, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = fibro_trajec)
LWC_D14_fibr_trajec <- FindVariableGenes(sub.fibro_trajec, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
length(x = LWC_D14_fibr_trajec@var.genes)
LWC_D14_fibr_trajec <- NormalizeData(LWC_D14_fibr_trajec, normalization.method = "LogNormalize", scale.factor = 10000)
LWC_D14_fibr_trajec <- ScaleData(LWC_D14_fibr_trajec, vars.to.regress = c("nUMI", "percent.mito"))
LWC_D14_fibr_trajec <- RunPCA(LWC_D14_fibr_trajec, pc.genes = LWC_D14_fibr_trajec@var.genes, pcs.compute = 50)
PCElbowPlot(LWC_D14_fibr_trajec, num.pc = 50) ## 20 significant PCs
LWC_D14_fibr_trajec <- FindClusters(LWC_D14_fibr_trajec, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
LWC_D14_fibr_trajec <- RunTSNE(LWC_D14_fibr_trajec, dims.use = 1:20, do.fast = T)
TSNEPlot(LWC_D14_fibr_trajec, do.label = TRUE, pt.size = 0.8)
LWC_D14_fibr_trajec <- RunDiffusion(LWC_D14_fibr_trajec,genes.use = LWC_D14_fibr_trajec@var.genes)
DMPlot(LWC_D14_fibr_trajec, do.label = TRUE)
load(file = "LWC_D14_fibr_trajec.Robj")

current.cluster.ids <- c(0, 1, 2, 3, 4)
new.cluster.ids <- c("Crabp1+ve Fibros", "Cd200+ve DS", "Mitotic Fibros", "S100a4+ve DSCs", 
                     "Rspo3+ve Condensate")
LWC_D14_fibr_trajec@ident <- plyr::mapvalues(x = LWC_D14_fibr_trajec@ident, from = current.cluster.ids, to = new.cluster.ids)
pdf("tsne_trial.pdf")
TSNEPlot(object = LWC_D14_fibr_trajec, do.label = TRUE, pt.size = 0.5)

####### SWNE on Synergy for LWC_D14_fibr_trajec #######
load(file = "LWC_D14_fibr_trajec_new.Robj")
var.genes <- LWC_D14_fibr_trajec_new@var.genes
length(var.genes)
cell.clusters <- LWC_D14_fibr_trajec_new@ident; names(cell.clusters) <- LWC_D14_fibr_trajec_new@cell.names;
levels(cell.clusters)
LWC_D14_fibr_trajec_new <- BuildSNN(LWC_D14_fibr_trajec_new, dims.use = 1:20, k.param = 20, prune.SNN = 1/20)
genes.embed <- c("Mki67", "Cd200", "Crabp1", "S100a4","Rspo3")
swne.embedding <- RunSWNE(LWC_D14_fibr_trajec_new, k = 10, var.genes = var.genes, genes.embed = genes.embed)
jpeg(file = "swne.embedding_LWC_D14_fibr_new", width = 25, height = 20, units = "cm", res = 500)
PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = cell.clusters,
         do.label = T, label.size = 3.5, pt.size = 1.5, show.legend = F,
         seed = 42)
dev.off()

TSNEPlot(object = LWC_D14_fibr_trajec_new, do.label = TRUE, pt.size = 0.5)
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Cd200", "S100a4", "Crabp1", "Rspo3"), reduction.use = "tsne", cols.use = c("yellow", "red"))

current.cluster.ids <- c(0, 1, 2, 3, 4)
new.cluster.ids <- c("Crabp1+ve Fibros", "Cd200+ve DS", "S100a4 DSC", "Crabp1+ve Fibros", 
                     "Rspo3+ve Condensate")
LWC_D14_fibr_trajec_new@ident <- plyr::mapvalues(x = LWC_D14_fibr_trajec_new@ident, from = current.cluster.ids, to = new.cluster.ids)
save(LWC_D14_fibr_trajec_new, file ="LWC_D14_fibr_trajec_new.Robj")
load("LWC_D14_fibr_trajec_new.Robj")

FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Mki67"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q50", min.cutoff = "q5")
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Pcna"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q50", min.cutoff = "q30")
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Cdkn1a", "Btg1"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q70", min.cutoff = "q10")

FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Fabp5"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q70", min.cutoff = "q10")
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Crabp1"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q70", min.cutoff = "q10")
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Ifitm3"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q70", min.cutoff = "q10")

VlnPlot(LWC_D14_fibr_trajec_new,features.plot = c("Mki67"))


############# Large Wound D14 Fibroblasts Scored for Cell Cycle: #############
LWC_D14_fibr_trajec_new_CC = CellCycleScoring(LWC_D14_fibr_trajec_new, s.genes = s.genes, 
                                              g2m.genes = g2m.genes, set.ident = TRUE)

jpeg(file = "LWC_D14_fibr_trajec_new_CC.jpeg", width = 12, height = 12, units = "cm", res = 500)
DMPlot(LWC_D14_fibr_trajec_new_CC, do.label = F, pt.size = 0.75)
dev.off()

TSNEPlot(LWC_D14_fibr_trajec_new_CC, do.label = F, pt.size = 0.8)

LWC_D14_fibr_trajec_new_CC_analysis = LWC_D14_fibr_trajec_new_CC@meta.data
write.csv(LWC_D14_fibr_trajec_new_CC_analysis,file = "LWC_D14_fibr_trajec_new_CC_analysis.csv")
############# Creating a stacked bar chart to display this:

#Create data
data_1 <- read.csv("LWC_D14_fibr_trajec_new_CC_analysis_processed.csv", stringsAsFactors = FALSE, header = F)
colnames(data_1)=c("Fibro","DSC","DS","Condensate")
rownames(data_1)=c("G1","S","G2/M")
data_1 =  data.matrix(data_1, rownames.force = NA)
View(data_1)

#create color palette:
library(RColorBrewer)
coul = c("#EE5D58", "#518AFE", "#28B022")
# Make a stacked barplot--> it will be in %!
jpeg(file = "data_1.jpeg", width = 13, height = 15, units = "cm", res = 500)
barplot(data_1, col=coul , border="black", space= c(1,1))
dev.off()

###################################### CELL CYCLE SCORING TO BE SAVED!! ######################################








#### Fibro Trajectory without Mitotic Fibros
fibro_trajec <- names(LWC_D14_fibr@ident[LWC_D14_fibr@ident %in% c(2, 3, 7)])
sub.fibro_trajec <- FilterCells(LWC_D14_fibr, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = fibro_trajec)
LWC_D14_fibr_trajec <- FindVariableGenes(sub.fibro_trajec, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
length(x = LWC_D14_fibr_trajec@var.genes)
LWC_D14_fibr_trajec <- NormalizeData(LWC_D14_fibr_trajec, normalization.method = "LogNormalize", scale.factor = 10000)
LWC_D14_fibr_trajec <- ScaleData(LWC_D14_fibr_trajec, vars.to.regress = c("nUMI", "percent.mito"))
LWC_D14_fibr_trajec <- RunPCA(LWC_D14_fibr_trajec, pc.genes = LWC_D14_fibr_trajec@var.genes, pcs.compute = 50)
PCElbowPlot(LWC_D14_fibr_trajec, num.pc = 50) ## 20 significant PCs
LWC_D14_fibr_trajec <- FindClusters(LWC_D14_fibr_trajec, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
LWC_D14_fibr_trajec <- RunTSNE(LWC_D14_fibr_trajec, dims.use = 1:20, do.fast = T)

load(file = "LWC_D14_fibr.Robj")
TSNEPlot(LWC_D14_fibr, do.label = TRUE, pt.size = 0.8)

### hello current
FeaturePlot(LWC_D14_fibr, features.plot = c("Crabp1", "Rspo3", "Sox18", "Top2a", "Mest", "Runx1"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(LWC_D14_fibr, features.plot = c("Fabp5"), reduction.use = "tsne", cols.use = c("yellow", "red"))


TSNEPlot(LWC_D14_fibr, do.label = TRUE, pt.size = 0.8)
current <- c(0,1,2,3,4,5,6,7,8,9)
new <- c("Lower Dermis LWC", #0
         "Upper Dermis LWC", #1
         "Upper Dermis LWC", #2
         "Upper Dermis LWC", #3
         "Mixed", #4
         "Lower Dermis LWC", #5
         "Lower Dermis LWC", #6
         "Upper Dermis LWC", #7
         "Mixed", #8
         "Upper Dermis LWC")


names(new) <- levels(LWC_D14_fibr)
LWC_D14_fibr <- RenameIdents(LWC_D14_fibr, new)
DimPlot(LWC_D14_fibr, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()

save(LWC_D14_fibr, file = "LWC_D14_fibr_defined.Robj")
load(file = "LWC_D14_fibr_defined.Robj")

LWC_D14_fibr@ident <- plyr::mapvalues(LWC_D14_fibr@ident, from = current, to = new)

jpeg(file = "LWC_D14_fibr_trial.jpeg", width = 15, height = 15, units = "cm", res = 500)
#TSNEPlot(LWC_D14_fibr, do.label = F, pt.size = 0.7, colors.use = c("#28B1B7", "#E95F5C", "black"))


TSNEPlot(LWC_D14_fibr)

dev.off()



##### Gene Expression Signatures Supplementary Table:
Lower_Dermal_LWC.markers <- FindMarkers(LWC_D14_fibr, ident.1 = c("Lower Dermis LWC"), min.pct = 0.20, test.use = "negbinom")
Lower_Dermal_LWC.markers <- Lower_Dermal_LWC.markers[order(Lower_Dermal_LWC.markers$avg_logFC),]
write.csv(Lower_Dermal_LWC.markers,file = "Lower_Dermal_LWC.markers.csv")
#
Upper_Dermal_LWC.markers <- FindMarkers(LWC_D14_fibr, ident.1 = c("Upper Dermis LWC"), min.pct = 0.20, test.use = "negbinom")
Upper_Dermal_LWC.markers <- Upper_Dermal_LWC.markers[order(Upper_Dermal_LWC.markers$avg_logFC),]
write.csv(Upper_Dermal_LWC.markers,file = "Upper_Dermal_LWC.markers.csv")
#
TSNEPlot(LWD14fibr_subset)
LWC.markers <- FindMarkers(LWD14fibr_subset, ident.1 = c("LWC"), ident.2 = c("LWP"), min.pct = 0.2, test.use = "negbinom")
LWC.markers <- LWC.markers[order(LWC.markers$avg_logFC),]
write.csv(LWC.markers,file = "LWC.markers_1.csv")
#
LWP.markers <- FindMarkers(LWD14fibr_subset, ident.1 = c("LWP"),ident.2 = c("LWC"), min.pct = 0.20, test.use = "negbinom")
LWP.markers <- LWP.markers[order(LWP.markers$avg_logFC),]
write.csv(LWP.markers,file = "LWP.markers_1.csv")
#
Neo_Condensate.markers <- FindMarkers(LWC_D14_fibr_trajec_new, ident.1 = c("Rspo3+ve Condensate"), min.pct = 0.10, test.use = "negbinom")
Neo_Condensate.markers <- Neo_Condensate.markers[order(Neo_Condensate.markers$avg_logFC),]
write.csv(Neo_Condensate.markers,file = "Neo_Condensate.markers.csv")
#
Neo_DSC.markers <- FindMarkers(LWC_D14_fibr_trajec_new, ident.1 = c("S100a4 DSCs"), min.pct = 0.20, test.use = "negbinom")
Neo_DSC.markers <- Neo_DSC.markers[order(Neo_DSC.markers$avg_logFC),]
write.csv(Neo_DSC.markers,file = "Neo_DSC.markers.csv")
#
Neo_CTS.markers <- FindMarkers(LWC_D14_fibr_trajec_new, ident.1 = c("S100a4 DSCs"), min.pct = 0.20, test.use = "negbinom")
Neo_CTS.markers <- Neo_CTS.markers[order(Neo_CTS.markers$avg_logFC),]
write.csv(Neo_CTS.markers,file = "Neo_CTS.markers.csv")


FeaturePlot(LWD14fibr_subset, features.plot = c("Acta2"), reduction.use = "tsne", cols.use = c("yellow", "red"))

##### 


load(file = "LWC_D14_fibr_defined.Robj")
library(Seurat)
VlnPlot(LWC_D14_fibr, features.plot = c("Crabp1"), do.sort = FALSE, ident.include = c("Lower Dermis LWC", "Upper Dermis LWC"), cols.use = c("#28B1B7", "#E95F5C"))
VlnPlot(LWC_D14_fibr, features.plot = c("Ly6a"), do.sort = FALSE, ident.include = c("Lower Dermis LWC", "Upper Dermis LWC"), cols.use = c("#28B1B7", "#E95F5C"))
VlnPlot(LWC_D14_fibr, features.plot = c("Fabp5"), do.sort = FALSE, ident.include = c("Lower Dermis LWC", "Upper Dermis LWC"), cols.use = c("#28B1B7", "#E95F5C"))

VlnPlot(LWC_D14_fibr, features.plot = c("Prss35"), do.sort = FALSE, ident.include = c("Lower Dermis LWC", "Upper Dermis LWC"), cols.use = c("#28B1B7", "#E95F5C"))

VlnPlot(LWC_D14_fibr, features.plot = c("Prss35"), do.sort = FALSE, ident.include = c("Lower Dermis LWC", "Upper Dermis LWC"), cols.use = c("#28B1B7", "#E95F5C"))



##### RNA/ATAC Validation:
VlnPlot(LWC_D14_fibr, features = c("Lum"), do.sort = FALSE, ident.include = c("Lower Dermis LWC", "Upper Dermis LWC"), cols.use = c("#28B1B7", "#E95F5C"))
VlnPlot(LWC_D14_fibr, features = c("Grem2"), do.sort = FALSE, ident.include = c("Lower Dermis LWC", "Upper Dermis LWC"), cols.use = c("#28B1B7", "#E95F5C"))

VlnPlot(LWC_D14_fibr, features = c("Murc"), do.sort = FALSE, ident.include = c("Lower Dermis LWC", "Upper Dermis LWC"), cols.use = c("#28B1B7", "#E95F5C"))
VlnPlot(LWC_D14_fibr, features = c("Prr16"), do.sort = FALSE, ident.include = c("Lower Dermis LWC", "Upper Dermis LWC"), cols.use = c("#28B1B7", "#E95F5C"))
VlnPlot(LWC_D14_fibr, features = c("Col3a1"), do.sort = FALSE, ident.include = c("Lower Dermis LWC", "Upper Dermis LWC"), cols.use = c("#28B1B7", "#E95F5C"))

VlnPlot(LWC_D14_fibr, features = c("Mest"), do.sort = FALSE, ident.include = c("Lower Dermis LWC", "Upper Dermis LWC"), cols.use = c("#28B1B7", "#E95F5C"))


VlnPlot(LWC_D14_fibr, features = c("Prss35"), do.sort = FALSE, ident = c("Lower Dermis LWC", "Upper Dermis LWC"), cols = c("#28B1B7", "#E95F5C"))


FeaturePlot(LWC_D14_fibr, features = c("Grem2"))

FeaturePlot(LWC_D14_fibr, features = c("Prss35"))
FeaturePlot(LWC_D14_fibr, features = c("Prss35"))


VlnPlot(LWC_D14_fibr, features = c("Prss35"), do.sort = FALSE, ident.include = c("Lower Dermis LWC", "Upper Dermis LWC"), cols.use = c("#28B1B7", "#E95F5C"))

features <- c("Crabp1", "Fabp5", "Prss35")
DotPlot(LWC_D14_fibr, features = features) + RotatedAxis()



jpeg(file = "Dlk1Crabp1_LW.jpeg", width = 15, height = 15, units = "cm", res = 500)
FeaturePlot(LWC_D14_fibr,
            features.plot = c("Crabp1", "Ly6a"),
            cols.use = c("grey", "blue", "red"), 
            reduction.use = "tsne",
            overlay = TRUE)
dev.off()





gene.set <- c("Ly6a", "Dlk1", "Mest")
mean.exp <- colMeans(x = LWC_D14_fibr@scale.data[gene.set, ], na.rm = TRUE)
if (all(names(x = mean.exp) == rownames(x = LWC_D14_fibr@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'LWC_D14_fibr@meta.data':\n", 
      "adding gene set mean expression values in 'LWC_D14_fibr@meta.data$gene.set.score'")
  LWC_D14_fibr@meta.data$gene.set.score <- mean.exp}

gene.set.2 <- c("Crabp1", "Fabp5", "Runx1")
mean.exp.2 <- colMeans(x = LWC_D14_fibr@scale.data[gene.set.2, ], na.rm = TRUE)
if (all(names(x = mean.exp.2) == rownames(x = LWC_D14_fibr@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'LWC_D14_fibr@meta.data':\n", 
      "adding gene set mean expression values in 'LWC_D14_fibr@meta.data$gene.set.score'")
  LWC_D14_fibr@meta.data$gene.set.score.2 <- mean.exp.2}


jpeg(file = "gene_sets_lwc_upper_lower.jpeg", width = 15, height = 15, units = "cm", res = 500)
FeaturePlot(LWC_D14_fibr, features.plot = c("gene.set.score", "gene.set.score.2"),
            cols.use = c("grey", "blue", "red"), 
            reduction.use = "tsne",
            overlay = TRUE)
dev.off()



###### Colocalization of Wnt and Sonic Hedgehog
gene.set <- c("Axin2", "Nkd1")
mean.exp <- colMeans(x = LWC_D14_fibr@scale.data[gene.set, ], na.rm = TRUE)
if (all(names(x = mean.exp) == rownames(x = LWC_D14_fibr@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'LWC_D14_fibr@meta.data':\n", 
      "adding gene set mean expression values in 'LWC_D14_fibr@meta.data$gene.set.score'")
  LWC_D14_fibr@meta.data$gene.set.score <- mean.exp}


gene.set.2 <- c("Gli1", "Ptch1")
mean.exp.2 <- colMeans(x = LWC_D14_fibr@scale.data[gene.set.2, ], na.rm = TRUE)
if (all(names(x = mean.exp.2) == rownames(x = LWC_D14_fibr@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'LWC_D14_fibr@meta.data':\n", 
      "adding gene set mean expression values in 'LWC_D14_fibr@meta.data$gene.set.score'")
  LWC_D14_fibr@meta.data$gene.set.score.2 <- mean.exp.2}

jpeg(file = "gene_sets_lwc_upper_lower.jpeg", width = 15, height = 15, units = "cm", res = 500)
FeaturePlot(LWC_D14_fibr, features.plot = c("Axin2", "Gli1"),
            cols.use = c("grey", "blue", "red"), 
            reduction.use = "tsne",
            overlay = TRUE)
dev.off()



load(file = "LWC_D14_fibr.Robj")

jpeg(file = "Axin2_LWC_PvsC.jpeg", width = 15, height = 15, units = "cm", res = 500)
FeaturePlot(LWC_D14_fibr, features.plot = c("Axin2"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q15", max.cutoff = "q70")
dev.off()

LWC_D14_fibr = UpdateSeuratObject(LWC_D14_fibr)
DimPlot(LWC_D14_fibr, features = c("Igfbp2"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q15", max.cutoff = "q70")

jpeg(file = "Igfbp2_LWC_D14_fibr.jpeg", width = 15, height = 15, units = "cm", res = 500)
FeaturePlot(LWC_D14_fibr, features = c("Igfbp2"),
            cols = c("yellow", "red"), max.cutoff = "q70")
dev.off()



FeaturePlot(LWC_D14_fibr, features = c("Grem2"),
            cols = c("yellow", "red"), max.cutoff = "q90")


FeaturePlot(LWC_D14_fibr, features = c("Dcn"),
            cols = c("yellow", "red"), max.cutoff = "q90")

FeaturePlot(LWC_D14_fibr, features = c("Casz1"),
            cols = c("yellow", "red"), max.cutoff = "q90")


VlnPlot(LWC_D14_fibr, features = c("Casz1"))

View(LWC_D14_fibr@meta.data)


FeaturePlot(LWD14fibr_subset, features.plot = c("Axin2"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10", max.cutoff = "q70")


FeaturePlot(LWC_D14_fibr, features.plot = c("Nkd1"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10", max.cutoff = "q70")
FeaturePlot(LWD14fibr_subset, features.plot = c("Nkd1"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10", max.cutoff = "q70")



FeaturePlot(LWC_D14_fibr, features.plot = c("Gli1"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10", max.cutoff = "q70")
FeaturePlot(LWD14fibr_subset, features.plot = c("Gli1"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10", max.cutoff = "q70")

FeaturePlot(LWC_D14_fibr, features.plot = c("Ptch1"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10", max.cutoff = "q70")
FeaturePlot(LWD14fibr_subset, features.plot = c("Ptch1"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10", max.cutoff = "q70")





FeaturePlot(LWC_D14_fibr, features.plot = c("Lef1"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10", max.cutoff = "q70")


FeaturePlot(LWC_D14_fibr, features.plot = c("Dkk1"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10", max.cutoff = "q70")





load("LWC_D14_fib.Robj")

jpeg(file = "Dlk1_LWC_D14_fibr_VlnPlot_lower_upper.jpeg", width = 15, height = 15, units = "cm", res = 500)
VlnPlot(LWC_D14_fibr, features.plot = c("Axin2"), do.sort = F, ident.include = c("Lower Dermis LWC", "Upper Dermis LWC"), cols.use = c("#28B1B7", "#E95F5C"))
dev.off()




LWC_D14_fibr_CC <- CellCycleScoring(LWC_D14_fibr, s.genes = s.genes, 
                                      g2m.genes = g2m.genes, set.ident = TRUE)
TSNEPlot(LWC_D14_fibr_CC, do.label = F, pt.size = 0.7, colors.use = c("#28B1B7", "#E95F5C", "black"))



gene.set <- c("Ly6a", "Dlk1", "Mest")
mean.exp <- colMeans(x = LWD14fibr_subset@scale.data[gene.set, ], na.rm = TRUE)
if (all(names(x = mean.exp) == rownames(x = LWD14fibr_subset@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'LWD14fibr_subset@meta.data':\n", 
      "adding gene set mean expression values in 'LWD14fibr_subset@meta.data$gene.set.score'")
  LWD14fibr_subset@meta.data$gene.set.score <- mean.exp}

gene.set.2 <- c("Crabp1", "Fabp5", "Runx1")
mean.exp.2 <- colMeans(x = LWD14fibr_subset@scale.data[gene.set.2, ], na.rm = TRUE)
if (all(names(x = mean.exp.2) == rownames(x = LWD14fibr_subset@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'LWD14fibr_subset@meta.data':\n", 
      "adding gene set mean expression values in 'LWD14fibr_subset@meta.data$gene.set.score'")
  LWD14fibr_subset@meta.data$gene.set.score.2 <- mean.exp.2}

jpeg(file = "LWCvsP_gene_sets.jpeg", width = 15, height = 15, units = "cm", res = 500)
FeaturePlot(LWD14fibr_subset,
            features.plot = c("gene.set.score", "gene.set.score.2"),
            cols.use = c("grey", "blue", "red"), 
            reduction.use = "tsne",
            overlay = TRUE)
dev.off()



jpeg(file = "Crabp1_LWCvsLWP.jpeg", width = 15, height = 15, units = "cm", res = 500)
FeaturePlot(LWD14fibr_subset, features.plot = c("Runx1"), reduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q20")
dev.off()



### 
FeaturePlot(LWC_D14_fibr, features.plot = c("Fabp5"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10", max.cutoff = "q90")

### SCENIC DATA Extraction
x = SubsetData(LWC_D14_fibr, subset.raw = TRUE)
save(x, file = "x.Robj")









LWC_D14_fibr <- BuildClusterTree(
  LWC_D14_fibr,
  pcs.use = 1:25,
  do.reorder = F,
  reorder.numeric = F,
  do.plot=F)

jpeg(file = "PlotClusterTree_Wound_Samples.jpeg", width = 15, height = 15, units = "cm", res = 500)
PlotClusterTree(LWC_D14_fibr)
dev.off()



LWC_D14_fibr_trajec_new = LWC_D14_fibr_trajec
LWC_D14_fibr_trajec_new = SubsetData(LWC_D14_fibr_trajec_new, subset.raw = TRUE)
save(LWC_D14_fibr_trajec_new, file = "LWC_D14_fibr_trajec_new.Robj")

LWC_D14_fibr_trajec_new <- RunDiffusion(LWC_D14_fibr_trajec_new,genes.use = LWC_D14_fibr_trajec@var.genes)
LWC_D14_fibr_trajec_new <- RunDiffusion(LWC_D14_fibr_trajec_new,genes.use = LWC_D14_fibr_trajec@var.genes, dims.use = 1:10)
LWC_D14_fibr_trajec_new <- RunDiffusion(LWC_D14_fibr_trajec_new,genes.use = LWC_D14_fibr_trajec@var.genes, dims.use = 1:10, q.use = 0.01, max.dim = 5)
TSNEPlot(LWC_D14_fibr_trajec_new, do.label = TRUE, pt.size = 0.8)

DMPlot(LWC_D14_fibr_trajec_new, do.label = TRUE)
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Crabp1", "Rspo3", "Sox18"), reduction.use = "dm", cols.use = c("yellow", "red"))

LWC_D14_fibr_trajec_new_CC <- CellCycleScoring(LWC_D14_fibr_trajec_new, s.genes = s.genes, 
                                        g2m.genes = g2m.genes, set.ident = TRUE)



jpeg(file = "LWC_D14_fibr_trajec_new_Fplot_WS_markers.jpeg", width = 20, height = 30, units = "cm", res = 500)
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Itm2a",
                                                       "Timp1", 
                                                       "S100a4", 
                                                       "Mif", 
                                                       "Pkm"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q10")
dev.off()



########## scATAC-Seq Validation ########## 
# Sema3a - chr5_13229634_13229956
# Btbd10 promoter - chr7_113377514_113378736
# Igfbp2 - chr1_72823921_72825039 - repressed
# Atp5g2 promoter - chr15_102692465_102694276
# Llgl1 - intermediady
# Sulf2 - Igfbp2-like NICE!
# Tmem119 - so so
# Afap1 - Igfbp2-like NICE!


load(file = "LWC_D14_fibr_trajec_new.Robj")
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Sema3a"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q10")
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Rapgef1"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q10")
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Btbd10"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q10")
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Abcc9"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q10")
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Lrrfip2"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q10")
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Cux1"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q10")

FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Igfbp2"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q10")
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Atp5g2"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q10")
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Phex"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q10")
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Dmrt2"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q10")

jpeg(file = "Enpp2_inductive_trajectory.jpeg", width = 15, height = 18, units = "cm", res = 500)
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Enpp2"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q10")
dev.off()



jpeg(file = "LWC_D14_fibr_trajec_new_Fplot_WS_markers_DP_r.jpeg", width = 20, height = 30, units = "cm", res = 500)
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Rspo3",
                                                       "Sox18", 
                                                       "Crabp1", 
                                                       "Cd200", 
                                                       "Stmn2",
                                                       "S100a4"), reduction.use = "dm", cols.use = c("yellow", "red"))
dev.off()

load(file = "LWC_D14_fibr_trajec_new.Robj")

#### PAGODA Export
j_j = SubsetData(LWC_D14_fibr_trajec_new, subset.raw = TRUE)
q_df = as.data.frame.array(j_j@data)
q_mx <- as.matrix(q_df)
save(q_mx, file = "q_mx.Robj")



FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Crabp1", "Rspo3", "Sox18", "Cd200", "S100a4", "Stmn2"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q10")

FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Eln"), reduction.use = "dm", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q10")


################# SCENIC Analysis on LWC_D14_fibr_trajec_new: #################
x = SubsetData(LWC_D14_fibr_trajec_new, subset.raw = TRUE)
TSNEPlot(x, do.label = TRUE, pt.size = 0.8)
save(x, file = "x.Robj")



































####### Merging Seurat Objects to Overlay Plcode and Condensate #######
load(file = "placode_transition_3.Robj")
TSNEPlot(placode_transition_3, do.label = TRUE, pt.size = 0.8)

placode_2step <- names(placode_transition_3@ident[placode_transition_3@ident %in% c("Placode Formation", "Placode Expansion")])
placode_2step <- FilterCells(placode_transition_3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = placode_2step)

TSNEPlot(placode_2step, do.label = TRUE, pt.size = 0.8)

LW14_plac_condensate = MergeSeurat(object1 = placode_2step, object2 =LWC_D14_fibr_trajec_new, add.cell.id1 = "Epi", add.cell.id2 = "Fibro")
LW14_plac_condensate <- FindVariableGenes(LW14_plac_condensate, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
length(x = LW14_plac_condensate@var.genes)
LW14_plac_condensate <- NormalizeData(LW14_plac_condensate, normalization.method = "LogNormalize", scale.factor = 10000)
LW14_plac_condensate <- ScaleData(LW14_plac_condensate, vars.to.regress = NULL)
LW14_plac_condensate <- RunPCA(LW14_plac_condensate, pc.genes = LW14_plac_condensate@var.genes, pcs.compute = 50)
PCElbowPlot(LW14_plac_condensate, num.pc = 50) ## 15 significant PCs
LW14_plac_condensate <- FindClusters(LW14_plac_condensate, reduction.type = "pca", dims.use = 1:20, resolution = 1, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
LW14_plac_condensate <- RunTSNE(LW14_plac_condensate, dims.use = 1:20, do.fast = T)
TSNEPlot(LW14_plac_condensate, do.label = TRUE, pt.size = 0.8)

FeaturePlot(LW14_plac_condensate, features.plot = c("Cd200",
                                                    "Ednrb",
                                                    "S100a4"), reduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q10")

current.cluster.ids <- c(0, 1, 2, 3, 4,5,6,7,8,9,10)
new.cluster.ids <- c("Crabp1+ Fibros", 
                     "Placode Formation", 
                     "Placode Formation", 
                     "Placode Formation", 
                     "Mitotic Placode",
                     "Cd200+ DS",
                     "Rspo3+ Condensate",
                     "Placode Formation",
                     "Crabp1+ Fibros",
                     "S100a4+ DSC",
                     "Placode Formation")
LW14_plac_condensate@ident <- plyr::mapvalues(x = LW14_plac_condensate@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(LW14_plac_condensate, do.label = TRUE, pt.size = 0.8)

LW14_plac_condensate = SubsetData(LW14_plac_condensate, subset.raw = TRUE)
LW14_plac_condensate <- StashIdent(LW14_plac_condensate, save.name = "Cell_names")

save(LW14_plac_condensate, file = "LW14_plac_condensate.Robj")
load(file = "LW14_plac_condensate.Robj")
CPD = as.data.frame.array(placode_2step@data)
CPD <- as.matrix(CPD)
write.csv(CPD, file = "CPD.csv")

CPD_2 = as.data.frame.array(LW14_plac_condensate@meta.data)
CPD_2 <- as.matrix(CPD_2)
write.csv(CPD_2, file = "CPD_2.csv")

install.packages("gProfileR")
library(gProfileR)

gene_list = rownames(CPD)
CPD_ortho_list = gorth(gene_list, source_organism = "mmusculus", target_organism = "hsapiens",
      region_query = F, numeric_ns = "", mthreshold = 1, filter_na = F,
      df = T)
write.csv(CPD_ortho_list, file = "CPD_ortho_list.csv")



# EXTRA
LWC_D14_fibr_trajec1.markers <- FindMarkers(LWC_D14_fibr_trajec, ident.1 = c(1), min.pct = 0.25, test.use = "negbinom")
LWC_D14_fibr_trajec1.markers <- LWC_D14_fibr_trajec1.markers[order(LWC_D14_fibr_trajec1.markers$avg_logFC),]
write.csv(LWC_D14_fibr_trajec1.markers,file = "LWC_D14_fibr_trajec1.markers.csv")

LWC_D14_fibr_trajec3.markers <- FindMarkers(LWC_D14_fibr_trajec, ident.1 = c(3), min.pct = 0.25, test.use = "negbinom")
LWC_D14_fibr_trajec3.markers <- LWC_D14_fibr_trajec3.markers[order(LWC_D14_fibr_trajec3.markers$avg_logFC),]
write.csv(LWC_D14_fibr_trajec3.markers,file = "LWC_D14_fibr_trajec3.markers.csv")


FeaturePlot(LWC_D14_fibr_trajec, features.plot = c("S100a4", "Timp1", "Cd200", "Rspo3", "Itm2a"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10", max.cutoff = "q90")

FeaturePlot(LWC_D14_fibr_trajec, features.plot = c("Mki67", "Igfbp4", "Crabp1", "Sox18"), reduction.use = "tsne", cols.use = c("yellow", "red"))
save(LWC_D14_fibr_trajec, file ="LWC_D14_fibr_trajec.Robj")
LWC_D14_fibr_trajec <- RunPHATE(object = LWC_D14_fibr_trajec)

LWC_D14_fibr_trajec <- RunDiffusion(LWC_D14_fibr_trajec,genes.use = LWC_D14_fibr_trajec@var.genes)
DMPlot(LWC_D14_fibr_trajec, do.label = TRUE)
load(file = "LWC_D14_fibr_trajec.Robj")



####################################### Large Wound Center vs Pheriphery #################################################################
sample_1 <- grep("-1", colnames(LWD14fibr_subset@data), value = T); length(sample_1)
sample_2 <- grep("-2", colnames(LWD14fibr_subset@data), value = T); length(sample_2)
sample_3 <- grep("-3", colnames(LWD14fibr_subset@data), value = T); length(sample_3)
sample_4 <- grep("-4", colnames(LWD14fibr_subset@data), value = T); length(sample_4)

s_tSNE <- as.data.frame(GetCellEmbeddings(LWD14fibr_subset, reduction.type = "tsne"))
s_tSNE[sample_1, "Sample"] <- "LWC D14"
s_tSNE[sample_2, "Sample"] <- "LWC D14"
s_tSNE[sample_3, "Sample"] <- "LWO D14"
s_tSNE[sample_4, "Sample"] <- "LWO D14"




################# SCENIC Data Extration - for sample colour overlaid #################
#######1 
LWD14fibr_subset_trial = LWD14fibr_subset
levels(LWD14fibr_subset_trial@meta.data$orig.ident) <- c("LWC D14", "LWO D14")
x = rownames(LWD14fibr_subset_trial@meta.data)
y = sample_4
for(i in 1:6235)
if (x[i] %in% y){LWD14fibr_subset_trial@meta.data$orig.ident[i] <- "LWO D14"}
View(LWD14fibr_subset_trial@meta.data)

LWD14fibr_subset_trial <- SetAllIdent(LWD14fibr_subset_trial, id = "orig.ident")
TSNEPlot(LWD14fibr_subset_trial)
save(LWD14fibr_subset_trial, file = "LWD14fibr_subset_trial.Robj")

x = SubsetData(LWD14fibr_subset_trial, subset.raw = TRUE)
TSNEPlot(x, do.label = TRUE, pt.size = 0.8)
save(x, file = "x.Robj")

#######2
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

jpeg(file = "D14LWC_LWP_SW8_UI_sample_of_interest_1_LWCvsP_seperated.jpeg", width = 15, height = 15, units = "cm", res = 500)
s_tSNE <- as.data.frame(GetCellEmbeddings(D14LWC_LWP_SW8_UI, reduction.type = "tsne"))
s_tSNE[sample_1, "Sample"] <- "LWC D14"
s_tSNE[sample_2, "Sample"] <- "LWC D14"
s_tSNE[sample_3, "Sample"] <- "LWP D14"
s_tSNE[sample_4, "Sample"] <- "LWP D14"
s_tSNE[sample_5, "Sample"] <- "SW D14"
s_tSNE[sample_6, "Sample"] <- "SW D14"
s_tSNE[Day18_Day8_Lib3_try1, "Sample"] <- "SW D8"
s_tSNE[sample_7, "Sample"] <- "UI"
ggplot(s_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample),size = 0.3) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()


D14LWC_LWP_SW8_UI_trial = D14LWC_LWP_SW8_UI
levels(D14LWC_LWP_SW8_UI_trial@meta.data$orig.ident) <- c("LWC D14", "LWP D14", "SW D14", "SW D8", "UI")
x = rownames(D14LWC_LWP_SW8_UI_trial@meta.data)

y = Day18_Day8_Lib3_try1
for(i in 1:13678)
  if (x[i] %in% y){D14LWC_LWP_SW8_UI_trial@meta.data$orig.ident[i] <- "SW D8"}

D14LWC_LWP_SW8_UI_trial <- SetAllIdent(D14LWC_LWP_SW8_UI_trial, id = "orig.ident")
TSNEPlot(D14LWC_LWP_SW8_UI_trial)

x = SubsetData(D14LWC_LWP_SW8_UI_trial, subset.raw = TRUE)
TSNEPlot(x, do.label = TRUE, pt.size = 0.8)
save(x, file = "x.Robj")


load(file = "D14LWC_LWP_SW8_UI.Robj")
D14LWC_LWP_SW8_UI_v3 = UpdateSeuratObject(D14LWC_LWP_SW8_UI)
VlnPlot(D14LWC_LWP_SW8_UI_v3, features = c("T"))


load(file = "LWD14fibr_subset.Robj")
LWD14fibr_subset_v3 = UpdateSeuratObject(LWD14fibr_subset)

jpeg(file = "LWD14fibr_subset_v3_Tagln.jpeg", width = 15, height = 15, units = "cm", res = 500)
VlnPlot(LWD14fibr_subset_v3, features = c("Tagln"))
dev.off()

pdf("Putative_DC_LWD14_fibros_tdT+_-_reclustered.pdf")


dev.off()

################# SCENIC Data Extration - for sample colour overlaid #################


jpeg(file = "Large_Wound_Center_vs_Pheriphery_tsne_categorized.jpeg", width = 15, height = 15, units = "cm", res = 500)
TSNEPlot(LWD14fibr_subset, do.label = TRUE, pt.size = 0.8)
dev.off()

jpeg(file = "Large_Wound_Center_vs_Pheriphery_tsne_categorized.jpeg", width = 15, height = 15, units = "cm", res = 500)
TSNEPlot(LWD14fibr_subset, do.label = TRUE, pt.size = 0.8)
dev.off()

LWP.markers <- FindMarkers(LWD14fibr_subset, ident.1 = c("LWP"), min.pct = 0.25, test.use = "negbinom")
LWP.markers <- LWP.markers[order(LWP.markers$avg_logFC),]
write.csv(LWP.markers,file = "LWP.markers.csv")

FeaturePlot(LWD14fibr_subset, features.plot = c("Rgs5"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(LWD14fibr_subset, features.plot = c("Cilp", "Crabp1"), reduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q15")

jpeg(file = "Dlk1Crabp1_LW.jpeg", width = 15, height = 15, units = "cm", res = 500)
FeaturePlot(LWD14fibr_subset,
            features.plot = c("Dlk1", "Crabp1"),
            cols.use = c("grey", "blue", "red"), 
            reduction.use = "tsne",
            overlay = TRUE)
dev.off()


jpeg(file = "Fabp5_LWCvsLWP.jpeg", width = 15, height = 15, units = "cm", res = 500)
VlnPlot(LWD14fibr_subset, features.plot = c("Fabp5"), do.sort = FALSE, ident.include = c("LWP", "LWC"), cols.use = c("#28B1B7", "#E95F5C"))
dev.off()

VlnPlot(LWD14fibr_subset, features.plot = c("Dlk1"))
VlnPlot(LWD14fibr_subset, features.plot = c("Ly6a"))
VlnPlot(LWD14fibr_subset, features.plot = c("Mest"))
VlnPlot(LWD14fibr_subset, features.plot = c("Blimp1"))



jpeg(file = "Hes1.jpeg", width = 15, height = 15, units = "cm", res = 500)
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Hes1"), reduction.use = "dm", cols.use = c("yellow", "red", "white"))
dev.off()



# EXTRA
pdf("Putative_DC_LWD14_fibros_tdT+_-_reclustered.pdf")
ggplot(s_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample),size = 0.3) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
TSNEPlot(LWD14fibr_subset, do.label = TRUE, pt.size = 0.8)
FeaturePlot(LWD14fibr_subset, features.plot = c("Crabp1"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(LWD14fibr_subset, features.plot = c("S100a4"), reduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q15")
FeaturePlot(LWD14fibr_subset, features.plot = c("Rspo3"), reduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q15")
FeaturePlot(LWD14fibr_subset, features.plot = c("Bmp4"), reTSNE_D14LWC_LWP_SW8_UI_sample_of_interest_1_for_refduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q15")
FeaturePlot(LWD14fibr_subset, features.plot = c("Frzb"), reduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q15")
FeaturePlot(LWD14fibr_subset, features.plot = c("Cd24a"), reduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q15")
FeaturePlot(LWD14fibr_subset, features.plot = c("Enpp2"), reduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q15")
FeaturePlot(LWD14fibr_subset, features.plot = c("Wif1"), reduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q15")
FeaturePlot(LWD14fibr_subset, features.plot = c("Prlr"), reduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q15")
FeaturePlot(LWD14fibr_subset, features.plot = c("Inhba"), reduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q15")
FeaturePlot(LWD14fibr_subset, features.plot = c("Fgf10"), reduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q15")
dev.off()
load(file ="LWD14fibr_subset.Robj")

TSNEPlot(LWD14fibr_subset, do.label = TRUE, pt.size = 0.8)

# Inserting Sample Identity for SCENIC:
View(LWD14fibr_subset@meta.data)
levels(LWD14fibr_subset@meta.data$orig.ident) <- c("LWC", "LWP")
microglia@meta.data[1:4568,]$orig.ident <- "LWC"
microglia@meta.data[4569:7167,]$orig.ident <- "LWP"







current.cluster.ids <- c(0, 1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15)
new.cluster.ids <- c("LWP",
                     "Mixed", 
                     "LWC",
                     "Mixed",
                     "LWC",
                     "Mixed",
                     "LWP",
                     "LWC",
                     "Mixed",
                     "LWC",
                     "LWC",
                     "Mixed",
                     "LWP",
                     "Mixed",
                     "Mixed",
                     "Mixed")
LWD14fibr_subset@ident <- plyr::mapvalues(x = LWD14fibr_subset@ident, from = current.cluster.ids, to = new.cluster.ids)

TSNEPlot(LWD14fibr_subset, do.label = TRUE, pt.size = 0.8)

LWD14fibr_subset.markers <- FindAllMarkers(LWD14fibr_subset, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
top10 <- LWD14fibr_subset.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

DoHeatmap(LWD14fibr_subset, genes.use = top10$gene, slim.col.label = T, remove.key = TRUE, disp.min = -2.0,
          disp.max = 3, col.low = "#5B7FAC",
          col.mid = "#FFF3BA",
          col.high = "#C6514E", 
          group.by = "ident")


VlnPlot(LWD14fibr_subset, features.plot = c("Runx1"))

FeaturePlot(LWD14fibr_subset, features.plot = c("Runx1"), reduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q15")


jpeg(file = "Crabp1_LWCvsLWP.jpeg", width = 15, height = 15, units = "cm", res = 500)
FeaturePlot(LWD14fibr_subset, features.plot = c("Runx1"), reduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q20")
dev.off()


jpeg(file = "Runx1_LWCvsLWP.jpeg", width = 15, height = 15, units = "cm", res = 500)
VlnPlot(LWD14fibr_subset, features.plot = c("Runx1"), do.sort = T, ident.include = c("LWC", "LWP"))
dev.off()


FeatureHeatmap(LWD14fibr_subset, features.plot = c("Runx1"), group.by = "stim", pt.size = 0.25, key.position = "top", max.exp = 3)
x = SubsetData(LWD14fibr_subset, subset.raw = TRUE)
TSNEPlot(x, do.label = TRUE, pt.size = 0.8)
save(x, file = "x.Robj")

save(LWD14fibr_subset, file = "LWD14fibr_subset.Robj")


#################### Cell cycle scoring of Large Wound Pheriphery and Center ####################
LWD14fibr_subset_CC <- CellCycleScoring(LWD14fibr_subset, s.genes = s.genes, 
                                           g2m.genes = g2m.genes, set.ident = TRUE)

TSNEPlot(LWD14fibr_subset_CC, do.label = F, pt.size = 0.3)
all_fibros_4_scenic_CC_analysis = all_fibros_4_scenic_CC@meta.data
write.csv(all_fibros_4_scenic_CC_analysis,file = "all_fibros_4_scenic_CC_analysis.csv")









LWC_LWP_exclusive_cells <- names(LWD14fibr_subset@ident[LWD14fibr_subset@ident %in% c("LWP", "LWC")])
LWC_LWP_exclusive_cells <- FilterCells(LWD14fibr_subset, subset.names = "nGene", low.thresholds = 20, high.thresholds = 200000, cells.use = LWC_LWP_exclusive_cells)
LWC_LWP_exclusive_cells <- FindVariableGenes(LWC_LWP_exclusive_cells, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
length(x = LWC_LWP_exclusive_cells@var.genes)
LWC_LWP_exclusive_cells <- NormalizeData(LWC_LWP_exclusive_cells, normalization.method = "LogNormalize", scale.factor = 10000)
LWC_LWP_exclusive_cells <- ScaleData(LWC_LWP_exclusive_cells, vars.to.regress = c("nUMI", "percent.mito"))
LWC_LWP_exclusive_cells <- RunPCA(LWC_LWP_exclusive_cells, pc.genes = LWC_LWP_exclusive_cells@var.genes, pcs.compute = 50)
PCElbowPlot(LWC_LWP_exclusive_cells, num.pc = 50)
LWC_LWP_exclusive_cells <- FindClusters(LWC_LWP_exclusive_cells, reduction.type = "pca", dims.use = 1:25, resolution = 0.6, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
LWC_LWP_exclusive_cells <- RunTSNE(LWC_LWP_exclusive_cells, dims.use = 1:25, do.fast = T)
TSNEPlot(LWC_LWP_exclusive_cells, do.label = TRUE, pt.size = 0.8)

sample_1 <- grep("-1", colnames(LWC_LWP_exclusive_cells@data), value = T); length(sample_1)
sample_2 <- grep("-2", colnames(LWC_LWP_exclusive_cells@data), value = T); length(sample_2)
sample_3 <- grep("-3", colnames(LWC_LWP_exclusive_cells@data), value = T); length(sample_3)
sample_4 <- grep("-4", colnames(LWC_LWP_exclusive_cells@data), value = T); length(sample_4)



jpeg(file = "Large_Wound_CvsP_exclusive_cells_tsne.jpeg", width = 15, height = 15, units = "cm", res = 500)
s_tSNE <- as.data.frame(GetCellEmbeddings(LWC_LWP_exclusive_cells, reduction.type = "tsne"))
s_tSNE[sample_1, "Sample"] <- "LWC D14"
s_tSNE[sample_2, "Sample"] <- "LWC D14"
s_tSNE[sample_3, "Sample"] <- "LWO D14"
s_tSNE[sample_4, "Sample"] <- "LWO D14"
ggplot(s_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample),size = 0.4) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()



current.cluster.ids <- c(0, 1, 2, 3, 4,5,6,7,8,9)
new.cluster.ids <- c("LWP",
                     "LWC", 
                     "LWC",
                     "Mixed",
                     "LWC",
                     "Mixed",
                     "LWP",
                     "LWC",
                     "Mixed",
                     "LWC",
                     "LWC",
                     "Mixed",
                     "LWP",
                     "Mixed",
                     "Mixed",
                     "Mixed")
LWD14fibr_subset@ident <- plyr::mapvalues(x = LWD14fibr_subset@ident, from = current.cluster.ids, to = new.cluster.ids)

TSNEPlot(LWC_LWP_exclusive_cells, do.label = TRUE, pt.size = 0.8)




FeaturePlot(LWC_LWP_exclusive_cells, features.plot = c("Runx1"), reduction.use = "tsne", cols.use = c("yellow", "red"))



LWC_LWP_exclusive_cells = SubsetData(LWC_LWP_exclusive_cells, subset.raw = TRUE)
save(LWC_LWP_exclusive_cells, file ="LWC_LWP_exclusive_cells.Robj")





####################################### END #######



LWD14fibr_subset10.markers <- FindMarkers(LWD14fibr_subset, ident.1 = c(10), min.pct = 0.25, test.use = "negbinom")
LWD14fibr_subset10.markers <- LWD14fibr_subset10.markers[order(LWD14fibr_subset10.markers$avg_logFC),]
write.csv(LWD14fibr_subset10.markers,file = "LWD14fibr_subset10.markers.csv")




####### Redoing Fibroblast Trajectory using tdTom Dataset  #######
wounds_D14 <- Read10X(data.dir = "/Users/SarthakSinha/Desktop/10X Bioinformatics/filtered_gene_bc_matrices_mex_d14/cellranger_mouse_genome_pluscretomssm2gfp201")
dense.size <- object.size(x = as.matrix(x = wounds_D14))
sparse.size <- object.size(x = wounds_D14)
dense.size/sparse.size
wound_D14 <- CreateSeuratObject(raw.data = wounds_D14, min.cells = 2, min.genes = 200, project = "10X_Wound")
mito.genes <- grep(pattern = "^mt-", x = rownames(x = wound_D14@data), value = TRUE)
percent.mito <- Matrix::colSums(wound_D14@raw.data[mito.genes, ])/Matrix::colSums(wound_D14@raw.data)
wound_D14 <- AddMetaData(object = wound_D14, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = wound_D14, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = wound_D14, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = wound_D14, gene1 = "nUMI", gene2 = "nGene")

wound_D14 <- FilterCells(object = wound_D14, subset.names = c("nGene", "percent.mito"), low.thresholds = c(300, -Inf), high.thresholds = c(6000, 0.6))
wound_D14 <- FilterCells(object = wound_D14, subset.names = c("nUMI"), low.thresholds = c(-Inf), high.thresholds = c(60000))
wound_D14 <- NormalizeData(object = wound_D14, normalization.method = "LogNormalize",scale.factor = 10000)
wound_D14 <- FindVariableGenes(object = wound_D14, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = wound_D14@var.genes)
wound_D14 <- ScaleData(object = wound_D14, vars.to.regress = c("nUMI", "percent.mito"))
wound_D14 <- RunPCA(object = wound_D14, pc.genes = wound_D14@var.genes, do.print = TRUE, pcs.print = 1:5,genes.print = 5)
PrintPCA(object = wound_D14, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
wound_D14 <- RunPCA(wound_D14, pc.genes = wound_D14@var.genes, pcs.compute = 200, do.print = T, pcs.print = 5, genes.print = 5)
PCElbowPlot(wound_D14, num.pc = 50)
wound_D14 <- FindClusters(wound_D14, reduction.type = "pca", dims.use = 1:25, resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
wound_D14 <- RunTSNE(wound_D14, dims.use = 1:25, do.fast = T)
TSNEPlot(wound_D14, do.label = TRUE, pt.size = 0.8)

save(wound_D14, file = "wound_D14.Robj")

sample_1 <- grep("-1", colnames(wound_D14@data), value = T); length(sample_1)
sample_2 <- grep("-2", colnames(wound_D14@data), value = T); length(sample_2)
sample_3 <- grep("-3", colnames(wound_D14@data), value = T); length(sample_3)
sample_4 <- grep("-4", colnames(wound_D14@data), value = T); length(sample_4)

FeaturePlot(fibro_recluster, features.plot = c("tdTomato"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(fibro_recluster, features.plot = c("Pdgfra"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(fibro_recluster, features.plot = c("Dpt"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(fibro_recluster, features.plot = c("Slc27a1"), reduction.use = "tsne", cols.use = c("yellow", "red"))

FeaturePlot(fibro_recluster, features.plot = c("Sox18"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(fibro_recluster, features.plot = c("Sox18"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10", max.cutoff = "q90")



fibro <- names(wound_D14@ident[wound_D14@ident %in% c(0,7,14,3,13, 20, 16,19,6)])
sub.fibro <- FilterCells(wound_D14, subset.names = "nGene", low.thresholds = 50, high.thresholds = 2000000, cells.use =fibro)
sub.fibro <- FindVariableGenes(sub.fibro, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
length(x = sub.fibro@var.genes)
fibro_recluster <- NormalizeData(sub.fibro, normalization.method = "LogNormalize", scale.factor = 10000)
fibro_recluster <- ScaleData(fibro_recluster, vars.to.regress = c("nUMI", "percent.mito"))
fibro_recluster <- RunPCA(fibro_recluster, pc.genes = fibro_recluster@var.genes, pcs.compute = 50)
PCElbowPlot(fibro_recluster, num.pc = 50)
fibro_recluster <- FindClusters(fibro_recluster, reduction.type = "pca", dims.use = 1:25, resolution = 1.5, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
fibro_recluster <- RunTSNE(fibro_recluster, dims.use = 1:25, do.fast = T)
TSNEPlot(fibro_recluster, do.label = TRUE, pt.size = 0.8)
save(fibro_recluster, file = "fibro_recluster_D14_tdT.Robj")


fibro_recluster_tdT <- RunDiffusion(fibro_recluster,genes.use = fibro_recluster@var.genes)
fibro_recluster <- RunDiffusion(fibro_recluster,genes.use = fibro_recluster@var.genes, dims.use = 1:10)
fibro_recluster_tdT <- RunDiffusion(fibro_recluster_tdT,genes.use = fibro_recluster_tdT@var.genes, dims.use = 1:10, q.use = 0.01, max.dim = 5)
TSNEPlot(fibro_recluster, do.label = TRUE, pt.size = 0.8)

FeaturePlot(fibro_recluster, features.plot = c("Crabp1"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(fibro_recluster, features.plot = c("S100a4"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(fibro_recluster, features.plot = c("Cd200"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(fibro_recluster, features.plot = c("Sox18"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(fibro_recluster, features.plot = c("Rspo3"), reduction.use = "tsne", cols.use = c("yellow", "red"))

fibro <- names(fibro_recluster@ident[fibro_recluster@ident %in% c(0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)])
sub.fibro <- FilterCells(fibro_recluster, subset.names = "nGene", low.thresholds = 0, high.thresholds = 2000000, cells.use =fibro)
sub.fibro <- FindVariableGenes(sub.fibro, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
length(x = sub.fibro@var.genes)
fibro_recluster <- NormalizeData(sub.fibro, normalization.method = "LogNormalize", scale.factor = 10000)
fibro_recluster <- ScaleData(fibro_recluster, vars.to.regress = c("nUMI", "percent.mito"))
fibro_recluster <- RunPCA(fibro_recluster, pc.genes = fibro_recluster@var.genes, pcs.compute = 50)
PCElbowPlot(fibro_recluster, num.pc = 50)
fibro_recluster <- FindClusters(fibro_recluster, reduction.type = "pca", dims.use = 1:25, resolution = 1.5, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
fibro_recluster <- RunTSNE(fibro_recluster, dims.use = 1:25, do.fast = T)
TSNEPlot(fibro_recluster, do.label = TRUE, pt.size = 0.8)
save(fibro_recluster, file = "fibro_recluster_D14_tdT.Robj")

fibro_recluster <- RunDiffusion(fibro_recluster,genes.use = fibro_recluster@var.genes)

DMPlot(fibro_recluster, do.label = TRUE)
FeaturePlot(fibro_recluster, features.plot = c("S100a4", "Cd200", "Sox18", "Crabp1"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(fibro_recluster, features.plot = c("tdTomato"), reduction.use = "tsne", cols.use = c("yellow", "red"))







sample_1 <- grep("-1", colnames(fibro_recluster@data), value = T); length(sample_1)
sample_2 <- grep("-2", colnames(fibro_recluster@data), value = T); length(sample_2)
sample_3 <- grep("-3", colnames(fibro_recluster@data), value = T); length(sample_3)
sample_4 <- grep("-4", colnames(fibro_recluster@data), value = T); length(sample_4)

s_tSNE <- as.data.frame(GetCellEmbeddings(fibro_recluster, reduction.type = "tsne"))
s_tSNE[sample_1, "Sample"] <- "LWC D14"
s_tSNE[sample_2, "Sample"] <- "LWC D14"
s_tSNE[sample_3, "Sample"] <- "LWO D14"
s_tSNE[sample_4, "Sample"] <- "LWO D14"
ggplot(s_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample),size = 0.3) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()

FeaturePlot(fibro_recluster, features.plot = c("Sox18"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(fibro_recluster, features.plot = c("Rspo3"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(fibro_recluster, features.plot = c("tdTomato"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(fibro_recluster, features.plot = c("S100a4"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(fibro_recluster, features.plot = c("Crabp1"), reduction.use = "tsne", cols.use = c("yellow", "red"))
FeaturePlot(fibro_recluster, features.plot = c("Cd200"), reduction.use = "tsne", cols.use = c("yellow", "red"))




library(RColorBrewer)
jpeg(file = "SCENIC_Contour.jpeg", width = 20, height = 20, units = "cm", res = 500)
dens2d <- bkde2D(tSNE_scenic$Y, 1)$fhat
image(dens2d, col=brewer.pal(10, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=12, drawlabels=F)
dev.off()

jpeg(file = "SCENIC_Creb3l1_Gli1_colocalize_cex_pos.jpeg", width = 20, height = 20, units = "cm", res = 500)
regulonNames <- c( "Creb3l1","Gli1")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.8, offColor="black")
text(-30,-25, attr(cellCol,"red"), col="red", cex=1.4, pos=2)
text(-30,-25-4, attr(cellCol,"blue"), col="blue", cex=1.4, pos=2)
dev.off()

jpeg(file = "SCENIC_Creb3l1_Irf7_colocalize.jpeg", width = 20, height = 20, units = "cm", res = 500)
regulonNames <- c( "Creb3l1","Irf7")
cellCol <- plotTsne_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)
text(-30,-25, attr(cellCol,"red"), col="red", cex=.7, pos=4)
text(-30,-25-4, attr(cellCol,"green"), col="green3", cex=.7, pos=4)
dev.off()



jpeg(file = "all_fibros_tSNE_seperated_by_sample_of_interest.jpeg", width = 15, height = 15, units = "cm", res = 500)
s_tSNE <- as.data.frame(GetCellEmbeddings(all_fibros, reduction.type = "tsne"))
s_tSNE[sample_1, "Sample"] <- "LWC D14"
s_tSNE[sample_2, "Sample"] <- "LWC D14"
s_tSNE[sample_3, "Sample"] <- "LWP D14"
s_tSNE[sample_4, "Sample"] <- "LWP D14"
s_tSNE[sample_5, "Sample"] <- "SW D14"
s_tSNE[sample_6, "Sample"] <- "SW D14"
s_tSNE[Day18_Day8_Lib3_try1, "Sample"] <- "SW D8"
s_tSNE[Day18_Day8_Lib4_try1, "Sample"] <- "SW D8"
s_tSNE[sample_7, "Sample"] <- "UI"
ggplot(s_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample),size = 0.3) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()




###### RNA Probe Verification
FeaturePlot(LWD14fibr_subset, features.plot = c("Prss35"), reduction.use = "tsne", cols.use = c("yellow", "red"))
VlnPlot(LWD14fibr_subset, features.plot = c("Prss35"))
FeaturePlot(LWC_D14_fibr_trajec_new, features.plot = c("Crabp1"), reduction.use = "dm", cols.use = c("yellow", "red"))
FeaturePlot(D14LWC_LWP_SW8_UI, features.plot = c("Crabp1","Prss35"), reduction.use = "tsne", cols.use = c("yellow", "red"))



#Hexb




#####
save(LWC_D14_fibr, file = "LWC_D14_fibr.latest.Robj")
jpeg(file = "LWC_D14_fibr_Prss35.jpeg", width = 15, height = 15, units = "cm", res = 500)
VlnPlot(LWC_D14_fibr, features = c("Prss35"), do.sort = FALSE, ident = c("Lower Dermis LWC", "Upper Dermis LWC"), cols = c("#28B1B7", "#E95F5C"))
dev.off()



