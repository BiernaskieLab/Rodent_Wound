#Running SCENIC on Day 14 LW vs Day 8 SW Cells: # RAN ON SYNERGY:

library(Seurat)
load(file = "wound_subset.Robj")

jpeg(file = "wound_subset_Dpt_Pdgfra.jpeg", width = 15, height = 15, units = "cm", res = 500)
FeaturePlot(wound_subset, features.plot = c("Dpt", "Pdgfra"), reduction.use = "tsne", cols.use = c("yellow", "red"))
dev.off()

jpeg(file = "wound_subset_Dpt_Pdgfra.jpeg", width = 25, height = 15, units = "cm", res = 500)
FeaturePlot(wound_subset, features.plot = c("Dpt", "Pdgfra", "Rspo3", "Sox18"), reduction.use = "tsne", cols.use = c("yellow", "red"))
dev.off()

# Clusters: 12,5,13,0,9,20,3,7

jpeg(file = "tsne_overall.jpeg", width = 15, height = 15, units = "cm", res = 500)
TSNEPlot(wound_subset, do.label = TRUE, pt.size = 0.8); dev.off()

all_fibros <- names(wound_subset@ident[wound_subset@ident %in% c(12,5,13,0,9,20,3,7)])
all_fibros <- FilterCells(wound_subset, subset.names = "nGene", low.thresholds = 1, high.thresholds = 200000, cells.use = all_fibros)
all_fibros <- FindVariableGenes(all_fibros, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
length(x = all_fibros@var.genes) #1290
all_fibros <- NormalizeData(all_fibros, normalization.method = "LogNormalize", scale.factor = 10000)
all_fibros <- ScaleData(all_fibros, vars.to.regress = c("nUMI", "percent.mito"))
all_fibros <- RunPCA(all_fibros, pc.genes = all_fibros@var.genes, pcs.compute = 50)
jpeg(file = "all_fibros_elbow.jpeg", width = 15, height = 15, units = "cm", res = 500)
PCElbowPlot(all_fibros, num.pc = 50); dev.off()
all_fibros <- FindClusters(all_fibros, reduction.type = "pca", dims.use = 1:25, resolution = 0.6, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
all_fibros <- RunTSNE(all_fibros, dims.use = 1:25, do.fast = T)
all_fibros = SubsetData(all_fibros, subset.raw = TRUE)
save(all_fibros, file ="all_fibros.Robj")

jpeg(file = "all_fibros_FP_regen_markers.jpeg", width = 25, height = 20, units = "cm", res = 500)
FeaturePlot(all_fibros, features.plot = c("Crabp1", "Rspo3", "Sox18", "S100a4"), reduction.use = "tsne", cols.use = c("yellow", "red")); dev.off()




sample_1 <- grep("-1", colnames(all_fibros@data), value = T); length(sample_1)
sample_2 <- grep("-2", colnames(all_fibros@data), value = T); length(sample_2)
sample_3 <- grep("-3", colnames(all_fibros@data), value = T); length(sample_3)
sample_4 <- grep("-4", colnames(all_fibros@data), value = T); length(sample_4)
sample_5 <- grep("-5", colnames(all_fibros@data), value = T); length(sample_5)
sample_6 <- grep("-6", colnames(all_fibros@data), value = T); length(sample_6)
sample_7 <- grep("-7", colnames(all_fibros@data), value = T); length(sample_7)
Day18_Day8_Lib1_try1 <- grep("-8", colnames(all_fibros@data), value = T); length(Day18_Day8_Lib1_try1)
Day18_Day8_Lib2_try1 <- grep("-9", colnames(all_fibros@data), value = T); length(Day18_Day8_Lib2_try1)
Day18_Day8_Lib3_try1 <- grep("-10", colnames(all_fibros@data), value = T); length(Day18_Day8_Lib3_try1)
Day18_Day8_Lib4_try1 <- grep("-11", colnames(all_fibros@data), value = T); length(Day18_Day8_Lib4_try1)

jpeg(file = "all_fibros_tSNE_seperated_by_sample_of_interest_1.jpeg", width = 15, height = 15, units = "cm", res = 500)
s_tSNE <- as.data.frame(GetCellEmbeddings(all_fibros, reduction.type = "tsne"))
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




# Seurat Obj 1:
D14LWC_LWP_SW8_UI_1 <- FilterCells(all_fibros, cells.use = c(sample_1, sample_2, sample_3, sample_4,sample_5,sample_6, sample_7), subset.names = "nGene", low.thresholds = 1, high.thresholds = 200000)
# Seurat Obj 2:
D14LWC_LWP_SW8_UI_2 <- FilterCells(all_fibros, cells.use = c(Day18_Day8_Lib3_try1, Day18_Day8_Lib4_try1), subset.names = "nGene", low.thresholds = 1, high.thresholds = 200000)
D14LWC_LWP_SW8_UI = MergeSeurat(D14LWC_LWP_SW8_UI_1, D14LWC_LWP_SW8_UI_2, add.cell.id1 = "r_1", add.cell.id2 = "r_2")

D14LWC_LWP_SW8_UI <- FindVariableGenes(D14LWC_LWP_SW8_UI, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
length(x = D14LWC_LWP_SW8_UI@var.genes) #1290
D14LWC_LWP_SW8_UI <- NormalizeData(D14LWC_LWP_SW8_UI, normalization.method = "LogNormalize", scale.factor = 10000)
D14LWC_LWP_SW8_UI <- ScaleData(D14LWC_LWP_SW8_UI, vars.to.regress = c("nUMI", "percent.mito"))
D14LWC_LWP_SW8_UI <- RunPCA(D14LWC_LWP_SW8_UI, pc.genes = D14LWC_LWP_SW8_UI@var.genes, pcs.compute = 50)
jpeg(file = "D14LWC_LWP_SW8_UI_elbow.jpeg", width = 15, height = 15, units = "cm", res = 500)
PCElbowPlot(D14LWC_LWP_SW8_UI, num.pc = 50); dev.off()
D14LWC_LWP_SW8_UI <- FindClusters(D14LWC_LWP_SW8_UI, reduction.type = "pca", dims.use = 1:25, resolution = 0.6, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
D14LWC_LWP_SW8_UI <- RunTSNE(D14LWC_LWP_SW8_UI, dims.use = 1:25, do.fast = T, check_duplicates = FALSE)
jpeg(file = "tsne_D14LWC_LWP_SW8_UI.jpeg", width = 15, height = 15, units = "cm", res = 500)
TSNEPlot(D14LWC_LWP_SW8_UI, do.label = TRUE, pt.size = 0.8); dev.off()

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



jpeg(file = "D14LWC_LWP_SW8_UI_sample_of_interest_1.jpeg", width = 15, height = 15, units = "cm", res = 500)
s_tSNE <- as.data.frame(GetCellEmbeddings(D14LWC_LWP_SW8_UI, reduction.type = "tsne"))
s_tSNE[sample_1, "Sample"] <- "LW D14"
s_tSNE[sample_2, "Sample"] <- "LW D14"
s_tSNE[sample_3, "Sample"] <- "LW D14"
s_tSNE[sample_4, "Sample"] <- "LW D14"
s_tSNE[sample_5, "Sample"] <- "SW D14"
s_tSNE[sample_6, "Sample"] <- "SW D14"
s_tSNE[Day18_Day8_Lib3_try1, "Sample"] <- "SW D8"
s_tSNE[sample_7, "Sample"] <- "UI"
ggplot(s_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample),size = 0.3) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()




save(D14LWC_LWP_SW8_UI, file ="D14LWC_LWP_SW8_UI.Robj")


jpeg(file = "D14LWC_LWP_SW8_UI_FP_regen_markers.jpeg", width = 25, height = 20, units = "cm", res = 500)
FeaturePlot(D14LWC_LWP_SW8_UI, features.plot = c("Crabp1", "Rspo3", "Sox18", "S100a4", "Igfbp4"), reduction.use = "tsne", cols.use = c("yellow", "red")); dev.off()


###################################### START FROM HERE
####### Importing all the fibroblast data onto local drive: #######
load(file = "D14LWC_LWP_SW8_UI.Robj")
rm(D14LWC_LWP_SW8_UI)
TSNEPlot(D14LWC_LWP_SW8_UI, do.label = TRUE, pt.size = 0.3)


LWD14fibr_subset_CC <- CellCycleScoring(LWD14fibr_subset, s.genes = s.genes, 
                                        g2m.genes = g2m.genes, set.ident = TRUE)

TSNEPlot(LWD14fibr_subset_CC, do.label = F, pt.size = 0.8)
DimPlot(LWD14fibr_subset_CC, reduction.use = 'phate')




