
tmux
bsub -Is bash
source activate Seurat_R
R
library(Seurat)
library(dplyr)

##### RUN ON Syngergy:
LW_R <- Read10X(data.dir = "/home/ssinha/10X_Cell_Ranger/agg_tdTCustomRefGenome_ALL_samples/outs/filtered_gene_bc_matrices_mex/cellranger_mouse_genome_pluscretomssm2gfp201")
dense.size <- object.size(x = as.matrix(x = LW_R))
sparse.size <- object.size(x = LW_R)
dense.size/sparse.size

wound_tdT <- CreateSeuratObject(raw.data = LW_R, min.cells = 3, min.genes = 300, 
                            project = "wound_tdT")

mito.genes <- grep(pattern = "^mt-", x = rownames(x = wound_tdT@data), value = TRUE)
percent.mito <- Matrix::colSums(wound_tdT@raw.data[mito.genes, ])/Matrix::colSums(wound_tdT@raw.data)
wound_tdT <- AddMetaData(wound_tdT, metadata = percent.mito, col.name = "percent.mito")

jpeg(file = "QC_wound_tdT.jpeg", width = 20, height = 15, units = "cm", res = 500)
VlnPlot(wound_tdT, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

wound_tdT <- FilterCells(wound_tdT, subset.names = c("nGene", "percent.mito"), 
                     low.thresholds = c(300, -Inf), high.thresholds = c(6000, 0.6))
wound_tdT <- FilterCells(wound_tdT, subset.names = c("nUMI"), 
                     low.thresholds = c(-Inf), high.thresholds = c(60000))

jpeg(file = "QC_wound_tdT_GenePlots.jpeg", width = 20, height = 15, units = "cm", res = 500)
par(mfrow = c(1, 2))
GenePlot(object = wound_tdT, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = wound_tdT, gene1 = "nUMI", gene2 = "nGene")
dev.off()

wound_tdT <- NormalizeData(wound_tdT, normalization.method = "LogNormalize", scale.factor = 10000)
wound_tdT <- FindVariableGenes(wound_tdT, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = F)
length(x = wound_tdT@var.genes) #1290

wound_tdT <- ScaleData(wound_tdT, vars.to.regress = c("nUMI", "percent.mito"))
wound_tdT <- RunPCA(wound_tdT, pc.genes = wound_tdT@var.genes, pcs.compute = 50)

jpeg(file = "QC_wound_tdT_PCElbowPlot.jpeg", width = 20, height = 15, units = "cm", res = 500)
PCElbowPlot(wound_tdT, num.pc = 50)
dev.off()

wound_tdT <- FindClusters(wound_tdT, reduction.type = "pca", dims.use = 1:28, resolution = 0.6, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
wound_tdT <- RunTSNE(wound_tdT, dims.use = 1:25, do.fast = T, check_duplicates = F)

jpeg(file = "TSNEPlot_wound_tdT.jpeg", width = 20, height = 15, units = "cm", res = 500)
TSNEPlot(wound_tdT, do.label = TRUE, pt.size = 0.8)
dev.off()

save(wound_tdT, file = "wound_tdT.Robj")

############################# CELL CYCLE ANALYSIS: ############################
all_fibros_4_scenic_CC <- CellCycleScoring(all_fibros_4_scenic, s.genes = s.genes, 
                                         g2m.genes = g2m.genes, set.ident = TRUE)

TSNEPlot(all_fibros_4_scenic_CC, do.label = F, pt.size = 0.3)


all_fibros_4_scenic_CC_analysis = all_fibros_4_scenic_CC@meta.data
write.csv(all_fibros_4_scenic_CC_analysis,file = "all_fibros_4_scenic_CC_analysis.csv")
############# Creating a stacked bar chart to display this:

#Create data
data_1 <- read.csv("woundtypes_CC_analysis_processed.csv", stringsAsFactors = FALSE, header = F)
colnames(data_1)=c("UI","SWD8","SWD14","LWD14")
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





load(file = "all_fibros_4_scenic.Robj")
all_fibros_4_scenic_v3 = UpdateSeuratObject(all_fibros_4_scenic)

all_fibros_4_scenic_v3[["percent.mt"]] <- PercentageFeatureSet(all_fibros_4_scenic_v3, pattern = "^mt-")

View(all_fibros_4_scenic_v3@meta.data)

all_fibros_4_scenic_v3[["percent.mt"]] <- PercentageFeatureSet(all_fibros_4_scenic_v3, pattern = "mt-")

Vlnplot_gg = VlnPlot(all_fibros_4_scenic_v3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

VlnPlot(all_fibros_4_scenic_v3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, cols = c("#6EA42A", "#E95F5C", "#28B1B7", "#9668AF"), idents = c("UI", "SW D8", "LW D14", "SW D14"), sort = F, adjust = 1.5)

jpeg(file = "QC_wound_tdT_fibros.jpeg", width = 20, height = 15, units = "cm", res = 500)
VlnPlot(all_fibros_4_scenic_v3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, cols = c("#6EA42A", "#E95F5C", "#28B1B7", "#9668AF"), idents = c("UI", "SW D8", "LW D14", "SW D14"), sort = F, adjust = 2)
dev.off()


library(RColorBrewer)
library(ggplot2)
myColors <- c("#6EA42A", "#E95F5C", "#28B1B7", "#9668AF")
names(myColors) <- levels(all_fibros_4_scenic_v3@active.ident)
colScale <- scale_colour_manual(values = myColors)
Vlnplot_gg_full <- Vlnplot_gg + colScale
Vlnplot_gg_full

View(all_fibros_4_scenic_v3@active.ident)

plot1 <- FeatureScatter(all_fibros_4_scenic_v3, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(all_fibros_4_scenic_v3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))




all_fibros_4_scenic_v3_subsetted <- subset(all_fibros_4_scenic_v3, percent.mt < 10)


TSNEPlot(all_fibros_4_scenic_v3, do.label = TRUE, pt.size = 0.3)
TSNEPlot(all_fibros_4_scenic_v3_subsetted, do.label = TRUE, pt.size = 0.3)


FeaturePlot(all_fibros_4_scenic, features.plot = c("Runx1"), reduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q90", min.cutoff = "q20")
VlnPlot(all_fibros_4_scenic, features.plot = c("Runx1"))

