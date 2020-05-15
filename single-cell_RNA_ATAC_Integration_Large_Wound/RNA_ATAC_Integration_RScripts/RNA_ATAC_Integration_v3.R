

library(Seurat)
library(ggplot2)
library(Matrix)

peaks <- Read10X_h5("filtered_peak_bc_matrix.h5")
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = "Mus_musculus.GRCm38.96.gtf", 
                                            seq.levels = c(1:19, "X", "Y"), upstream = 2000, verbose = TRUE)
fibro_center_atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "10x_ATAC")

fibro_center_atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
meta <- read.table("singlecell.csv", sep = ",", header = TRUE, row.names = 1, 
                   stringsAsFactors = FALSE)
meta <- meta[colnames(fibro_center_atac), ]
fibro_center_atac <- AddMetaData(fibro_center_atac, metadata = meta)
fibro_center_atac <- subset(fibro_center_atac, subset = nCount_ATAC > 1000)
fibro_center_atac$tech <- "atac"

View(fibro_center_atac@meta.data)

DefaultAssay(fibro_center_atac) <- "ACTIVITY"
fibro_center_atac <- FindVariableFeatures(fibro_center_atac)
fibro_center_atac <- NormalizeData(fibro_center_atac)
fibro_center_atac <- ScaleData(fibro_center_atac)


DefaultAssay(fibro_center_atac) <- "ATAC"
VariableFeatures(fibro_center_atac) <- names(which(Matrix::rowSums(fibro_center_atac) > 100))
ElowPlot(fibro_center_atac)
fibro_center_atac <- RunLSI(fibro_center_atac, n = 20, scale.max = NULL)
fibro_center_atac <- RunUMAP(fibro_center_atac, reduction = "lsi", dims = 1:20)

load(file = "fibro_center_rna.Robj")
fibro_center_rna$tech <- "rna"
fibro_center_rna$celltype = fibro_center_rna@active.ident


p1 <- DimPlot(fibro_center_atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
p2 <- DimPlot(LW_center_rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))


load(file = "LW_center_rna.Robj")
load(file = "fibro_center_atac_LATEST_peakcorrected.Robj")

LW_center_rna$tech <- "rna"
LW_center_rna$celltype = LW_center_rna@active.ident


p1 <- DimPlot(fibro_center_atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
p2 <- DimPlot(LW_center_rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))


FeaturePlot(fibro_center_atac, features = "chr5:75149771-75157509", max.cutoff = 500) # Pdgfra
FeaturePlot(fibro_center_atac, features = "Crabp1", max.cutoff = 500) # Pdgfra
FeaturePlot(fibro_center_atac, features = "Rspo3", max.cutoff = 500) # Pdgfra
FeaturePlot(fibro_center_atac, features = "Sox18", max.cutoff = 500) # Pdgfra


#### Current:
DimPlot(fibro_center_atac, reduction = "umap")

fibro_center_atac <- ScaleData(fibro_center_atac, features = all.genes)
fibro_center_atac <- RunPCA(fibro_center_atac, features = VariableFeatures(object = fibro_center_atac))
fibro_center_atac <- FindNeighbors(fibro_center_atac, dims = 1:20)
fibro_center_atac <- FindClusters(fibro_center_atac, resolution = 0.5)





View(fibro_center_atac@assays$ATAC@counts@Dimnames)

save(fibro_center_atac, file = "fibro_center_atac_LATEST_peakcorrected.Robj")


ta <- FindTransferAnchors(reference = LW_center_rna, query = fibro_center_atac, features = VariableFeatures(object = LW_center_rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

View(ta)
LW_center_rna$celltype = LW_center_rna@active.ident
celltype.predictions <- TransferData(anchorset = ta, refdata = LW_center_rna@active.ident, 
                                     weight.reduction = fibro_center_atac[["lsi"]])



### Trying with fibro_center_rna:
ta <- FindTransferAnchors(reference = fibro_center_rna, 
                          query = fibro_center_atac, 
                          features = VariableFeatures(object = fibro_center_rna), 
                          reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")


fibro_center_rna$celltype = fibro_center_rna@active.ident
celltype.predictions <- TransferData(anchorset = ta, refdata = fibro_center_rna@active.ident, 
                                     weight.reduction = fibro_center_atac[["lsi"]])



### Trying to make the two datasets - rna and atac - comparable!
fibro_center_atac <- NormalizeData(fibro_center_atac, normalization.method = "LogNormalize", scale.factor = 10000)
fibro_center_atac <- FindVariableFeatures(fibro_center_atac, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(fibro_center_atac)
all.genes <- rownames(fibro_center_atac)
fibro_center_atac <- ScaleData(fibro_center_atac)
fibro_center_atac <- RunPCA(fibro_center_atac, features = VariableFeatures(object = fibro_center_atac))
ElbowPlot(fibro_center_atac)
fibro_center_atac <- FindNeighbors(fibro_center_atac, dims = 1:20)
fibro_center_atac <- FindClusters(fibro_center_atac, resolution = 0.5)
fibro_center_atac <- RunUMAP(fibro_center_atac, dims = 1:20)
DimPlot(fibro_center_atac, reduction = "umap")

FeaturePlot(fibro_center_atac, features = "Pdgfra", max.cutoff = 500)
FeaturePlot(fibro_center_atac, features = "Crabp1", max.cutoff = 500)
FeaturePlot(fibro_center_atac, features = "Rspo3", max.cutoff = 500)



fibro_center_atac <- AddMetaData(fibro_center_atac, metadata = meta)
fibro_center_atac <- subset(fibro_center_atac, subset = nCount_ATAC > 2500)
fibro_center_atac$tech <- "atac"
DefaultAssay(fibro_center_atac) <- "ACTIVITY"
fibro_center_atac <- FindVariableFeatures(fibro_center_atac)
fibro_center_atac <- NormalizeData(fibro_center_atac)
fibro_center_atac <- ScaleData(fibro_center_atac)
DefaultAssay(fibro_center_atac) <- "ATAC"
VariableFeatures(fibro_center_atac) <- names(which(Matrix::rowSums(fibro_center_atac) > 100))
fibro_center_atac <- RunLSI(fibro_center_atac, n = 20, scale.max = NULL)
fibro_center_atac <- RunUMAP(fibro_center_atac, reduction = "lsi", dims = 1:20)


pdf(file = "fibro_center_atac_pdgfra")
DimPlot(fibro_center_atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
FeaturePlot(fibro_center_atac, features = "Pdgfra")
dev.off()


fibro_center_rna
fibro_center_rna$tech <- "rna"
p1 <- DimPlot(fibro_center_atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
p2 <- DimPlot(fibro_center_rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))

all.genes <- rownames(fibro_center_atac)
fibro_center_atac <- ScaleData(fibro_center_atac, features = all.genes)
fibro_center_atac <- RunPCA(fibro_center_atac, features = VariableFeatures(object = fibro_center_atac))
fibro_center_atac <- FindNeighbors(fibro_center_atac, dims = 1:20)
fibro_center_atac <- FindClusters(fibro_center_atac, resolution = 0.5)

DimPlot(fibro_center_atac, reduction = "umap")+ ggtitle("scATAC-seq")
FeaturePlot(fibro_center_atac, features = "Pdgfra")
FeaturePlot(fibro_center_atac, features = "Dpt")



p1 <- DimPlot(fibro_center_atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
p2 <- DimPlot(fibro_center_rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))


# save(fibro_center_rna, file = "fibro_center_rna.Robj")
# save(fibro_center_atac, file = "fibro_center_atac.Robj")
# 
# fibro_center_rna$tech <- "rna"
# fibro_center_atac <- RunLSI(fibro_center_atac, n = 50, scale.max = NULL)
transfer.anchors <- FindTransferAnchors(reference = fibro_center_rna,
                                        query = fibro_center_atac, 
                                        features = VariableFeatures(object = fibro_center_rna), 
                                        reference.assay = "RNA", 
                                        query.assay = "ACTIVITY", 
                                        reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = fibro_center_rna@active.ident, 
                                     weight.reduction = fibro_center_atac[["lsi"]])

fibro_center_atac <- AddMetaData(fibro_center_atac, metadata = celltype.predictions)
hist(fibro_center_atac$prediction.score.max)
abline(v = 0.5, col = "red")


table(fibro_center_atac$prediction.score.max > 0.5)

fibro_center_atac.filtered <- subset(fibro_center_atac, subset = prediction.score.max > 0.5)
fibro_center_atac.filtered$predicted.id <- factor(fibro_center_atac.filtered$predicted.id, levels = levels(fibro_center_rna))  # to make the colors match
p1 <- DimPlot(fibro_center_atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(fibro_center_rna, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
  NoLegend()
CombinePlots(plots = list(p1, p2))


fibro_center_atac.filtered <- subset(fibro_center_atac, subset = prediction.score.max > 0.5)
fibro_center_atac.filtered$predicted.id <- factor(fibro_center_atac.filtered$predicted.id, levels = levels(fibro_center_rna))  # to make the colors match
p1 <- DimPlot(fibro_center_atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(fibro_center_rna, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
  NoLegend()
CombinePlots(plots = list(p1, p2))


DimPlot(fibro_center_atac, reduction = "umap") + ggtitle("scATAC-seq")


fibro_center_atac <- subset(fibro_center_atac, subset = prediction.score.max > 0.5)
fibro_center_atac <- RunLSI(fibro_center_atac, n = 50, scale.max = NULL)
fibro_center_atac <- RunUMAP(fibro_center_atac, reduction = "lsi", dims = 1:50)


transfer.anchors <- FindTransferAnchors(reference = fibro_center_rna, 
                                        query = fibro_center_atac, 
                                        features = VariableFeatures(object = fibro_center_rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")


celltype.predictions <- TransferData(anchorset = transfer.anchors, 
                                     refdata = fibro_center_rna$celltype, 
                                     weight.reduction = fibro_center_atac[["lsi"]])
fibro_center_atac <- AddMetaData(fibro_center_atac, metadata = celltype.predictions)

hist(fibro_center_atac$prediction.score.max)
abline(v = 0.5, col = "red")
table(fibro_center_atac$prediction.score.max > 0.5)


fibro_center_atac.filtered <- subset(fibro_center_atac, subset = prediction.score.max > 0.5)
fibro_center_atac.filtered$predicted.id <- factor(fibro_center_atac.filtered$predicted.id, levels = levels(fibro_center_rna))  # to make the colors match
p1 <- DimPlot(fibro_center_atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(fibro_center_rna, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
  NoLegend()
CombinePlots(plots = list(p1, p2))





genes.use <- VariableFeatures(fibro_center_rna)
refdata <- GetAssayData(fibro_center_rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = fibro_center_atac[["lsi"]])

# this line adds the imputed data matrix to the pbmc.atac object
fibro_center_atac[["RNA"]] <- imputation


coembed <- merge(x = fibro_center_rna, y = fibro_center_atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
ElbowPlot(coembed)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)

p1 <- DimPlot(coembed, group.by = "tech")
p2 <- DimPlot(coembed, group.by = "celltype", repel = T)
CombinePlots(list(p1, p2))


save(coembed, file = "coembed_updated_v3.Robj")
load(file = "coembed_updated_v3.Robj")

jpeg(file = "coembed_groupByCellType.jpeg", width = 20, height = 15, units = "cm", res = 500)
DimPlot(coembed, group.by = "celltype", repel = T)
dev.off()


DimPlot(coembed, split.by = "tech", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()



RidgePlot(coembed, features = "Adamts16", group.by = "celltype")

View(coembed@meta.data)
write.csv(coembed@meta.data, file = "coembed.meta.data.csv")










coembed$blacklist_region_fragments[is.na(coembed$blacklist_region_fragments)] <- 0
FeaturePlot(coembed, features = "blacklist_region_fragments", min.cutoff = 200, max.cutoff = 500)





fibro_center_atac <- AddMetaData(fibro_center_atac, metadata = celltype.predictions)

View(fibro_center_atac@meta.data)

hist(fibro_center_atac$prediction.score.max)
abline(v = 0.5, col = "red")
table(fibro_center_atac$prediction.score.max > 0.5)

fibro_center_atac.filtered <- subset(fibro_center_atac, subset = prediction.score.max > 0.5)
fibro_center_atac.filtered$predicted.id <- factor(fibro_center_atac.filtered$predicted.id, levels = levels(LW_center_rna))  # to make the colors match

p1 <- DimPlot(fibro_center_atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(LW_center_rna, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
  NoLegend()
CombinePlots(plots = list(p1, p2))





transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

##### Inquiring FeaturePlot from coembed data:

jpeg(file = "Crabp1_coembed.jpeg", width = 15, height = 15, units = "cm", res = 500)
FeaturePlot(coembed, features = "Crabp1", cols = c("yellow", "red"))
dev.off()

jpeg(file = "Rspo3_coembed.jpeg", width = 15, height = 15, units = "cm", res = 500)
FeaturePlot(coembed, features = "Rspo3", cols = c("yellow", "red"))
dev.off()

jpeg(file = "S100a4_coembed.jpeg", width = 15, height = 15, units = "cm", res = 500)
FeaturePlot(coembed, features = "S100a4", cols = c("yellow", "red"))
dev.off()

jpeg(file = "Cd200_coembed.jpeg", width = 15, height = 15, units = "cm", res = 500)
FeaturePlot(coembed, features = "Cd200", cols = c("yellow", "red"))
dev.off()



save(coembed, file = "coembed.Robj")


FeaturePlot(coembed, features = "chr16:92425483-92426458")

View(coembed@assays$ATAC@counts@Dimnames[[1]])

View(coembed@meta.data)

View(coembed@meta.data)


FeaturePlot(coembed, features = "chr1:3670740-3672492")
FeaturePlot(coembed, features = "Crabp1")

coembed_strip = coembed

coembed_strip@assays$RNA@scale.data = matrix(0,nrow = 1,ncol = 1)
coembed_strip@assays$RNA@var.features = matrix(0,nrow = 1,ncol = 1)
coembed_strip@meta.data = data.frame("CLEAR" = 1:1)
coembed_strip@assays$RNA@meta.features = data.frame("CLEAR" = 1:1)
coembed_strip@assays$RNA@counts = matrix(0,nrow = 1,ncol = 1)
coembed_strip@reductions$pca = matrix(0,nrow = 1,ncol = 1)
C = c("clear")
coembed_strip@active.ident = factor(C, levels= c("Clear"))
coembed_strip@assays$RNA@scale.data = matrix(0,nrow = 1,ncol = 1)
coembed_strip@assays$ATAC@meta.features = data.frame("CLEAR" = 1:1)
coembed_strip@assays$ATAC@scale.data = matrix(0,nrow = 1,ncol = 1)
coembed_strip@commands$ScaleData.RNA = matrix(0,nrow = 1,ncol = 1)
coembed_strip@commands$RunPCA.RNA = matrix(0,nrow = 1,ncol = 1)
coembed_strip@commands$RunUMAP.RNA.pca = matrix(0,nrow = 1,ncol = 1)
coembed_strip@assays$ATAC@counts = matrix(0,nrow = 1,ncol = 1)

coembed_strip@assays$RNA@x = data.frame("CLEAR" = 1:1)


FeaturePlot(coembed_strip, features = "chr1:3670740-3672492")
FeaturePlot(coembed_strip, features = "Crabp1")
  
save(coembed_strip, file = "coembed_strip.Robj")


load(file = "all_genes.Robj"); load(file = "all_chromatin.Robj")
combined_gene_chromatin = matrix(nrow = 114915, ncol = 1)
combined_gene_chromatin[1:20418] = all_genes
combined_gene_chromatin[20419:114915] = all_chromatin

View(combined_gene_chromatin)
save(combined_gene_chromatin, file = "combined_gene_chromatin.Robj")




#####################################################

coverageplot

coembed
