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

DefaultAssay(fibro_center_atac) <- "ACTIVITY"
fibro_center_atac <- FindVariableFeatures(fibro_center_atac)
fibro_center_atac <- NormalizeData(fibro_center_atac)
fibro_center_atac <- ScaleData(fibro_center_atac)


DefaultAssay(fibro_center_atac) <- "ATAC"
VariableFeatures(fibro_center_atac) <- names(which(Matrix::rowSums(fibro_center_atac) > 100))
fibro_center_atac <- RunLSI(fibro_center_atac, n = 50, scale.max = NULL)
fibro_center_atac <- RunUMAP(fibro_center_atac, reduction = "lsi", dims = 1:20)

p1 <- DimPlot(fibro_center_atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
p2 <- DimPlot(LW_center_rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))
