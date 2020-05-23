library(Seurat)
library(dplyr)

# Loading the Day8,14,18 dataset
wounds <- Read10X(data.dir = "/Users/SarthakSinha/Dropbox/Sarthak/Large_Wound_Alignments/Agg_Day8_Day14_Day18_Wounding/outs/filtered_gene_bc_matrices_mex/mm10/")



# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = wounds))
sparse.size <- object.size(x = wounds)
dense.size/sparse.size
# memory saving of 5.9 bytes when using sparse matrices

# My Seurat Object is called 'wound (pbmc)', and the raw data it takes is called 'wounds (pbmc.data)'
wound <- CreateSeuratObject(raw.data = wounds, min.cells = 3, min.genes = 400, 
                           project = "10X_Wound")

# "Low quality cells (<400 genes/cell and <3 cells/gene) were excluded from the overall experiment." - Shah Cell 18.

# QC and selecting cells for further analysis
mito.genes <- grep(pattern = "^mt-", x = rownames(x = wound@data), value = TRUE)
percent.mito <- Matrix::colSums(wound@raw.data[mito.genes, ])/Matrix::colSums(wound@raw.data)

wound <- AddMetaData(object = wound, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = wound, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# I have filtered out cells that have unique gene counts over 5,500 or less than 400, as done by Shah Cell 2018
wound <- FilterCells(object = wound, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(400, -Inf), high.thresholds = c(5500, 0.6))
wound <- FilterCells(object = wound, subset.names = c("nUMI"), 
                     low.thresholds = c(-Inf), high.thresholds = c(50000))

sample_1 <- grep("-1", colnames(wound@data), value = T); length(sample_1)
sample_2 <- grep("-2", colnames(wound@data), value = T); length(sample_2)
sample_3 <- grep("-3", colnames(wound@data), value = T); length(sample_3)
sample_4 <- grep("-4", colnames(wound@data), value = T); length(sample_4)
sample_5 <- grep("-5", colnames(wound@data), value = T); length(sample_5)
sample_6 <- grep("-6", colnames(wound@data), value = T); length(sample_6)
sample_7 <- grep("-7", colnames(wound@data), value = T); length(sample_7)
Day18_Day8_Lib1_try1 <- grep("-8", colnames(wound@data), value = T); length(Day18_Day8_Lib1_try1)
Day18_Day8_Lib2_try1 <- grep("-9", colnames(wound@data), value = T); length(Day18_Day8_Lib2_try1)
Day18_Day8_Lib3_try1 <- grep("-10", colnames(wound@data), value = T); length(Day18_Day8_Lib3_try1)
Day18_Day8_Lib4_try1 <- grep("-11", colnames(wound@data), value = T); length(Day18_Day8_Lib4_try1)

save(wound, file = "wound.Robj")
load(file = "wound.Robj")

####### QUALITY CHECK FOR WOUNDS ########
# Saving VlnPlot looking at "nGene", "nUMI", "percent.mito"
jpeg(file = "VlnPlot looking at nGene, nUMI, percent.mito", width = 30, height = 15, units = "cm", res = 500)
VlnPlot(object = wound_neg_0.3, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()
##

jpeg(file = "VlnPlot looking at nGene, nUMI, percent.mito", width = 30, height = 15, units = "cm", res = 500)
VlnPlot(object = wound_neg_0.3, features.plot = c("nGene"), nCol = 3)
dev.off()


jpeg(file = "VlnPlot looking a nUMI", width = 30, height = 15, units = "cm", res = 500)
VlnPlot(object = wound_neg_0.3, features.plot = c("nUMI"), nCol = 3)
dev.off()

jpeg(file = "VlnPlot looking a nUMI", width = 30, height = 15, units = "cm", res = 500)
VlnPlot(object = wound_neg_0.3, features.plot = c("percent.mito"), nCol = 3)
dev.off()


# Examining GenePlots
par(mfrow = c(1, 2))
GenePlot(object = wound, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = wound, gene1 = "nUMI", gene2 = "nGene")


# Saving GenePlot looking at "nGene", "nUMI", "percent.mito"
jpeg(file = "GenePlot looking at nGene, nUMI, percent.mito", width = 20, height = 15, units = "cm", res = 500)
GenePlot(object = wound, gene1 = "nUMI", gene2 = "percent.mito")
dev.off()
jpeg(file = "GenePlot looking at nGene, nUMI, percent.mito_2", width = 20, height = 15, units = "cm", res = 500)
GenePlot(object = wound, gene1 = "nUMI", gene2 = "nGene")
dev.off()

####### QUALITY CHECK FOR WOUNDS END ########


####### FILTERING OUT THE TD- cells ########
wound_neg <- FilterCells(wound, cells.use = c(sample_2, sample_4, sample_6, Day18_Day8_Lib4_try1), 
                 subset.names = "nGene", low.thresholds = 400, high.thresholds = 5500)







save(wound_neg, file = "wound_neg.Robj")

####### FILTERING OUT THE TD- cells ########



# Saving VlnPlots for nGene, nUMI, percent.mito AFTER Filteration - for tdTom-ve Samples
jpeg(file = "VlnPlot looking at nGene, nUMI, percent.mito - AFTER Filteration - neg samples", width = 30, height = 15, units = "cm", res = 500)
VlnPlot(object = wound_neg, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

# Normalization of Filtered Data Points
wound_neg <- NormalizeData(object = wound_neg, normalization.method = "LogNormalize", 
                      scale.factor = 10000)

# Detection of variable genes across the single cells
wound_neg <- FindVariableGenes(object = wound_neg, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = wound_neg@var.genes) #2494

# Scaling our data to remove unwanted sources of variation.
wound_neg <- ScaleData(object = wound_neg, vars.to.regress = c("nUMI", "percent.mito"))


# Perform linear dimensional reduction
wound_neg <- RunPCA(object = wound_neg, pc.genes = wound_neg@var.genes, do.print = TRUE, pcs.print = 1:5, 
               genes.print = 5)
PrintPCA(object = wound_neg, pcs.print = 1:5, genes.print = 5, use.full = FALSE)

####### PAUSE

saveRDS(wound, file = "/Users/SarthakSinha/Desktop/Seurat_Wounding/wound.rds")
saveRDS(wound_neg, file = "/Users/SarthakSinha/Desktop/Seurat_Wounding/wound_neg.rds")


### Resume

# Visualizing Heterogenity with Seurat tools to identify which PCs to include
VizPCA(object = wound_neg, pcs.use = 1:2)
PCAPlot(object = wound_neg, dim.1 = 1, dim.2 = 2)
wound_neg <- ProjectPCA(object = wound_neg, do.print = FALSE)


# Determine statistically significant principal components
wound_neg <- RunPCA(wound_neg, pc.genes = wound_neg@var.genes, pcs.compute = 200, do.print = T, pcs.print = 5, genes.print = 5)
PCHeatmap(wound_neg, pc.use = 20:40, cells.use = 100, do.balanced = TRUE, label.columns = F, use.full = F) ## To be saved
wound_neg <- JackStraw(wound_neg, num.pc = 40, num.replicate = 100, prop.freq = 0.10) 
JackStrawPlot(wound_neg_0.3, PCs = 1:40) ## To be saved
PCElbowPlot(wound_neg_0.3, num.pc = 50) ## To be saved
save(wound_neg, file = "wound_neg.Robj")

wound_neg_1 <- RunPCA(wound_neg_1, pc.genes = wound_neg_1@var.genes, pcs.compute = 50)
 ## Save this plot

#################Cell clustering - PCs 1:20#################################
## PC's 1 to 28 were used based on PC scoring metrics calculated above using PCElbowPlot analysis
wound_neg <- FindClusters(wound_neg, reduction.type = "pca", dims.use = 1:28, resolution = 0.1, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
wound_neg <- StashIdent(wound_neg, save.name = "ClusterNames_0.1")
wound_neg <- RunTSNE(wound_neg, dims.use = 1:28, do.fast = T)
save(wound_neg, file = "wound_neg.Robj")

wound_neg_1 <- FindClusters(wound_neg, reduction.type = "pca", dims.use = 1:28, resolution = 0.2, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
wound_neg_0.3 <- FindClusters(wound_neg, reduction.type = "pca", dims.use = 1:28, resolution = 0.3, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

jpeg(file = "TSNEPlot of Neg Cells", width = 20, height = 15, units = "cm", res = 500)
TSNEPlot(wound_neg_0.3, do.label = TRUE, pt.size = 0.8)
dev.off()

save(wound_neg_0.3, file = "wound_neg_0.3.Robj")
load(file = "wound_neg_0.3.Robj")

### Pericyte markers
pericyte.markers <- FindMarkers(wound_neg_0.3, ident.1 = "Pericytes", min.pct = 0.25, test.use = "negbinom")
pericyte.markers <- pericyte.markers[order(pericyte.markers$avg_logFC),]
write.csv(pericyte.markers,file = "pericyte.markers.markers")


pdf("pericyte_genes.pdf")
TSNEPlot(wound_neg_0.3, do.label = TRUE, pt.size = 0.8)
FeaturePlot(wound_neg_0.3, features.plot = c("Sox18"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Acta2"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Tagln"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Myl9"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Ccl2"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Sparcl1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Mylk"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Mustn1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Cxcl1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Tpm2"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Mgp"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Il6"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Igfbp5"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Nrg1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()


pdf('Inductive Fibro in tdTom-ve Cell Population.pdf')
TSNEPlot(wound_neg_0.3, do.label = T, pt.size = 0.5)
FeaturePlot(wound_neg_0.3, features.plot = c("Pdgfra"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Dpt"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Sox18"), reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()

pdf('Inductive Fibro in tdTom-ve Cell Population.pdf')
TSNEPlot(wound_neg_0.3, do.label = T, pt.size = 0.5)
FeaturePlot(wound_neg_0.3, features.plot = c("Pdgfra"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Dpt"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Sox18"), reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()



### CURRENT Endothelial markers
endothelial.markers <- FindMarkers(wound_neg_0.3, ident.1 = "Endothelial cells", min.pct = 0.25, test.use = "negbinom")
endothelial.markers <- endothelial.markers[order(endothelial.markers$avg_logFC),]
write.csv(endothelial.markers,file = "endothelial.markers.csv")
#
pdf("endo_genes.pdf")
TSNEPlot(wound_neg_0.3, do.label = TRUE, pt.size = 0.8)
FeaturePlot(wound_neg_0.3, features.plot = c("Fabp4"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Tm4sf1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Gng11"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Ctla2a"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Cav1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Cldn5"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Egfl7"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Pecam1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Cdh5"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Ecscr"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Sparcl1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Cd93"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Rnd1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Myct1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()



### Isolate Pericytes to look at pericytes seperated by samples:

## Filtering Pericytes
pericytes_cells <- names(wound_neg_0.3@ident[wound_neg_0.3@ident %in% c("Pericytes")])
sub.pericytes <- FilterCells(wound_neg_0.3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = pericytes_cells)
TSNEPlot(sub.pericytes, pt.size = 0.5)


### Isolating Epithelial cells:






## Subdividing Macrophases into more clusters
sub.pericytes_0.6 <- FindClusters(sub.pericytes, reduction.type = "pca", dims.use = 1:25, resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
TSNEPlot(sub.pericytes_0.6, pt.size = 0.5)




## subsetting data by Library ID's
sample_2 <- grep("-2", colnames(wound_neg_1@data), value = T); length(sample_2)
sample_4 <- grep("-4", colnames(wound_neg_1@data), value = T); length(sample_4)
sample_6 <- grep("-6", colnames(wound_neg_1@data), value = T); length(sample_6)
Day18_Day8_Lib4_try1 <- grep("-11", colnames(wound_neg_1@data), value = T); length(Day18_Day8_Lib4_try1)


## graphing tSNE plot by sample number
s_tSNE <- as.data.frame(GetCellEmbeddings(wound_neg_1, reduction.type = "tsne"))
s_tSNE[Day18_Day8_Lib4_try1, "Sample"] <- "Day 8 tdT- SW"
s_tSNE[sample_6, "Sample"] <- "Day 14 tdT- SW"
s_tSNE[sample_4, "Sample"] <- "Day 14 tdT- LWO"
s_tSNE[sample_2, "Sample"] <- "Day 14 tdT- LWC"

ggplot(s_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()


########## Schwann Cells ######
Schwann_genes = c("Gpm6b", "Sox10", "Plp1", "S100b")
jpeg("Schwann Cells", width = 20, height = 20, units = "cm", res = 500)
FeaturePlot(wound_neg_0.3, features.plot = Schwann_genes, reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()

FeaturePlot(wound_neg_0.3, features.plot = c("Scx"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(wound_neg_0.3, features.plot = c("Fabp5"), reduction.use = "tsne", cols.use = c("grey", "blue"))


## 9 is Schwann Cells

########## Macrophases ######

Macrophages_genes = c("Cd68", "Adgre1", "Fcgr1", "Itgam")
jpeg("Macrophages", width = 20, height = 20, units = "cm", res = 500)
FeaturePlot(wound_neg_0.3, features.plot = Macrophages_genes, reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()

## 3 and 5 are Macrophases

########## T Cells ######
Tcell_genes = c("Cd4", "Cxcr6", "Ctla2a", "Cd3g")
jpeg("T cells", width = 20, height = 20, units = "cm", res = 500)
FeaturePlot(wound_neg_0.3, features.plot = Tcell_genes, reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()

## 8 is T Cells

########## Endothelial Cells ######

Endothelial_genes = c("Pecam1","Egfl7", "Ecscr", "Cdh5")
jpeg("Endothelial_cells", width = 20, height = 20, units = "cm", res = 500)
FeaturePlot(wound_neg_0.3, features.plot = Endothelial_genes, reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()

## 16 is Endothelial cells


########## Neutrophils ######

Neutrophil_genes = c("Cd177","Lcn2","Ly6g","Ngp")
jpeg("Neutrophils cell", width = 20, height = 20, units = "cm", res = 500)
FeaturePlot(wound_neg_0.3, features.plot = Neutrophil_genes, reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()

## 18 is Neutrophils



########## Langerhans ######

Langerhans_genes = c("Cd207","H2-M5","Cd207","H2-M5")
jpeg("Langerhans", width = 20, height = 20, units = "cm", res = 500)
FeaturePlot(wound_neg_0.3, features.plot = Langerhans_genes, reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()

## 6 is Langerhans


########## Pericytes ######

Pericytes_genes = c("Pdgfrb","Notch3","Des","Acta2")
jpeg("Pericytes", width = 20, height = 20, units = "cm", res = 500)
FeaturePlot(wound_neg_0.3, features.plot = Pericytes_genes, reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()

## 14 is Pericytes

########## Skeletal Muscle Progenitors ######

Skeletal_genes = c("Chrna1","Cdh15","Pax7","Myf5")
jpeg("Skeletal Muscle Progenitors", width = 20, height = 20, units = "cm", res = 500)
FeaturePlot(wound_neg_0.3, features.plot = Skeletal_genes, reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()

FeaturePlot(wound_neg_0.3, features.plot = c("Prg4"), reduction.use = "tsne", cols.use = c("grey", "blue"))


## 17 is Skeletal Muscle Progenitors
FeaturePlot(wound_neg_0.3, features.plot = c(""), reduction.use = "tsne", cols.use = c("grey", "blue"))


########## Dendritic cells ######

Dendritic_genes = c("Siglech","Cox6a2","Siglech","Cox6a2")
jpeg("Dendritic_cells", width = 20, height = 20, units = "cm", res = 500)
FeaturePlot(wound_neg_0.3, features.plot =Dendritic_genes, reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()

## 20 is Dendritic cells


########## Melanocytes ######

Melanocytes_genes = c("Mitf","Dct","Mitf","Dct")
jpeg("Melanocytes_genes", width = 20, height = 20, units = "cm", res = 500)
FeaturePlot(wound_neg_0.3, features.plot =Melanocytes_genes, reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()

###

########## tdTom-ve Fibroblasts ######

Fibro_genes = c("Col6a2","Mmp2","Col5a1","Pdgfra")
jpeg("Fibro", width = 20, height = 20, units = "cm", res = 500)
FeaturePlot(wound_neg_0.3, features.plot =Fibro_genes, reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()

## 10 and 7 are tdTom-ve Fibroblasts

########## keratinocytes ######
keratinocytes_genes = c("Cdh1","Krt14","Krt5","Krt28")
jpeg("keratinocytes", width = 20, height = 20, units = "cm", res = 500)
FeaturePlot(wound_neg_0.3, features.plot =keratinocytes_genes, reduction.use = "tsne", cols.use = c("grey", "blue"))
dev.off()
## 11, 0, 2, 4, 12, 1, 13 Keratinocytes

Epi13.markers <- FindMarkers(wound_neg_0.3, ident.1 = 13, min.pct = 0.25, test.use = "negbinom")
print(x = head(x = Epi13.markers, n = 100))


# 15??? #
# 19??? #
# 20??? #

Cluster15.markers <- FindMarkers(wound_neg_0.3, ident.1 = 15, min.pct = 0.25, test.use = "negbinom")
Cluster19.markers <- FindMarkers(wound_neg_0.3, ident.1 = 19, min.pct = 0.25, test.use = "negbinom")
Cluster20.markers <- FindMarkers(wound_neg_0.3, ident.1 = 19, min.pct = 0.25, test.use = "negbinom")

print(x = head(x = Cluster15.markers, n = 100))
Cluster15.markers <- Cluster15.markers[order(Cluster15.markers$avg_logFC),]
write.csv(Cluster15.markers,file = "Cluster15.markers.csv")

print(x = head(x = Cluster19.markers, n = 100))
Cluster19.markers <- Cluster19.markers[order(Cluster19.markers$avg_logFC),]
write.csv(Cluster19.markers,file = "Cluster19.markers.csv")

print(x = head(x = Cluster20.markers, n = 100))
Cluster20.markers <- Cluster20.markers[order(Cluster20.markers$avg_logFC),]
write.csv(Cluster20.markers,file = "Cluster20.markers.csv")


### Differentiating Melanocytes ###
Melan_diff.markers <- FindMarkers(wound_neg_0.3, ident.1 = 19, ident.2 = 20, min.pct = 0.25, test.use = "negbinom")
Melan_diff.markers <- Melan_diff.markers[order(Melan_diff.markers$avg_logFC),]
write.csv(Melan_diff.markers,file = "Melan_diff.markers")

### Differentiating Macrophases ###
Macro_diff.markers <- FindMarkers(wound_neg_0.3, ident.1 = 3, ident.2 = 5, min.pct = 0.25, test.use = "negbinom")
Macro_diff.markers <- Macro_diff.markers[order(Macro_diff.markers$avg_logFC),]
write.csv(Macro_diff.markers,file = "Macro_diff.markers")



FeaturePlot(wound_neg_0.3, features.plot =c("Tfap2a"), reduction.use = "tsne", cols.use = c("grey", "blue"))

### Visualizing Samples on the same tSNE:
wound_neg_0.3_tSNE <- as.data.frame(GetCellEmbeddings(wound_neg_0.3, reduction.type = "tsne", dims.use = 1:2))
wound_neg_0.3_tSNE[Day18_Day8_Lib4_try1, "Sample"] <- "Day 8 tdT- SW"
wound_neg_0.3_tSNE[sample_6, "Sample"] <- "Day 14 tdT- SW"
wound_neg_0.3_tSNE[sample_4, "Sample"] <- "Day 14 tdT- LWO"
wound_neg_0.3_tSNE[sample_2, "Sample"] <- "Day 14 tdT- LWC"

ggplot(wound_neg_0.3_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()


########## Standerdized Distribution Visualization ######

### Large wound Center Cell Types
lwc <- FilterCells(wound_neg_0.3, cells.use = c(sample_2), 
                   subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
lwc <- FindVariableGenes(lwc, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
lwc <- NormalizeData(object = lwc, normalization.method = "LogNormalize", 
                     scale.factor = 10000)
lwc <- FindVariableGenes(object = lwc, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
lwc <- ScaleData(object = lwc, vars.to.regress = c("nUMI", "percent.mito"))
lwc <- RunPCA(lwc, pc.genes = q@var.genes, pcs.compute = 50)
PCElbowPlot(lwc, num.pc = 50)
lwc <- FindClusters(lwc, reduction.type = "pca", dims.use = 1:28, resolution = 1.5, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
lwc <- StashIdent(lwc, save.name = "ClusterNames_0.1_lwc")
lwc <- RunTSNE(lwc, dims.use = 1:28, do.fast = T)
#TSNEPlot(lwc)

### Large wound Outside Cell Types
lwo <- FilterCells(wound_neg_0.3, cells.use = c(sample_4), 
                   subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
lwo <- FindVariableGenes(lwo, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)

lwo <- NormalizeData(object = lwo, normalization.method = "LogNormalize", 
                     scale.factor = 10000)
lwo <- FindVariableGenes(object = lwo, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
lwo <- ScaleData(object = lwo, vars.to.regress = c("nUMI", "percent.mito"))
lwo <- RunPCA(lwo, pc.genes = lwo@var.genes, pcs.compute = 50)
PCElbowPlot(lwo, num.pc = 50)

lwo <- FindClusters(lwo, reduction.type = "pca", dims.use = 1:28, resolution = 1.5, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
lwo <- StashIdent(lwo, save.name = "ClusterNames_1.5_lwo")
lwo <- RunTSNE(lwo, dims.use = 1:28, do.fast = T)
TSNEPlot(object = lwo, do.label = TRUE, pt.size = 0.5)


## For JOVE
lwo_jove <- FindClusters(lwo, reduction.type = "pca", dims.use = 1:20, resolution = 4, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
lwo_jove <- StashIdent(lwo_jove, save.name = "ClusterNames_1.5_lwo")
lwo_jove <- RunTSNE(lwo_jove, dims.use = 1:20, do.fast = T)
TSNEPlot(object = lwo_jove, do.label = TRUE, pt.size = 0.5)
FeaturePlot(lwo_jove, features.plot = c("Sox10"), reduction.use = "tsne", cols.use = c("grey", "blue"))

FeaturePlot(lwo_jove, features.plot = c("Sox10"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(lwo_jove, features.plot = c("Cd68"), no.legend = FALSE, 
            min.cutoff = "q10", max.cutoff = "q90", dark.theme = TRUE)

### Small wound Cell Types
sw <- FilterCells(wound_neg_0.3, cells.use = c(sample_6), 
                  subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
sw <- FindVariableGenes(sw, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
sw <- NormalizeData(object = sw, normalization.method = "LogNormalize", 
                    scale.factor = 10000)
sw <- FindVariableGenes(object = sw, mean.function = ExpMean, dispersion.function = LogVMR, 
                        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
sw <- ScaleData(object = sw, vars.to.regress = c("nUMI", "percent.mito"))
sw <- RunPCA(sw, pc.genes = sw@var.genes, pcs.compute = 50)
PCElbowPlot(sw, num.pc = 50)
sw <- FindClusters(sw, reduction.type = "pca", dims.use = 1:28, resolution = 1.5, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
sw <- StashIdent(sw, save.name = "ClusterNames_1.5_lwo")
sw <- RunTSNE(sw, dims.use = 1:28, do.fast = T)
#TSNEPlot(object = sw, do.label = TRUE, pt.size = 0.5)

### Small wound Day8 Cell Types
sw8 <- FilterCells(wound_neg_0.3, cells.use = c(Day18_Day8_Lib4_try1), 
                   subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
sw8 <- FindVariableGenes(sw8, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
sw8 <- NormalizeData(object = sw8, normalization.method = "LogNormalize", 
                     scale.factor = 10000)
sw8 <- FindVariableGenes(object = sw8, mean.function = ExpMean, dispersion.function = LogVMR, 
                         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
sw8 <- ScaleData(sw8, vars.to.regress = c("nUMI", "percent.mito"))
sw8 <- RunPCA(sw8, pc.genes = sw8@var.genes, pcs.compute = 50)
PCElbowPlot(sw8, num.pc = 50)
sw8 <- FindClusters(sw8, reduction.type = "pca", dims.use = 1:28, resolution = 1.5, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
sw8 <- StashIdent(sw8, save.name = "ClusterNames_1.5_lwo")
sw8 <- RunTSNE(sw8, dims.use = 1:28, do.fast = T)
#TSNEPlot(object = sw8, do.label = TRUE, pt.size = 0.5)



### Choosing cells randomly from Samples 2, 4, 6, Lib 4 to show distribution with even number of cells
random_lwc <- sample(1:2683, 1700, replace = FALSE)
random_lwo <- sample(1:4565, 1700, replace = FALSE)
random_sw <- sample(1:3476, 1700, replace = FALSE)
random_sw8 <- sample(1:1787, 1700, replace = FALSE) 

TSNEPlot(wound_neg_0.3)


########## Standerdized tSNEs for Figure 1  ######
jpeg(file = "LWC", width = 18, height = 15, units = "cm", res = 500)
TSNEPlot(wound_neg_0.3, cells.use = lwc@cell.names[random_lwc])
dev.off()

jpeg(file = "LWO", width = 18, height = 15, units = "cm", res = 500)
TSNEPlot(wound_neg_0.3, cells.use = lwo@cell.names[random_lwo])
dev.off()

jpeg(file = "SW14", width = 18, height = 15, units = "cm", res = 500)
TSNEPlot(wound_neg_0.3, cells.use = sw@cell.names[random_sw])
dev.off()

jpeg(file = "SW8", width = 18, height = 15, units = "cm", res = 500)
TSNEPlot(wound_neg_0.3, cells.use = sw8@cell.names[random_sw8])
dev.off()

########## Labeling Clusters - changing ident ######
current <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
new <- c("Epidermal KCS",
         "Epidermal KCS", 
         "Epidermal KCS", 
         "M1 Macrophages", 
         "Epidermal KCS", 
         "M2 Macrophages", 
         "Langerhans", 
         "Fibroblasts", 
         "T Cells", 
         "Schwann Cells",
         "Fibroblasts",
         "Epidermal KCS",
         "Epidermal KCS",
         "Epidermal KCS",
         "Pericytes", 
         "WAT macrophases/TAM", 
         "Endothelial cells", 
         "Skeletal Muscle Progenitors", 
         "Neutrophils", 
         "Malignant Melanocytes", 
         "Melanocytes")
wound_neg_0.3@ident <- plyr::mapvalues(wound_neg_0.3@ident, from = current, to = new)
TSNEPlot(wound_neg_0.3, do.label = TRUE, pt.size = 0.5)


########## Re-Labeling Clusters - changing ident ######
current <- c("Epidermal KCS",
         "Epidermal KCS", 
         "Epidermal KCS", 
         "Macrophages", 
         "Epidermal KCS", 
         "Macrophages", 
         "Langerhans", 
         "Fibroblasts", 
         "T Cells", 
         "Schwann Cells",
         "Fibroblasts",
         "Epidermal KCS",
         "Epidermal KCS",
         "Epidermal KCS",
         "Pericytes", 
         "Macrophages", 
         "Endothelial cells", 
         "Skeletal Muscle Progenitors", 
         "Neutrophils", 
         "Melanoblast", 
         "Melanocytes")
new <- c("Epidermal KCS",
             "Epidermal KCS", 
             "Epidermal KCS", 
             "Macrophages", 
             "Epidermal KCS", 
             "Macrophages", 
             "Langerhans", 
             "Fibroblasts", 
             "T Cells", 
             "Schwann Cells",
             "Fibroblasts",
             "Epidermal KCS",
             "Epidermal KCS",
             "Epidermal KCS",
             "Pericytes", 
             "Macrophages", 
             "Endothelial cells", 
             "Skeletal Muscle Progenitors", 
             "Neutrophils", 
             "Melanocytes", 
             "Dendritic Cells")

wound_neg_0.3@ident <- plyr::mapvalues(wound_neg_0.3@ident, from = current, to = new)
TSNEPlot(wound_neg_0.3, do.label = TRUE, pt.size = 0.5)






TSNEPlot(wound_neg_0.3, do.label = TRUE, pt.size = 0.5)
TSNEPlot(wound_neg_0.3, do.label = TRUE, cells.use = lwc@cell.names[random_lwc])
TSNEPlot(wound_neg_0.3, do.label = TRUE, cells.use = lwo@cell.names[random_lwo])
TSNEPlot(wound_neg_0.3, do.label = TRUE, cells.use = sw@cell.names[random_sw])
TSNEPlot(wound_neg_0.3, do.label = TRUE, cells.use = sw8@cell.names[random_sw8])

save(wound_neg_0.3, file = "wound_neg_0.3.Robj")
save(sw8, file = "sw8.Robj")
save(sw, file = "sw.Robj")
save(lwo, file = "lwo.Robj")
save(lwc, file = "lwc.Robj")


load(file = "sw8.Robj")
load(file = "sw.Robj")
load(file = "lwo.Robj")
load(file = "lwc.Robj")
load(file = "wound_neg_0.3.Robj")



#### For Figure ###

TSNEPlot(wound_neg_0.3, do.label = TRUE, pt.size = 0.5)
TSNEPlot(wound_neg_0.3, do.label = TRUE, cells.use = lwc@cell.names[random_lwc])
TSNEPlot(wound_neg_0.3, do.label = TRUE, cells.use = lwo@cell.names[random_lwo])
TSNEPlot(wound_neg_0.3, do.label = TRUE, cells.use = sw@cell.names[random_sw])
TSNEPlot(wound_neg_0.3, do.label = TRUE, cells.use = sw8@cell.names[random_sw8])


jpeg(file = "tSNE", width = 25, height = 20, units = "cm", res = 500)
TSNEPlot(wound_neg_0.3, pt.size = 0.5)
dev.off()


########## Heatmap for top 10 markers ######

wound_neg_0.3.markers <- FindAllMarkers(object = wound_neg_0.3, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
wound_neg_0.3.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
#
top10 <- wound_neg_0.3.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
top5 <- wound_neg_0.3.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
top3 <- wound_neg_0.3.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)
DoHeatmap(object = wound_neg_0.3, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

jpeg(file = "tSNE", width = 25, height = 20, units = "cm", res = 500)
DoHeatmap(object = wound_neg_0.3, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()



########## Supplementary Figure 1 Panels ######
jpeg(file = "Heatmap", width = 25, height = 20, units = "cm", res = 500)
DoHeatmap(object = wound_neg_0.3, genes.use = top3$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()

jpeg(file = "ElbowPlot", width = 25, height = 20, units = "cm", res = 500)
PCElbowPlot(wound_neg_0.3, num.pc = 50) ## To be saved
dev.off()

jpeg(file = "HeatMap", width = 25, height = 20, units = "cm", res = 500)
PCHeatmap(wound_neg_0.3, pc.use = 20:30, cells.use = 100, do.balanced = TRUE, label.columns = F, use.full = F) ## To be saved
dev.off()

jpeg(file = "GenePlot_1", width = 25, height = 20, units = "cm", res = 500)
GenePlot(wound_neg_0.3, gene1 = "nUMI", gene2 = "percent.mito", col.use = "Black")
dev.off()
jpeg(file = "GenePlot_2", width = 25, height = 20, units = "cm", res = 500)
GenePlot(wound_neg_0.3, gene1 = "nUMI", gene2 = "nGene", col.use = "Black")
dev.off()

# Neutrophil Markers for IPA
Neutro.markers <- FindMarkers(wound_neg_0.3, ident.1 = c("Neutrophils"), min.pct = 0.25, test.use = "negbinom")
Neutro.markers <- Neutro.markers[order(Neutro.markers$avg_logFC),]
write.csv(Neutro.markers,file = "Neutro.markers.csv")
summary(wound_neg_0.3@ident)
summary(wound_neg_0.3$sample1)

# Langerhand Markers for Wisoo
Langer.markers <- FindMarkers(wound_neg_0.3, ident.1 = c("Langerhans"), min.pct = 0.25, test.use = "negbinom")
Langer.markers <- Langer.markers[order(Langer.markers$avg_logFC),]
write.csv(Langer.markers,file = "Langerhan_markers.csv")
summary(wound_neg_0.3@ident)




## Filtering Macrophases
macrophase_cells <- names(wound_neg_0.3@ident[wound_neg_0.3@ident %in% c("Macrophages")])
sub_sub.macro <- FilterCells(wound_neg_0.3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = macrophase_cells)
TSNEPlot(sub_sub.macro, pt.size = 0.5)
## Subdividing Macrophases into more clusters
sub_sub.macro_0.6 <- FindClusters(sub_sub.macro, reduction.type = "pca", dims.use = 1:28, resolution = 0.6, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)
TSNEPlot(sub_sub.macro_0.6, pt.size = 0.5)
### Choosing macrophases randomly from Samples 2, 4, 6, Lib 4 to show standerdized distribution
random_macro_lwc <- sample(1:775, 170, replace = FALSE)
random_macro_lwo <- sample(1:683, 170, replace = FALSE)
random_macro_sw <- sample(1:175, 170, replace = FALSE)
random_macro_sw8 <- sample(1:237, 170, replace = FALSE)
sample_2 <- grep("-2", colnames(sub_sub.macro_0.6@data), value = T); length(sample_2)
sample_4 <- grep("-4", colnames(sub_sub.macro_0.6@data), value = T); length(sample_4)
sample_6 <- grep("-6", colnames(sub_sub.macro_0.6@data), value = T); length(sample_6)
sample_11 <- grep("-11", colnames(sub_sub.macro_0.6@data), value = T); length(sample_11)

TSNEPlot(sub_sub.macro_0.6, cells.use = sample_2[random_macro_lwc])
TSNEPlot(sub_sub.macro_0.6, cells.use = sample_4[random_macro_lwo])
TSNEPlot(sub_sub.macro_0.6, cells.use = sample_6[random_macro_sw])
TSNEPlot(sub_sub.macro_0.6, cells.use = sample_11[random_macro_sw8])

standerdized_macro <- c(sample_2[random_macro_lwc], sample_4[random_macro_lwo], sample_6[random_macro_sw], sample_11[random_macro_sw8])
standerdized_macro_filter <- FilterCells(sub_sub.macro_0.6, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = standerdized_macro)

TSNEPlot(standerdized_macro_filter)

macro_tSNE_standerized <- as.data.frame(GetCellEmbeddings(standerdized_macro_filter, reduction.type = "tsne"))
macro_tSNE_standerized[sample_2, "Sample"] <- "tdT- Large Wound Center"
macro_tSNE_standerized[sample_4, "Sample"] <- "tdT- Large Wound Outside"
macro_tSNE_standerized[sample_6, "Sample"] <- "tdT- Small Wound 14"
macro_tSNE_standerized[sample_11, "Sample"] <- "tdT- Small Wound 8"
### STANDERDIZED Macrophase Graph ## SAVE!
jpeg(file = "tsne by sample #.jpeg", width = 30, height = 20, units = "cm", res = 500)
ggplot(macro_tSNE_standerized, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()


VlnPlot(sub_sub.macro_0.6, features.plot = c("Serpinb6b"))


########## Standerdized tSNEs for Figure 1  ######
jpeg(file = "LWC", width = 18, height = 15, units = "cm", res = 500)
TSNEPlot(sub_sub.macro_0.6, cells.use = lwc@cell.names[random_macro_lwc])
dev.off()


## Filtering EPIDERMAL KERATINOCYTES
epidermal_cells <- names(wound_neg_0.3@ident[wound_neg_0.3@ident %in% c("Epidermal KCS")])
sub_sub.epi <- FilterCells(wound_neg_0.3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = epidermal_cells)
TSNEPlot(sub_sub.epi, pt.size = 0.5)

## Filtering Fibroblasts
fibroblasts_cells <- names(wound_neg_0.3@ident[wound_neg_0.3@ident %in% c("Fibroblasts")])
sub_sub.fibroblasts <- FilterCells(wound_neg_0.3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = fibroblasts_cells)
TSNEPlot(sub_sub.fibroblasts, pt.size = 0.5)
## Filtering Pericytes
pericytes_cells <- names(wound_neg_0.3@ident[wound_neg_0.3@ident %in% c("Pericytes")])
sub_sub.pericytes <- FilterCells(wound_neg_0.3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = pericytes_cells)
TSNEPlot(sub_sub.pericytes, pt.size = 0.5)
## Filtering Langerhans
langerhans_cells <- names(wound_neg_0.3@ident[wound_neg_0.3@ident %in% c("Langerhans")])
sub_sub.langerhans <- FilterCells(wound_neg_0.3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = langerhans_cells)
TSNEPlot(sub_sub.langerhans, pt.size = 0.5)
## Filtering T Cells
T_cells <- names(wound_neg_0.3@ident[wound_neg_0.3@ident %in% c("T Cells")])
sub_sub.Tcells <- FilterCells(wound_neg_0.3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = T_cells)
TSNEPlot(sub_sub.Tcells, pt.size = 0.5)
## Filtering Schwann Cells
Schwann_cells <- names(wound_neg_0.3@ident[wound_neg_0.3@ident %in% c("Schwann Cells")])
sub_sub.Schwann_cells <- FilterCells(wound_neg_0.3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = Schwann_cells)
TSNEPlot(sub_sub.Schwann_cells, pt.size = 0.5)
## Filtering Endothelial cells
Endothelial_cells <- names(wound_neg_0.3@ident[wound_neg_0.3@ident %in% c("Endothelial cells")])
sub_sub.Endothelial_cells <- FilterCells(wound_neg_0.3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = Endothelial_cells)
TSNEPlot(sub_sub.Endothelial_cells, pt.size = 0.5)
## Filtering Skeletal Muscle Progenitors
skeletal_cells <- names(wound_neg_0.3@ident[wound_neg_0.3@ident %in% c("Skeletal Muscle Progenitors")])
sub_sub.skeletal <- FilterCells(wound_neg_0.3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = skeletal_cells)
TSNEPlot(sub_sub.skeletal, pt.size = 0.5)
## Filtering Neutrophils
Neutrophils_cells <- names(wound_neg_0.3@ident[wound_neg_0.3@ident %in% c("Neutrophils")])
sub_sub.Neutrophils <- FilterCells(wound_neg_0.3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = Neutrophils_cells)
TSNEPlot(sub_sub.Neutrophils, pt.size = 0.5)
## Filtering Melanocytes
Melanocytes_cells <- names(wound_neg_0.3@ident[wound_neg_0.3@ident %in% c("Melanocytes")])
sub_sub.Melanocytes <- FilterCells(wound_neg_0.3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = Melanocytes_cells)
TSNEPlot(sub_sub.Melanocytes, pt.size = 0.5)
## Filtering Dendritic Cells
Dendritic_cells <- names(wound_neg_0.3@ident[wound_neg_0.3@ident %in% c("Dendritic Cells")])
sub_sub.Dendritic <- FilterCells(wound_neg_0.3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = Dendritic_cells)
TSNEPlot(sub_sub.Dendritic, pt.size = 0.5)


# Assigning sample #s:
sample_2 <- grep("-2", colnames(wound_neg_0.3@data), value = T); length(sample_2)
sample_4 <- grep("-4", colnames(wound_neg_0.3@data), value = T); length(sample_4)
sample_6 <- grep("-6", colnames(wound_neg_0.3@data), value = T); length(sample_6)
sample_11 <- grep("-11", colnames(wound_neg_0.3@data), value = T); length(sample_11)

## LOOKING AT tSNEs BY SAMPLES - Global, pooled ##
s_tSNE <- as.data.frame(GetCellEmbeddings(wound_neg_0.3, reduction.type = "tsne"))
s_tSNE[sample_2, "Sample"] <- "tdT- Large Wound Center"
s_tSNE[sample_4, "Sample"] <- "tdT- Large Wound Outside"
s_tSNE[sample_6, "Sample"] <- "tdT- Small Wound 14"
s_tSNE[sample_11, "Sample"] <- "tdT- Small Wound 8"

## LOOKING AT tSNEs BY SAMPLES - Macrophases, pooled ##
macro_tSNE_sample <- as.data.frame(GetCellEmbeddings(sub_sub.macro, reduction.type = "tsne"))
macro_tSNE_sample[sample_2, "Sample"] <- "tdT- Large Wound Center"
macro_tSNE_sample[sample_4, "Sample"] <- "tdT- Large Wound Outside"
macro_tSNE_sample[sample_6, "Sample"] <- "tdT- Small Wound 14"
macro_tSNE_sample[sample_11, "Sample"] <- "tdT- Small Wound 8"

jpeg(file = "tsne by sample #.jpeg", width = 30, height = 20, units = "cm", res = 500)
ggplot(macro_tSNE_sample, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()


## LOOKING AT tSNEs BY SAMPLES - Epidermal Keratinocytes, pooled ##
epi_tSNE_sample <- as.data.frame(GetCellEmbeddings(sub_sub.epi, reduction.type = "tsne"))
epi_tSNE_sample[sample_2, "Sample"] <- "tdT- Large Wound Center"
epi_tSNE_sample[sample_4, "Sample"] <- "tdT- Large Wound Outside"
epi_tSNE_sample[sample_6, "Sample"] <- "tdT- Small Wound 14"
epi_tSNE_sample[sample_11, "Sample"] <- "tdT- Small Wound 8"

jpeg(file = "tsne by sample #.jpeg", width = 30, height = 20, units = "cm", res = 500)
ggplot(epi_tSNE_sample, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()

## LOOKING AT tSNEs BY SAMPLES - Fibroblasts, pooled ##
fib_tSNE_sample <- as.data.frame(GetCellEmbeddings(sub_sub.fibroblasts, reduction.type = "tsne"))
fib_tSNE_sample[sample_2, "Sample"] <- "tdT- Large Wound Center"
fib_tSNE_sample[sample_4, "Sample"] <- "tdT- Large Wound Outside"
fib_tSNE_sample[sample_6, "Sample"] <- "tdT- Small Wound 14"
fib_tSNE_sample[sample_11, "Sample"] <- "tdT- Small Wound 8"

jpeg(file = "tsne by sample #.jpeg", width = 30, height = 20, units = "cm", res = 500)
ggplot(fib_tSNE_sample, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()

## LOOKING AT tSNEs BY SAMPLES - Pericytes, pooled ##
per_tSNE_sample <- as.data.frame(GetCellEmbeddings(sub_sub.pericytes, reduction.type = "tsne"))
per_tSNE_sample[sample_2, "Sample"] <- "tdT- Large Wound Center"
per_tSNE_sample[sample_4, "Sample"] <- "tdT- Large Wound Outside"
per_tSNE_sample[sample_6, "Sample"] <- "tdT- Small Wound 14"
per_tSNE_sample[sample_11, "Sample"] <- "tdT- Small Wound 8"

jpeg(file = "tsne by sample #.jpeg", width = 30, height = 20, units = "cm", res = 500)
ggplot(per_tSNE_sample, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()


## LOOKING AT tSNEs BY SAMPLES - Langerhans, pooled ##
lan_tSNE_sample <- as.data.frame(GetCellEmbeddings(sub_sub.langerhans, reduction.type = "tsne"))
lan_tSNE_sample[sample_2, "Sample"] <- "tdT- Large Wound Center"
lan_tSNE_sample[sample_4, "Sample"] <- "tdT- Large Wound Outside"
lan_tSNE_sample[sample_6, "Sample"] <- "tdT- Small Wound 14"
lan_tSNE_sample[sample_11, "Sample"] <- "tdT- Small Wound 8"

jpeg(file = "tsne by sample #.jpeg", width = 30, height = 20, units = "cm", res = 500)
ggplot(lan_tSNE_sample, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()

## LOOKING AT tSNEs BY SAMPLES - T Cells, pooled ##
T_tSNE_sample <- as.data.frame(GetCellEmbeddings(sub_sub.Tcells, reduction.type = "tsne"))
T_tSNE_sample[sample_2, "Sample"] <- "tdT- Large Wound Center"
T_tSNE_sample[sample_4, "Sample"] <- "tdT- Large Wound Outside"
T_tSNE_sample[sample_6, "Sample"] <- "tdT- Small Wound 14"
T_tSNE_sample[sample_11, "Sample"] <- "tdT- Small Wound 8"

jpeg(file = "tsne by sample #.jpeg", width = 30, height = 20, units = "cm", res = 500)
ggplot(T_tSNE_sample, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()

## LOOKING AT tSNEs BY SAMPLES - Schwann Cells, pooled ##
sch_tSNE_sample <- as.data.frame(GetCellEmbeddings(sub_sub.Schwann_cells, reduction.type = "tsne"))
sch_tSNE_sample[sample_2, "Sample"] <- "tdT- Large Wound Center"
sch_tSNE_sample[sample_4, "Sample"] <- "tdT- Large Wound Outside"
sch_tSNE_sample[sample_6, "Sample"] <- "tdT- Small Wound 14"
sch_tSNE_sample[sample_11, "Sample"] <- "tdT- Small Wound 8"

jpeg(file = "tsne by sample #.jpeg", width = 30, height = 20, units = "cm", res = 500)
ggplot(sch_tSNE_sample, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()

## LOOKING AT tSNEs BY SAMPLES - Endothelial cells, pooled ##
end_tSNE_sample <- as.data.frame(GetCellEmbeddings(sub_sub.Endothelial_cells, reduction.type = "tsne"))
end_tSNE_sample[sample_2, "Sample"] <- "tdT- Large Wound Center"
end_tSNE_sample[sample_4, "Sample"] <- "tdT- Large Wound Outside"
end_tSNE_sample[sample_6, "Sample"] <- "tdT- Small Wound 14"
end_tSNE_sample[sample_11, "Sample"] <- "tdT- Small Wound 8"

jpeg(file = "tsne by sample #.jpeg", width = 30, height = 20, units = "cm", res = 500)
ggplot(end_tSNE_sample, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()


## LOOKING AT tSNEs BY SAMPLES - Skeletal Muscle Progenitors, pooled ##
ske_tSNE_sample <- as.data.frame(GetCellEmbeddings(sub_sub.skeletal, reduction.type = "tsne"))
ske_tSNE_sample[sample_2, "Sample"] <- "tdT- Large Wound Center"
ske_tSNE_sample[sample_4, "Sample"] <- "tdT- Large Wound Outside"
ske_tSNE_sample[sample_6, "Sample"] <- "tdT- Small Wound 14"
ske_tSNE_sample[sample_11, "Sample"] <- "tdT- Small Wound 8"

jpeg(file = "tsne by sample #.jpeg", width = 30, height = 20, units = "cm", res = 500)
ggplot(ske_tSNE_sample, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()

## LOOKING AT tSNEs BY SAMPLES - Neutrophils, pooled ##
neu_tSNE_sample <- as.data.frame(GetCellEmbeddings(sub_sub.Neutrophils, reduction.type = "tsne"))
neu_tSNE_sample[sample_2, "Sample"] <- "tdT- Large Wound Center"
neu_tSNE_sample[sample_4, "Sample"] <- "tdT- Large Wound Outside"
neu_tSNE_sample[sample_6, "Sample"] <- "tdT- Small Wound 14"
neu_tSNE_sample[sample_11, "Sample"] <- "tdT- Small Wound 8"

jpeg(file = "tsne by sample #.jpeg", width = 30, height = 20, units = "cm", res = 500)
ggplot(neu_tSNE_sample, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()

## LOOKING AT tSNEs BY SAMPLES - Melanocytes, pooled ##
mel_tSNE_sample <- as.data.frame(GetCellEmbeddings(sub_sub.Melanocytes, reduction.type = "tsne"))
mel_tSNE_sample[sample_2, "Sample"] <- "tdT- Large Wound Center"
mel_tSNE_sample[sample_4, "Sample"] <- "tdT- Large Wound Outside"
mel_tSNE_sample[sample_6, "Sample"] <- "tdT- Small Wound 14"
mel_tSNE_sample[sample_11, "Sample"] <- "tdT- Small Wound 8"

jpeg(file = "tsne by sample #.jpeg", width = 30, height = 20, units = "cm", res = 500)
ggplot(mel_tSNE_sample, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()


## LOOKING AT tSNEs BY SAMPLES - Dendritic Cells, pooled ##
den_tSNE_sample <- as.data.frame(GetCellEmbeddings(sub_sub.Dendritic, reduction.type = "tsne"))
den_tSNE_sample[sample_2, "Sample"] <- "tdT- Large Wound Center"
den_tSNE_sample[sample_4, "Sample"] <- "tdT- Large Wound Outside"
den_tSNE_sample[sample_6, "Sample"] <- "tdT- Small Wound 14"
den_tSNE_sample[sample_11, "Sample"] <- "tdT- Small Wound 8"

jpeg(file = "tsne by sample #.jpeg", width = 30, height = 20, units = "cm", res = 500)
ggplot(den_tSNE_sample, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()

### Saving and Loading Each Population Files ### 
save(sub_sub.Dendritic, file = "sub_sub.Dendritic.Robj")
save(sub_sub.Endothelial_cells, file = "sub_sub.Endothelial_cells.Robj")
save(sub_sub.epi, file = "sub_sub.epi.Robj")
save(sub_sub.fibroblasts, file = "sub_sub.fibroblasts.Robj")
save(sub_sub.langerhans, file = "sub_sub.langerhans.Robj")
save(sub_sub.macro, file = "sub_sub.macro.Robj")
save(sub_sub.Melanocytes, file = "sub_sub.Melanocytes.Robj")
save(sub_sub.Neutrophils, file = "sub_sub.Neutrophils.Robj")
save(sub_sub.pericytes, file = "sub_sub.pericytes.Robj")
save(sub_sub.Schwann_cells, file = "sub_sub.Schwann_cells.Robj")
save(sub_sub.skeletal, file = "sub_sub.skeletal.Robj")
save(sub_sub.Tcells, file = "sub_sub.Tcells.Robj")
#
load(file = "sub_sub.Dendritic.Robj")
load(file = "sub_sub.Endothelial_cells.Robj")
load(file = "sub_sub.epi.Robj")
load(file = "sub_sub.fibroblasts.Robj")
load(file = "sub_sub.langerhans.Robj")
load(file = "sub_sub.macro.Robj")
load(file = "sub_sub.Melanocytes.Robj")
load(file = "sub_sub.Neutrophils.Robj")
load(file = "sub_sub.pericytes.Robj")
load(file = "sub_sub.Schwann_cells.Robj")
load(file = "sub_sub.skeletal.Robj")
load(file = "sub_sub.Tcells.Robj")
### Saving and Loading Files ### 


#### Finding a pan-Schwann Cell marker to do further peudotime analysis in Monocle:
Schwann.markers <- FindMarkers(wound_neg_0.3, ident.1 = c("Schwann Cells"), min.pct = 0.25, test.use = "negbinom")
Schwann.markers <- Schwann.markers[order(Schwann.markers$avg_logFC),]
write.csv(Schwann.markers,file = "Schwann.markers.csv")
FeaturePlot(wound_neg_0.3, features.plot = c("Plp1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
#### Plp1 is a pan-Schwann Cell Specific Marker


#### Finding a pan-T Cell marker to do further peudotime analysis in Monocle:
T_cell.markers <- FindMarkers(wound_neg_0.3, ident.1 = c("T Cells"), min.pct = 0.25, test.use = "negbinom")
T_cell.markers <- T_cell.markers[order(T_cell.markers$avg_logFC),]
write.csv(T_cell.markers,file = "T_cell.markers.csv")
FeaturePlot(wound_neg_0.3, features.plot = c(""), reduction.use = "tsne", cols.use = c("grey", "blue"))
####  is a pan-Schwann Cell Specific Marker


#### Finding a pan-T Cell marker to do further peudotime analysis in Monocle:
T_cell.markers <- FindMarkers(wound_neg_0.3, ident.1 = c("T Cells"), min.pct = 0.25, test.use = "negbinom")
T_cell.markers <- T_cell.markers[order(T_cell.markers$avg_logFC),]
write.csv(T_cell.markers,file = "T_cell.markers.csv")
FeaturePlot(wound_neg_0.3, features.plot = c("Cd3d"), reduction.use = "tsne", cols.use = c("grey", "blue"))
####  is a pan-Schwann Cell Specific Marker Cd3d

#### Finding a pan-T Cell marker to do further peudotime analysis in Monocle:
T_cell.markers <- FindMarkers(wound_neg_0.3, ident.1 = c("T Cells"), min.pct = 0.25, test.use = "negbinom")
T_cell.markers <- T_cell.markers[order(T_cell.markers$avg_logFC),]
write.csv(T_cell.markers,file = "T_cell.markers.csv")
FeaturePlot(wound_neg_0.3, features.plot = c("Cd3d"), reduction.use = "tsne", cols.use = c("grey", "blue"))
####  is a pan-Schwann Cell Specific Marker Cd3d


top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)



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

## Looking at macrophases at a Higher Resolution to isolate regional specific clusters
macro_recluster_4 <- FindClusters(macro_recluster, reduction.type = "pca", dims.use = 1:22, resolution = 4, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
macro_recluster_4 <- StashIdent(macro_recluster_4, save.name = "ClusterNames_macro_recluster_4")
macro_recluster_4 <- RunTSNE(macro_recluster_4, dims.use = 1:22, do.fast = T)

jpeg(file = "macro.jpeg", width = 15, height = 15, units = "cm", res = 500)
TSNEPlot(macro_recluster_4, do.label = FALSE, pt.size = 0.8)
dev.off()

## Save and Load
save(macro_recluster, file = "macro_recluster.Robj")
save(macro_recluster_4, file = "macro_recluster_4.Robj")
load(file = "macro_recluster.Robj")
load(file = "macro_recluster_4.Robj")
##

# LWC (3, 0) Macro Markers
LWC_macro.markers <- FindMarkers(macro_recluster_4, ident.1 = c(3, 0), min.pct = 0.25, test.use = "negbinom")
LWC_macro.markers <- LWC_macro.markers[order(LWC_macro.markers$avg_logFC),]
write.csv(LWC_macro.markers,file = "LWC_macro_markers.csv")
# LWO (2) Macro Markers
LWO_macro.markers <- FindMarkers(macro_recluster_4, ident.1 = c(2), min.pct = 0.25, test.use = "negbinom")
LWO_macro.markers <- LWO_macro.markers[order(LWO_macro.markers$avg_logFC),]
write.csv(LWO_macro.markers,file = "LWO_macro.markers.csv")
# SW8 (17, 13) Macro Markers
SW8_macro.markers <- FindMarkers(macro_recluster_4, ident.1 = c(17, 13), min.pct = 0.25, test.use = "negbinom")
SW8_macro.markers <- SW8_macro.markers[order(SW8_macro.markers$avg_logFC),]
write.csv(SW8_macro.markers,file = "SW8_macro.markers.csv")
# SW14 (15) Macro Markers
SW14_macro.markers <- FindMarkers(macro_recluster_4, ident.1 = c(15), min.pct = 0.25, test.use = "negbinom")
SW14_macro.markers <- SW14_macro.markers[order(SW14_macro.markers$avg_logFC),]
write.csv(SW14_macro.markers,file = "SW14_macro.markers.csv")

FeaturePlot(macro_recluster, features.plot = c("Fn1"), reduction.use = "tsne", cols.use = c("grey", "blue"))

macro_recluster_4@meta.data[1:776,]$sample <- "LWC"





FeatureHeatmap(macro.combined, features.plot = c("Tnf"), group.by = "stim", pt.size = 0.5, key.position = "top", 
               max.exp = 3)


# Looking at Macrophase Clusters
or_tSNE <- as.data.frame(GetCellEmbeddings(macro_recluster, reduction.type = "tsne"))
or_tSNE[sample_2, "Sample"] <- "D14 LWC"
or_tSNE[sample_4, "Sample"] <- "D14 LWP"
or_tSNE[sample_6, "Sample"] <- "D14 SW"
or_tSNE[sample_11, "Sample"] <- "D8 SW"

ggplot(or_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()

### Choosing 2,723 cells randomly from Samples 2, 4, and 6 to show distribution with even number of cells
random_lwc <- sample(1:775, 600, replace = FALSE)
random_lwo <- sample(1:682, 600, replace = FALSE)
random_sw8 <- sample(1:175, 175, replace = FALSE)
random_sw14 <- sample(1:237, 175, replace = FALSE)



# Assigning sample #s:
sample_2 <- grep("-2", colnames(macro_recluster@data), value = T); length(sample_2)
sample_4 <- grep("-4", colnames(macro_recluster@data), value = T); length(sample_4)
sample_6 <- grep("-6", colnames(macro_recluster@data), value = T); length(sample_6)
sample_11 <- grep("-11", colnames(macro_recluster@data), value = T); length(sample_11)



## Preparing Macrophases by sample:
lwc_macro <- FilterCells(macro_recluster, cells.use = c(sample_2), subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
lwc_macro <- NormalizeData(lwc_macro)
lwc_macro <- ScaleData(lwc_macro, display.progress = F)
lwc_macro@meta.data$stim <- "LWC"

lwo_macro <- FilterCells(macro_recluster, cells.use = c(sample_4), subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
lwo_macro <- NormalizeData(lwo_macro)
lwo_macro <- ScaleData(lwo_macro, display.progress = F)
lwo_macro@meta.data$stim <- "LWO"

sw8_macro <- FilterCells(macro_recluster, cells.use = c(sample_6), subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
sw8_macro <- NormalizeData(sw8_macro)
sw8_macro <- ScaleData(sw8_macro, display.progress = F)
sw8_macro@meta.data$stim <- "SW8"

sw14_macro <- FilterCells(macro_recluster, cells.use = c(sample_11), subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
sw14_macro <- NormalizeData(sw14_macro)
sw14_macro <- ScaleData(sw14_macro, display.progress = F)
sw14_macro@meta.data$stim <- "SW14"



current.cluster.ids <- c(3, 0, 2, 17, 13, 15)
new.cluster.ids <- c("LWC", "LWC", "LWO", "SW8", "SW8", "SW14")
macro_recluster_4@ident <- plyr::mapvalues(macro_recluster_4@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(macro_recluster_4, do.label = TRUE, pt.size = 0.5)

current.cluster.ids <- c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 18, 19, 20, 21, 22, 23)
new.cluster.ids <- c("Mix", "Mix", "Mix", "Mix", "Mix", "Mix", "Mix", "Mix", "Mix", "Mix", "Mix", "Mix", "Mix", "Mix", "Mix","Mix", "Mix", "Mix")
macro_recluster_4@ident <- plyr::mapvalues(macro_recluster_4@ident, from = current.cluster.ids, to = new.cluster.ids)
jpeg(file = "macro_of_interest.jpeg", width = 15, height = 15, units = "cm", res = 500)
TSNEPlot(macro_recluster_4, do.label = FALSE, pt.size = 0.5)
dev.off()

markers.to.plot <- c("Tnf","Tnfaip2", "Slpi", "Ptgs2", 
                     "Clec4d", "Thbs1", "F10", "Sod2", 
                     "Clec4e", "Cxcl3", "Sepp1", "Apoe", "Ms4a7", "C1qa",
                     "Ccl8", "Ccl12", "C1qc", "C1qb", "Pf4", "Maf",
                     "Ifrd1", "Mmp12", "Pxdc1", "Pgf", "Eif5", "Ptgr1", "Hilpda",
                     "Ndrg1", "Rgcc", "Gm26917", "Retnla", "Lyz1", "Ear2",
                     "Fcrls", "Ccl6", "Hspb1", "Clec4b1", "Pltp", "Phlda1", "Lpl")

jpeg(file = "macro_of_interest_genes.jpeg", width = 40, height = 15, units = "cm", res = 500)
sdp <- SplitDotPlotGG(macro_recluster_4, 
                      genes.plot = rev(markers.to.plot),
                      cols.use = c("#F85665", "khaki4","#0EC068","#01A9EF","#DA4DF3"),
                      x.lab.rot = T, 
                      plot.legend = T, 
                      dot.scale = 8, 
                      do.return = T, 
                      grouping.var = "ident")
dev.off()


markers.to.plot_trial <- c("Tnf", "Il23a", "Csf3", "Il6", "Ccl4", 
                           "Il10", "Cx3cr1", "Ccr2", "Map3k5", "Igf1", "Fgf2", "Mrc1")
sdp_trial <- SplitDotPlotGG(macro_recluster_4, 
                      genes.plot = rev(markers.to.plot_trial),
                      cols.use = c("#F85665", "khaki4","#0EC068","#01A9EF","#DA4DF3"),
                      x.lab.rot = T, 
                      plot.legend = T, 
                      dot.scale = 8, 
                      do.return = T, 
                      grouping.var = "ident")






###CURRENT





# Gene selection for input to CCA
ctrl <- FindVariableGenes(ctrl, do.plot = F)
stim <- FindVariableGenes(stim, do.plot = F)
g.1 <- head(rownames(ctrl@hvg.info), 1000)
g.2 <- head(rownames(stim@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(ctrl@scale.data))
genes.use <- intersect(genes.use, rownames(stim@scale.data))







a_mac_recluster <- RunPCA(lwc_macro, lwo_macro, sw8_macro, sw14_macro,  pc.genes = macro_recluster@var.genes, pcs.compute = 50)
PCElbowPlot(macro_recluster, num.pc = 50) ## 22 significant PCs
macro_recluster <- FindClusters(macro_recluster, reduction.type = "pca", dims.use = 1:22, resolution = 1.5, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
macro_recluster <- StashIdent(macro_recluster, save.name = "ClusterNames_macro_recluster")
macro_recluster <- RunTSNE(macro_recluster, dims.use = 1:22, do.fast = T)
TSNEPlot(macro_recluster, do.label = TRUE, pt.size = 0.8)





lwc_macro <- FindVariableGenes(lwc_macro, do.plot = F)
lwo_macro <- FindVariableGenes(lwo_macro, do.plot = F)
g.1 <- head(rownames(lwc_macro@hvg.info), 1000)
g.2 <- head(rownames(lwo_macro@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(lwc_macro@scale.data))
genes.use <- intersect(genes.use, rownames(lwo_macro@scale.data))

macro.combined <- RunCCA(lwc_macro, lwo_macro, genes.use = genes.use, num.cc = 30)


p1 <- DimPlot(object = macro.combined, reduction.use = "cca", group.by = "stim", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = macro.combined, features.plot = "CC1", group.by = "stim", 
              do.return = TRUE)
plot_grid(p1, p2)

PrintDim(macro.combined, reduction.type = "cca", dims.print = 1:2, 
         genes.print = 10)

p3 <- MetageneBicorPlot(macro.combined, grouping.var = "stim", dims.eval = 1:30, 
                        display.progress = FALSE)



macro.combined <- AlignSubspace(macro.combined, reduction.type = "cca", grouping.var = "stim", 
                                 dims.align = 1:25)

p1 <- VlnPlot(object = macro.combined, features.plot = "ACC1", group.by = "stim", 
              do.return = TRUE)
p2 <- VlnPlot(object = macro.combined, features.plot = "ACC2", group.by = "stim", 
              do.return = TRUE)
plot_grid(p1, p2)

macro.combined <- RunTSNE(macro.combined, reduction.use = "cca.aligned", dims.use = 1:25, do.fast = T)
macro.combined <- FindClusters(macro.combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:25)

p1 <- TSNEPlot(macro.combined, do.return = T, pt.size = 0.5, group.by = "stim")
p2 <- TSNEPlot(macro.combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)



markers.to.plot <- c("Tnf")

sdp <- SplitDotPlotGG(macro.combined, 
                      genes.plot = rev(markers.to.plot), 
                      cols.use = c("blue","red"), 
                      x.lab.rot = T, 
                      plot.legend = T, 
                      dot.scale = 8, 
                      do.return = T, 
                      grouping.var = "stim")


LabelPoint <- function(plot, genes, exp.mat, adj.x.t = 0, adj.y.t = 0, adj.x.s = 0, 
                       adj.y.s = 0, text.size = 2.5, segment.size = 0.1) {
  for (i in genes) {
    x1 <- exp.mat[i, 1]
    y1 <- exp.mat[i, 2]
    plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t, 
                            label = i, size = text.size)
    plot <- plot + annotate("segment", x = x1 + adj.x.s, xend = x1, y = y1 + 
                              adj.y.s, yend = y1, size = segment.size)
  }
  return(plot)
}

LabelUR <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.r.t = 0.15, adj.u.s = 0.05, 
                    adj.r.s = 0.05, ...) {
  return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = adj.r.t, 
                    adj.y.s = adj.u.s, adj.x.s = adj.r.s, ...))
}

LabelUL <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05, 
                    adj.l.s = 0.05, ...) {
  return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = -adj.l.t, 
                    adj.y.s = adj.u.s, adj.x.s = -adj.l.s, ...))
}



macro.cells <- SubsetData(macro.combined, subset.raw = T)
macro.cells <- SetAllIdent(macro.cells, id = "stim")
avg.macro.cells <- log1p(AverageExpression(macro.cells, show.progress = FALSE))
avg.macro.cells$gene <- rownames(avg.macro.cells)

macro.mono <- SubsetData(macro.combined, subset.raw = T)
macro.mono <- SetAllIdent(macro.mono, id = "stim")
avg.macro.mono <- log1p(AverageExpression(macro.mono, show.progress = FALSE))
avg.macro.mono$gene <- rownames(avg.macro.mono)

genes.to.label1 = c("Tnf")
genes.to.label2 = c("Tnf")
genes.to.label3 = c("Tnf")
p1 <- ggplot(avg.macro.cells, aes(LWC, LWO)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelUR(p1, genes = c(genes.to.label1, genes.to.label2), avg.macro.cells, 
              adj.u.t = 0.3, adj.u.s = 0.23)
p1 <- LabelUL(p1, genes = genes.to.label3, avg.macro.cells, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25)
p2 <- ggplot(avg.macro.mono, aes(LWC, LWO)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelUR(p2, genes = c(genes.to.label1, genes.to.label3), avg.macro.mono, 
              adj.u.t = 0.3, adj.u.s = 0.23)
p2 <- LabelUL(p2, genes = genes.to.label2, avg.macro.mono, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25)
plot_grid(p1, p2)

FeatureHeatmap(macro.combined, features.plot = c("Tnf"), group.by = "stim", pt.size = 0.5, key.position = "top", 
               max.exp = 3)






sw14_macro <- FilterCells(macro_recluster, cells.use = c(sample_6), subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
sw8_macro <- FilterCells(macro_recluster, cells.use = c(sample_11), subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)

#### SAVE!!!
jpeg(file = "macro_recluster_lwc", width = 15, height = 15, units = "cm", res = 500)
TSNEPlot(macro_recluster, cells.use = lwc_macro@cell.names[random_lwc])
dev.off()

jpeg(file = "macro_recluster_lwo", width = 15, height = 15, units = "cm", res = 500)
TSNEPlot(macro_recluster, cells.use = lwo_macro@cell.names[random_lwo])
dev.off()

jpeg(file = "macro_recluster", width = 30, height = 15, units = "cm", res = 500)
TSNEPlot(macro_recluster, cells.use = sw14_macro@cell.names[random_sw8])
dev.off()

jpeg(file = "macro_recluster", width = 30, height = 15, units = "cm", res = 500)
TSNEPlot(macro_recluster, cells.use = sw8_macro@cell.names[random_sw14])
dev.off()

jpeg(file = "macro_recluster", width = 15, height = 15, units = "cm", res = 500)
TSNEPlot(macro_recluster, do.label = FALSE, pt.size = 0.8)
dev.off()


jpeg(file = "macro_recluster_by_sample", width = 15, height = 15, units = "cm", res = 500)
ggplot(or_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()








length(lwc);length(lwo);length(sw);length(sw8)
random_lwc<- sample(1:775, 600, replace = FALSE)
random_m2 <- sample(1:682, 600, replace = FALSE) 
random_m3 <- sample(1:175, 1500, replace = FALSE)
random_m4 <- sample(1:236, 1500, replace = FALSE)

# Not working yet
TSNEPlot(macro_recluster, cells.use = lwc@cell.names[random_m1])
TSNEPlot(macro_recluster, cells.use = lwo@cell.names[random_m2])
TSNEPlot(macro_recluster, cells.use = sw@cell.names[random_m1])
TSNEPlot(macro_recluster, cells.use = sw8@cell.names[random_m1])

TSNEPlot(macro_recluster, cells.use = lwc@cell.names[random_lwc])
TSNEPlot(macro_recluster, cells.use = lwo@cell.names[random_lwo])





########## T cells from wound_neg_0.3 and Re-Clustering ######
tcells <- names(wound_neg_0.3@ident[wound_neg_0.3@ident %in% c("T Cells")])
sub.t <- FilterCells(wound_neg_0.3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = tcells)
t_recluster <- FilterCells(sub.t, cells.use = c(sample_2, sample_4, sample_6, sample_11), subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
t_recluster <- FindVariableGenes(t_recluster, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
length(x = t_recluster@var.genes)
t_recluster <- NormalizeData(t_recluster, normalization.method = "LogNormalize", scale.factor = 10000)
t_recluster <- ScaleData(t_recluster, vars.to.regress = c("nUMI", "percent.mito"))
t_recluster <- RunPCA(t_recluster, pc.genes = t_recluster@var.genes, pcs.compute = 50)
PCElbowPlot(t_recluster, num.pc = 50) ## 22 significant PCs
t_recluster <- FindClusters(t_recluster, reduction.type = "pca", dims.use = 1:22, resolution = 1.5, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
t_recluster <- StashIdent(t_recluster, save.name = "ClusterNames_t_recluster")
t_recluster <- RunTSNE(t_recluster, dims.use = 1:22, do.fast = T)
TSNEPlot(t_recluster, do.label = TRUE, pt.size = 0.8)


# Looking at T Clusters
t_tSNE <- as.data.frame(GetCellEmbeddings(t_recluster, reduction.type = "tsne"))
t_tSNE[sample_2, "Sample"] <- "D14 LWC"
t_tSNE[sample_4, "Sample"] <- "D14 LWP"
t_tSNE[sample_6, "Sample"] <- "D14 SW"
t_tSNE[sample_11, "Sample"] <- "D8 SW"

ggplot(t_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
FeaturePlot(t_recluster, features.plot = c("Fgf9"), reduction.use = "tsne", cols.use = c("grey", "blue"))


########## Schwann cells from wound_neg_0.3 and Re-Clustering ######
schwann_cells <- names(wound_neg_0.3@ident[wound_neg_0.3@ident %in% c("Schwann Cells")])
sub.sc <- FilterCells(wound_neg_0.3, subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000, cells.use = schwann_cells)
sc_recluster <- FilterCells(sub.sc, cells.use = c(sample_2, sample_4, sample_6, sample_11), subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
sc_recluster <- FindVariableGenes(sc_recluster, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 8, y.cutoff = 1, do.plot = TRUE)
length(x = sc_recluster@var.genes)
sc_recluster <- NormalizeData(sc_recluster, normalization.method = "LogNormalize", scale.factor = 10000)
sc_recluster <- ScaleData(sc_recluster, vars.to.regress = c("nUMI", "percent.mito"))
sc_recluster <- RunPCA(sc_recluster, pc.genes = sc_recluster@var.genes, pcs.compute = 50)
PCElbowPlot(sc_recluster, num.pc = 50) ## 22 significant PCs
sc_recluster <- FindClusters(sc_recluster, reduction.type = "pca", dims.use = 1:22, resolution = 1.5, print.output = 1, save.SNN = TRUE, force.recalc = TRUE)
sc_recluster <- StashIdent(sc_recluster, save.name = "ClusterNames_sc_recluster")
sc_recluster <- RunTSNE(sc_recluster, dims.use = 1:22, do.fast = T)
TSNEPlot(sc_recluster, do.label = TRUE, pt.size = 0.8)

# Looking at T Clusters
sc_tSNE <- as.data.frame(GetCellEmbeddings(sc_recluster, reduction.type = "tsne"))
sc_tSNE[sample_2, "Sample"] <- "D14 LWC"
sc_tSNE[sample_4, "Sample"] <- "D14 LWP"
sc_tSNE[sample_6, "Sample"] <- "D14 SW"
sc_tSNE[sample_11, "Sample"] <- "D8 SW"



ggplot(sc_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample)) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
FeaturePlot(t_recluster, features.plot = c("Fgf9"), reduction.use = "tsne", cols.use = c("grey", "blue"))









#################Canonical Correlation Analysis#################################

lwc_wound <- FilterCells(wound_neg_0.3, cells.use = c(sample_2), subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
lwc_wound <- NormalizeData(lwc_wound)
lwc_wound <- ScaleData(lwc_wound, display.progress = F)
lwc_wound@meta.data$stim <- "LWC"

lwo_wound <- FilterCells(wound_neg_0.3, cells.use = c(sample_4), subset.names = "nGene", low.thresholds = 200, high.thresholds = 20000)
lwo_wound <- NormalizeData(lwo_wound)
lwo_wound <- ScaleData(lwo_wound, display.progress = F)
lwo_wound@meta.data$stim <- "LWO"


lwc_wound <- FindVariableGenes(lwc_wound, do.plot = F)
lwo_wound <- FindVariableGenes(lwo_wound, do.plot = F)
g.1 <- head(rownames(lwc_wound@hvg.info), 1000)
g.2 <- head(rownames(lwo_wound@hvg.info), 1000)
genes.use <- unique(c(g.1, g.2))
genes.use <- intersect(genes.use, rownames(lwc_wound@scale.data))
genes.use <- intersect(genes.use, rownames(lwo_wound@scale.data))

wound.combined <- RunCCA(lwc_wound, lwo_wound, genes.use = genes.use, num.cc = 30)
p1 <- DimPlot(wound.combined, reduction.use = "cca", group.by = "stim", pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(wound.combined, features.plot = "CC1", group.by = "stim", do.return = TRUE)
plot_grid(p1, p2)
PrintDim(object = wound.combined, reduction.type = "cca", dims.print = 1:2, 
         genes.print = 10)
p3 <- MetageneBicorPlot(wound.combined, grouping.var = "stim", dims.eval = 1:30, display.progress = TRUE)
wound.combined <- AlignSubspace(wound.combined, reduction.type = "cca", grouping.var = "stim", dims.align = 1:25)
p1 <- VlnPlot(wound.combined, features.plot = "ACC1", group.by = "stim", 
              do.return = TRUE)
p2 <- VlnPlot(wound.combined, features.plot = "ACC2", group.by = "stim", 
              do.return = TRUE)
plot_grid(p1, p2)


wound.combined <- RunTSNE(wound.combined, reduction.use = "cca.aligned", dims.use = 1:25, do.fast = T)
wound.combined <- FindClusters(wound.combined, reduction.type = "cca.aligned", resolution = 0.6, dims.use = 1:25)

p1 <- TSNEPlot(wound.combined, do.return = T, pt.size = 0.5, group.by = "stim")
p2 <- TSNEPlot(wound.combined, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)

FeaturePlot(wound.combined, features.plot = c("Krt14", "Krt5", "Cdh1"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
FeaturePlot(wound.combined, features.plot = c("Cd68"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
FeaturePlot(wound.combined, features.plot = c("Cd207"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
FeaturePlot(wound.combined, features.plot = c("Sox10"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
FeaturePlot(wound.combined, features.plot = c("Ngp"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
FeaturePlot(wound.combined, features.plot = c("Cdh15"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
FeaturePlot(wound.combined, features.plot = c("Cd3g"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
FeaturePlot(wound.combined, features.plot = c("Il22ra2"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)

clus10.markers <- FindConservedMarkers(wound.combined, ident.1 = 10, grouping.var = "stim", 
                                   print.bar = FALSE)
head(clus10.markers)

## Save and Load
save(wound.combined, file = "wound.combined.Robj")
save(macro.mono, file = "macro.mono.Robj")
save(t.cells, file = "wound.combined.Robj")
load(file = "wound.combined.Robj")
##


current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 15, 17, 18)
new.cluster.ids <- c("Epi", "Epi", "Macro", "Macro", "T_cells", "Schwann", "Macro", "Fibro", 
                     "Epi", "Fibro", "Langerhans", "Fibro", "Epi", "Neutrophils")
wound.combined@ident <- plyr::mapvalues(wound.combined@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(wound.combined, do.label = TRUE, pt.size = 0.5)

## Looking at T Cells and Macro
t.cells <- SubsetData(wound.combined, ident.use = "T_cells", subset.raw = T)
t.cells <- SetAllIdent(t.cells, id = "stim")
avg.t.cells <- log1p(AverageExpression(t.cells, show.progress = FALSE))
avg.t.cells$gene <- rownames(avg.t.cells)

macro.mono <- SubsetData(wound.combined, ident.use = "Macro", subset.raw = T)
macro.mono <- SetAllIdent(macro.mono, id = "stim")
avg.macro.mono <- log1p(AverageExpression(macro.mono, show.progress = FALSE))
avg.macro.mono$gene <- rownames(avg.macro.mono)

genes.to.label1 = c("Tnf", "Sod2", "Mmp9")
genes.to.label2 = c("Fgf9")
genes.to.label3 = c("Sepp1")
p1 <- ggplot(avg.t.cells, aes(LWC, LWO)) + geom_point() + ggtitle("T Cells")
p1 <- LabelUR(p1, genes = c(genes.to.label1, genes.to.label2), avg.t.cells, 
              adj.u.t = 0.3, adj.u.s = 0.23)
p1 <- LabelUL(p1, genes = genes.to.label3, avg.t.cells, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25)
p2 <- ggplot(avg.macro.mono, aes(LWC, LWO)) + geom_point() + ggtitle("Macrophases")
p2 <- LabelUR(p2, genes = c(genes.to.label1, genes.to.label3), avg.macro.mono, 
              adj.u.t = 0.3, adj.u.s = 0.23)
p2 <- LabelUL(p2, genes = genes.to.label2, avg.macro.mono, adj.u.t = 0.5, adj.u.s = 0.4, 
              adj.l.t = 0.25, adj.l.s = 0.25)
plot_grid(p1, p2)



wound.combined@meta.data$celltype.stim <- paste0(wound.combined@ident, "_", 
                                                 wound.combined@meta.data$stim)
wound.combined <- StashIdent(wound.combined, save.name = "celltype")
wound.combined <- SetAllIdent(wound.combined, id = "celltype.stim")
macro.response <- FindMarkers(wound.combined, ident.1 = "Macro_LWC_LWC", ident.2 = "Macro_LWO_LWO", print.bar = FALSE)
macro.response <- macro.response[order(macro.response$avg_logFC),]
write.csv(macro.response,file = "macro.response.csv")

jpeg(file = "FeatureHeatmap_Macro", width = 15, height = 40, units = "cm", res = 500)
FeatureHeatmap(macro.mono, features.plot = c("Tnf", "Cxcl1", "Cxcl2", "Il10", "Ifitm1", "Tnfaip2", "Cxcl3", "Il1b", "Clec4e", "Sod2"), group.by = "stim", pt.size = 0.25, key.position = "top", max.exp = 3)
dev.off()

jpeg(file = "FeatureHeatmap_wound.combined", width = 15, height = 40, units = "cm", res = 500)
FeatureHeatmap(wound.combined, features.plot = c("Tnf", "Cxcl1", "Cxcl2", "Il10", "Ifitm1", "Tnfaip2", "Cxcl3", "Il1b", "Clec4e", "Sod2"), group.by = "stim", pt.size = 0.25, key.position = "top", max.exp = 3)
dev.off()



### obtain the data matrix and subset out only for cells of interest from sub_sub.macro
## For PAGODA:
LWC_14 <- grep("-2", colnames(sub_sub.macro@data), value = T); length(LWC_14)
LWO_14 <- grep("-4", colnames(sub_sub.macro@data), value = T); length(LWO_14)
SW_14 <- grep("-6", colnames(sub_sub.macro@data), value = T); length(SW_14)
SW_8 <- grep("-11", colnames(sub_sub.macro@data), value = T); length(SW_8)

cells <- c(LWC_14, LWO_14, SW_14 , SW_8)


q_df = as.data.frame.array(sub_sub.macro@data)
q_mx <- as.matrix(q_df)

s_mx_cells <- s_mx[,cells]


save(q_mx, file = "q_mx.Robj")
save(s_mx_cells, file = "s_mx_cells.Robj")

load("s_mx_prajay.Robj")




FeaturePlot(wound_neg_0.3, features.plot = c("Cd68", "Ccr2"), reduction.use = "tsne", cols.use = c("grey", "green", "red"), overlay = T, no.legend = F)

pdf("trial.pdf")
RidgePlot(wound_neg_0.3, features.plot = c("Cd68"))
dev.off()



### Trying out PHATE:
install.packages("phateR", lib="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
library(phateR)
load(file = "epi_recluster_LW.Robj")
epi_recluster_LW = RunPHATE(epi_recluster_LW, k = 3, t = 12)
DimPlot(epi_recluster_LW, reduction.use = 'phate')
DimPlot(epi_recluster_LW, reduction.use = 'phate')


FeaturePlot(epi_recluster_LW, features.plot = c("Ifitm3"), reduction.use = "phate", cols.use = c("grey", "blue"))


save(epi_recluster_LW, file = "epi_recluster_LW.Robj")

pmbc_small <- RefinedMapping(pbmc_small, genes.use=pbmc_small@var.genes)


ColorTSNESplit(epi_recluster_LW, 1, color1 = "red", color2 = "blue")



########################## Data Export for Cell Phone Database: ##########################
load(file = "epi_recluster_LW.Robj")
View(wound_neg_0.3@data)
epi_recluster_LW_subset <- SubsetData(epi_recluster_LW, subset.raw = T)
epi_recluster_LW_CPD = as.data.frame.array(epi_recluster_LW_subset@raw.data)
epi_recluster_LW_CPD <- as.matrix(epi_recluster_LW_CPD)
write.csv(epi_recluster_LW_CPD, file = "epi_recluster_LW_CPD.csv")
q_mx <- as.matrix(q_df)
save(q_mx, file = "q_mx.Robj")

load(file = "wound_neg_0.3.Robj")
wound_neg_0.3 <- StashIdent(wound_neg_0.3, save.name = "Cell_names")
CPD = as.data.frame.array(wound_neg_0.3@data)
CPD <- as.matrix(CPD)
write.csv(CPD, file = "CPD.csv")

CPD_2 = as.data.frame.array(wound_neg_0.3@meta.data)
CPD_2 <- as.matrix(CPD_2)
write.csv(CPD_2, file = "CPD_2.csv")

######### Day 8 SW tdTom+ve and tdTom-ve #########
load(file = "wound.Robj")
sample_1 <- grep("-1", colnames(wound@data), value = T); length(sample_1)
sample_2 <- grep("-2", colnames(wound@data), value = T); length(sample_2)
sample_3 <- grep("-3", colnames(wound@data), value = T); length(sample_3)
sample_4 <- grep("-4", colnames(wound@data), value = T); length(sample_4)
sample_5 <- grep("-5", colnames(wound@data), value = T); length(sample_5)
sample_6 <- grep("-6", colnames(wound@data), value = T); length(sample_6)
sample_7 <- grep("-7", colnames(wound@data), value = T); length(sample_7)
Day18_Day8_Lib1_try1 <- grep("-8", colnames(wound@data), value = T); length(Day18_Day8_Lib1_try1)
Day18_Day8_Lib2_try1 <- grep("-9", colnames(wound@data), value = T); length(Day18_Day8_Lib2_try1)
Day18_Day8_Lib3_try1 <- grep("-10", colnames(wound@data), value = T); length(Day18_Day8_Lib3_try1)
Day18_Day8_Lib4_try1 <- grep("-11", colnames(wound@data), value = T); length(Day18_Day8_Lib4_try1)

D8_SW <- FilterCells(wound, cells.use = c(Day18_Day8_Lib3_try1, Day18_Day8_Lib4_try1), subset.names = "nGene", low.thresholds = 400, high.thresholds = 5500)
D8_SW <- NormalizeData(object = D8_SW, normalization.method = "LogNormalize", scale.factor = 10000)
D8_SW <- FindVariableGenes(D8_SW, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = D8_SW@var.genes) #2548
D8_SW <- ScaleData(D8_SW, vars.to.regress = c("nUMI", "percent.mito"))
D8_SW <- RunPCA(object = D8_SW, pc.genes = D8_SW@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 50)
PrintPCA(object = D8_SW, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(D8_SW, pcs.use = 1:2)
PCAPlot(D8_SW, dim.1 =1, dim.2 = 2)
PCHeatmap(D8_SW, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCElbowPlot(D8_SW, num.pc = 50) # 25 PCs
D8_SW <- FindClusters(D8_SW, reduction.type = "pca", dims.use = 1:25, resolution = 0.6, print.output = 0, save.SNN = TRUE)
D8_SW <- StashIdent(D8_SW, save.name = "D8_SW_0.6")
D8_SW <- RunTSNE(D8_SW, dims.use = 1:25, do.fast = T)
TSNEPlot(D8_SW, do.label = TRUE, pt.size = 0.8)
## 0 5 10 Epi; 12, 11 macrophases; 3 langerhans; 9, 2, 7, 1, 4, 13 fibros; 8 endothelial cells; 16 Schwann; 
## 14 Tcells; 6 pericytes; 
FeaturePlot(D8_SW, features.plot = c("Krt14"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D8_SW, features.plot = c("Krt28"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D8_SW, features.plot = c("Cd68"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D8_SW, features.plot = c("Cd207"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D8_SW, features.plot = c("Pdgfra"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D8_SW, features.plot = c("Pecam1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D8_SW, features.plot = c("Sox10"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D8_SW, features.plot = c("Cd3g"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D8_SW, features.plot = c("Acta2"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D8_SW, features.plot = c("Notch3"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D8_SW, features.plot = c("Pecam1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D8_SW, features.plot = c("Cox6a2"), reduction.use = "tsne", cols.use = c("grey", "blue"))
current <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
new <- c("Epithelial cells",
         "Fibroblasts", 
         "Fibroblasts", 
         "Langerhans", 
         "Fibroblasts", 
         "Epithelial cells",
         "Pericytes",
         "Fibroblasts", 
         "Endothelial Cells",
         "Fibroblasts",
         "Epithelial cells",
         "Macrophases",
         "Macrophases",
         "Fibroblasts",
         "T cells",
         "Endothelial Cells",
         "Schwann Cells")
D8_SW@ident <- plyr::mapvalues(D8_SW@ident, from = current, to = new)
TSNEPlot(D8_SW, do.label = TRUE, pt.size = 0.5)
D8_SW <- StashIdent(D8_SW, save.name = "Cell_names")
D8_SW <- SubsetData(D8_SW, subset.raw = T)
save(D8_SW, file = "D8_SW.Robj")

CPD_D8_SW = as.data.frame.array(D8_SW@data)
CPD_D8_SW <- as.matrix(CPD_D8_SW)
write.csv(CPD_D8_SW, file = "CPD_D8_SW.csv")

CPDmeta_D8_SW <- as.data.frame.array(D8_SW@meta.data)
CPDmeta_D8_SW <- as.matrix(CPDmeta_D8_SW)
write.csv(CPDmeta_D8_SW, file = "CPDmeta_D8_SW.csv")

######### CELLPHONE: Day 14 SW tdTom+ve and tdTom-ve ######
D14_SW <- FilterCells(wound, cells.use = c(sample_5, sample_6), subset.names = "nGene", low.thresholds = 400, high.thresholds = 5500)
D14_SW <- NormalizeData(object = D14_SW, normalization.method = "LogNormalize", scale.factor = 10000)
D14_SW <- FindVariableGenes(D14_SW, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = D14_SW@var.genes) #2548
D14_SW <- ScaleData(D14_SW, vars.to.regress = c("nUMI", "percent.mito"))
D14_SW <- RunPCA(object = D14_SW, pc.genes = D14_SW@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 50)
PrintPCA(object = D14_SW, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(D14_SW, pcs.use = 1:2)
PCAPlot(D14_SW, dim.1 =1, dim.2 = 2)
PCHeatmap(D14_SW, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCElbowPlot(D14_SW, num.pc = 50) # 25 PCs
D14_SW <- FindClusters(D14_SW, reduction.type = "pca", dims.use = 1:25, resolution = 0.6, print.output = 0, save.SNN = TRUE)
D14_SW <- StashIdent(D14_SW, save.name = "D14_SW_0.6")
D14_SW <- RunTSNE(D14_SW, dims.use = 1:25, do.fast = T)
TSNEPlot(D14_SW, do.label = TRUE, pt.size = 0.8)
## 11 2 3 10 4 7 9 16 Epi; 6 macrophases; 1, 0, 5, 8 fibros; 18 15 endothelial cells; 14 Schwann; 
## 17 Tcells; 12 pericytes; 13 SKM Progenitors; 
FeaturePlot(D14_SW, features.plot = c("Krt14"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_SW, features.plot = c("Krt5"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_SW, features.plot = c("Krt28"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_SW, features.plot = c("Cdh1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_SW, features.plot = c("Cd68"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_SW, features.plot = c("Cd207"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_SW, features.plot = c("Pdgfra"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_SW, features.plot = c("Dpt"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_SW, features.plot = c("Pecam1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_SW, features.plot = c("Sox10"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_SW, features.plot = c("Cd3g"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_SW, features.plot = c("Acta2"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_SW, features.plot = c("Notch3"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_SW, features.plot = c("Pecam1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_SW, features.plot = c("Cox6a2"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_SW, features.plot = c("Mitf"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_SW, features.plot = c("Pax7"), reduction.use = "tsne", cols.use = c("grey", "blue"))

current <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
new <- c("Fibroblasts",
         "Fibroblasts",
         "Epithelial cells",
         "Epithelial cells",
         "Epithelial cells",
         "Fibroblasts",
         "Macrophases",
         "Epithelial cells",
         "Fibroblasts",
         "Epithelial cells",
         "Epithelial cells",
         "Epithelial cells",
         "Pericytes",
         "SKM Progenitors",
         "Schwann Cells",
         "Endothelial Cells",
         "Epithelial cells",
         "T cells",
         "Endothelial Cells")
D14_SW@ident <- plyr::mapvalues(D14_SW@ident, from = current, to = new)
TSNEPlot(D14_SW, do.label = TRUE, pt.size = 0.5)
D14_SW <- StashIdent(D14_SW, save.name = "Cell_names")
D14_SW <- SubsetData(D14_SW, subset.raw = T)
save(D14_SW, file = "D14_SW.Robj")

CPD_D14_SW = as.data.frame.array(D14_SW@data)
CPD_D14_SW <- as.matrix(CPD_D14_SW)
write.csv(CPD_D14_SW, file = "CPD_D14_SW.csv")

CPDmeta_D14_SW <- as.data.frame.array(D14_SW@meta.data)
CPDmeta_D14_SW <- as.matrix(CPDmeta_D14_SW)
write.csv(CPDmeta_D14_SW, file = "CPDmeta_D14_SW.csv")

### Day 14 LWCenter tdTom+ve and tdTom-ve
D14_LWC <- FilterCells(wound, cells.use = c(sample_1, sample_2), subset.names = "nGene", low.thresholds = 400, high.thresholds = 5500)
D14_LWC <- NormalizeData(object = D14_LWC, normalization.method = "LogNormalize", scale.factor = 10000)
D14_LWC <- FindVariableGenes(D14_LWC, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = D14_LWC@var.genes) #2548
D14_LWC <- ScaleData(D14_LWC, vars.to.regress = c("nUMI", "percent.mito"))
D14_LWC <- RunPCA(object = D14_LWC, pc.genes = D14_LWC@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 50)
PrintPCA(object = D14_LWC, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(D14_LWC, pcs.use = 1:2)
PCAPlot(D14_LWC, dim.1 =1, dim.2 = 2)
PCHeatmap(D14_LWC, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCElbowPlot(D14_LWC, num.pc = 50) # 25 PCs
D14_LWC <- FindClusters(D14_LWC, reduction.type = "pca", dims.use = 1:25, resolution = 0.5, print.output = 0, save.SNN = TRUE)
D14_LWC <- StashIdent(D14_LWC, save.name = "D14_LWC_0.5")
D14_LWC <- RunTSNE(D14_LWC, dims.use = 1:25, do.fast = T)
TSNEPlot(D14_LWC, do.label = TRUE, pt.size = 0.8)
#pause
## 1 9 12 6 Epi; 2 4(Langerhan) macrophases; 14 8 0 3 5 fibros; 10 15 endothelial cells; 13 Schwann; 
## 11 Tcells; 7 pericytes; SKM Progenitors; 16 Neutrophils
FeaturePlot(D14_LWC, features.plot = c("Krt14"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Krt5"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Krt28"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Cdh1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Cd68"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Cd207"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("H2-M5"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Pdgfra"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Dpt"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Pecam1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Sox10"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Cd3g"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Acta2"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Notch3"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Pecam1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Cox6a2"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Mitf"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Ngp"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Ly75"), reduction.use = "tsne", cols.use = c("grey", "blue"))

current <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
new <- c("Fibroblasts",
         "Epithelial cells",
         "Macrophases",
         "Fibroblasts",
         "Langerhans",
         "Fibroblasts",
         "Epithelial cells",
         "Pericytes",
         "Fibroblasts",
         "Epithelial cells",
         "Endothelial Cells",
         "T cells",
         "Epithelial cells",
         "Schwann Cells",
         "Fibroblasts",
         "Endothelial Cells",
         "Neutrophils")
D14_LWC@ident <- plyr::mapvalues(D14_LWC@ident, from = current, to = new)
TSNEPlot(D14_LWC, do.label = TRUE, pt.size = 0.5)
D14_LWC <- StashIdent(D14_LWC, save.name = "Cell_names")
D14_LWC <- SubsetData(D14_LWC, subset.raw = T)
save(D14_LWC, file = "D14_LWC.Robj")

CPD_D14_LWC = as.data.frame.array(D14_LWC@data)
CPD_D14_LWC <- as.matrix(CPD_D14_LWC)
write.csv(CPD_D14_LWC, file = "CPD_D14_LWC.csv")

CPDmeta_D14_LWC <- as.data.frame.array(D14_LWC@meta.data)
CPDmeta_D14_LWC <- as.matrix(CPDmeta_D14_LWC)
write.csv(CPDmeta_D14_LWC, file = "CPDmeta_D14_LWC.csv")

### Day 14 LWPheriphy tdTom+ve and tdTom-ve
load(file = "wound.Robj")
D14_LWP <- FilterCells(wound, cells.use = c(sample_3, sample_4), subset.names = "nGene", low.thresholds = 400, high.thresholds = 5500)
D14_LWP <- NormalizeData(object = D14_LWP, normalization.method = "LogNormalize", scale.factor = 10000)
D14_LWP <- FindVariableGenes(D14_LWP, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = D14_LWP@var.genes) #2548
D14_LWP <- ScaleData(D14_LWP, vars.to.regress = c("nUMI", "percent.mito"))
D14_LWP <- RunPCA(object = D14_LWP, pc.genes = D14_LWP@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 50)
PrintPCA(object = D14_LWP, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(D14_LWP, pcs.use = 1:2)
PCAPlot(D14_LWP, dim.1 =1, dim.2 = 2)
PCHeatmap(D14_LWP, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCElbowPlot(D14_LWP, num.pc = 50) # 25 PCs
D14_LWP <- FindClusters(D14_LWP, reduction.type = "pca", dims.use = 1:25, resolution = 0.5, print.output = 0, save.SNN = TRUE)
D14_LWP <- StashIdent(D14_LWP, save.name = "D14_LWP_0.5")
D14_LWP <- RunTSNE(D14_LWP, dims.use = 1:25, do.fast = T)
TSNEPlot(D14_LWP, do.label = TRUE, pt.size = 0.8)
save(D14_LWP, file = "D14_LWP.Robj")
#pause
## 1 9 12 6 Epi; 2 4(Langerhan) macrophases; 14 8 0 3 5 fibros; 10 15 endothelial cells; 13 Schwann; 
## 11 Tcells; 7 pericytes; SKM Progenitors; 16 Neutrophils
FeaturePlot(D14_LWC, features.plot = c("Krt14"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Krt5"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Krt28"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Cdh1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Cd68"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Cd207"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("H2-M5"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Pdgfra"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Dpt"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Pecam1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Sox10"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Cd3g"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Acta2"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Notch3"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Pecam1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Cox6a2"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Mitf"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Ngp"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_LWC, features.plot = c("Ly75"), reduction.use = "tsne", cols.use = c("grey", "blue"))

current <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
new <- c("Fibroblasts",
         "Epithelial cells",
         "Macrophases",
         "Fibroblasts",
         "Langerhans",
         "Fibroblasts",
         "Epithelial cells",
         "Pericytes",
         "Fibroblasts",
         "Epithelial cells",
         "Endothelial Cells",
         "T cells",
         "Epithelial cells",
         "Schwann Cells",
         "Fibroblasts",
         "Endothelial Cells",
         "Neutrophils")
D14_LWC@ident <- plyr::mapvalues(D14_LWC@ident, from = current, to = new)
TSNEPlot(D14_LWC, do.label = TRUE, pt.size = 0.5)
D14_LWC <- StashIdent(D14_LWC, save.name = "Cell_names")
D14_LWC <- SubsetData(D14_LWC, subset.raw = T)
save(D14_LWC, file = "D14_LWC.Robj")

CPD_D14_LWC = as.data.frame.array(D14_LWC@data)
CPD_D14_LWC <- as.matrix(CPD_D14_LWC)
write.csv(CPD_D14_LWC, file = "CPD_D14_LWC.csv")

CPDmeta_D14_LWC <- as.data.frame.array(D14_LWC@meta.data)
CPDmeta_D14_LWC <- as.matrix(CPDmeta_D14_LWC)
write.csv(CPDmeta_D14_LWC, file = "CPDmeta_D14_LWC.csv")


### Day 14 LW tdTom+ve and tdTom-ve

#### Examing Output from Force Cells:
force_cell_D14 <- Read10X(data.dir = "/Users/SarthakSinha/Desktop/Negative Cell Analysis/D14_LWneg_Cells_FC_5000cells/Agg/filtered_gene_bc_matrices_mex/mm10")
dense.size <- object.size(x = as.matrix(x = force_cell_D14))
sparse.size <- object.size(x = force_cell_D14)
dense.size/sparse.size
D14_neg_5000 <- CreateSeuratObject(raw.data = force_cell_D14, min.cells = 2, min.genes = 150, 
                           project = "FC_5000")
mito.genes <- grep(pattern = "^mt-", x = rownames(x = D14_neg_5000@data), value = TRUE)
percent.mito <- Matrix::colSums(D14_neg_5000@raw.data[mito.genes, ])/Matrix::colSums(D14_neg_5000@raw.data)
D14_neg_5000 <- AddMetaData(object = D14_neg_5000, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = D14_neg_5000, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
par(mfrow = c(1, 2))
GenePlot(object = D14_neg_5000, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = D14_neg_5000, gene1 = "nUMI", gene2 = "nGene")

D14_neg_5000 <- NormalizeData(D14_neg_5000, normalization.method = "LogNormalize", scale.factor = 10000)
D14_neg_5000 <- FindVariableGenes(D14_neg_5000, mean.function = ExpMean, dispersion.function = LogVMR, 
                               x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(D14_neg_5000@var.genes) #1358

D14_neg_5000 <- ScaleData(object = D14_neg_5000, vars.to.regress = c("nUMI", "percent.mito"))

D14_neg_5000 <- RunPCA(object = D14_neg_5000, pc.genes = D14_neg_5000@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PrintPCA(D14_neg_5000, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(D14_neg_5000, pcs.use = 1:2)
PCAPlot(D14_neg_5000, dim.1 =1, dim.2 = 2)
PCHeatmap(D14_neg_5000, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
D14_neg_5000 <- RunPCA(D14_neg_5000, pc.genes = D14_neg_5000@var.genes, pcs.compute = 200, do.print = T, pcs.print = 5, genes.print = 5)


### current
PCHeatmap(D14_neg_5000, pc.use = 1:30, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
PCElbowPlot(D14_neg_5000, num.pc = 50) # 15 PCs

D14_neg_5000 <- FindClusters(D14_neg_5000, reduction.type = "pca", dims.use = 1:26, 
                          resolution = 0.6, print.output = 0, save.SNN = TRUE)
D14_neg_5000 <- StashIdent(D14_neg_5000, save.name = "D14_neg_5000_ClusterRes_0.6")
D14_neg_5000 <- RunTSNE(D14_neg_5000, dims.use = 1:26, do.fast = T)
TSNEPlot(D14_neg_5000, do.label = TRUE, pt.size = 0.8)


FeaturePlot(D14_neg_5000, features.plot = c("Ly6g"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_neg_5000, features.plot = c("Ngp"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(D14_neg_5000, features.plot = c("Csf3r"), reduction.use = "tsne", cols.use = c("grey", "blue")) ## Cd114 gene
FeaturePlot(D14_neg_5000, features.plot = c("Cst7"), reduction.use = "tsne", cols.use = c("grey", "blue")) ## Cd114 gene
FeaturePlot(D14_neg_5000, features.plot = c("Elane"), reduction.use = "tsne", cols.use = c("grey", "blue")) ## Cd114 gene
FeaturePlot(D14_neg_5000, features.plot = c("Wnt3"), reduction.use = "tsne", cols.use = c("grey", "blue")) ## Cd114 gene
FeaturePlot(D14_neg_5000, features.plot = c("Cd68"), reduction.use = "tsne", cols.use = c("grey", "blue")) ## Cd114 gene
FeaturePlot(D14_neg_5000, features.plot = c("Sox10"), reduction.use = "tsne", cols.use = c("grey", "blue")) ## Cd114 gene


FeaturePlot(D14_neg_5000, features.plot = c("mt-Nd1"), reduction.use = "tsne", cols.use = c("grey", "blue")) ## Cd114 gene
FeaturePlot(D14_neg_5000, features.plot = c("Fabp5"), reduction.use = "tsne", cols.use = c("grey", "blue")) ## Cd114 gene


### Cluster 5 marker
clu5_5000.markers <- FindMarkers(D14_neg_5000, ident.1 = c(5), min.pct = 0.25, test.use = "negbinom")
clu5_5000.markers <- clu5_5000.markers[order(clu5_5000.markers$avg_logFC),]
write.csv(clu5_5000.markers,file = "clu5_5000.markers.csv")



Center <- grep("-1", colnames(D14_neg_5000@data), value = T); length(Center)
Pheriphery <- grep("-2", colnames(D14_neg_5000@data), value = T); length(Pheriphery)

or_tSNE <- as.data.frame(GetCellEmbeddings(D14_neg_5000, reduction.type = "tsne"))
or_tSNE[Center, "Sample"] <- "Center_cells"
or_tSNE[Pheriphery, "Sample"] <- "Pheriphery_cells"
ggplot(or_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample), size=1.20) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()


jpeg(file = "Microglia_state1_1.jpeg", width = 20, height = 15, units = "cm", res = 500)
ggplot(or_tSNE, aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(color = Sample), size=1.20) + xlab("tSNE 1") + ylab("tSNE 2") + theme_classic()
dev.off()


save(D14_neg_5000, file = "D14_neg_5000.Robj")





#################Uninjuired Charecterization#######################
load(file = "UI_tdTom.Robj")
TSNEPlot(UI_tdTom, do.label = TRUE, pt.size = 0.8)
UI_tdTom = StashIdent(UI_tdTom, save.name = "CellType_Class1")
View(UI_tdTom@meta.data)
save(UI_tdTom, file = "UI_tdTom.Robj")
UI_tdTom = SetAllIdent(UI_tdTom, id = "res.1.5")

TSNEPlot(UI_tdTom, do.label = TRUE, pt.size = 0.8)
FeaturePlot(UI_tdTom, features.plot = c("Dpp4"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Dlk1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Ly6a"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("tdTomato"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Fabp4"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Pdgfra"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Dpt"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Pecam1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Cdh5"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Acta2"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Notch3"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Crabp1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Col11a1"), reduction.use = "tsne", cols.use = c("grey", "blue"))


# 1, 15 - Endothelial
# 0, 6, 9 - Papillary_Dermal_Fibros
# 5,3,4,10,17,7,2,14 - Reticular_Dermal_Fibros
# 8 - Dermal_Sheath
# 11 - Dermal_Papilla
# 13, 16 - Pericytes
# 12 - Adipocytes


current <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
new <- c("Papillary_Dermal_Fibros", #0
         "Endothelial", #1
         "Reticular_Dermal_Fibros", #2
         "Reticular_Dermal_Fibros", #3
         "Reticular_Dermal_Fibros", #4
         "Reticular_Dermal_Fibros", #5
         "Papillary_Dermal_Fibros", #6 
         "Reticular_Dermal_Fibros", #7
         "Dermal_Sheath", #8
         "Papillary_Dermal_Fibros", #9
         "Reticular_Dermal_Fibros", #10
         "Dermal_Papilla", #11
         "Adipocytes", #12
         "Pericytes", #13
         "Reticular_Dermal_Fibros", #14
         "Endothelial", #15
         "Pericytes", #16
         "Reticular_Dermal_Fibros")
UI_tdTom@ident <- plyr::mapvalues(UI_tdTom@ident, from = current, to = new)
TSNEPlot(UI_tdTom, do.label = TRUE, pt.size = 0.5)
UI_tdTom = StashIdent(UI_tdTom, save.name = "CellType_Final_Classification")

UI_meta = UI_tdTom@meta.data
# save(UI_meta, file = "wound.Robj")
write.csv(UI_meta,file = "UI_meta.csv")

table(UI_tdTom@ident)
cell.num <- table(UI_tdTom@ident)
ClusterLabels = paste("Cluster", names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)



pdf("UI_tdTom_Quantification.pdf")
TSNEPlot(UI_tdTom, do.label = TRUE, pt.size = 0.5)
TSNEPlot(object = UI_tdTom, do.return = T) +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "t-SNE 1",
       y = "t-SNE 2")
FeaturePlot(UI_tdTom, features.plot = c("Dpp4"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Dlk1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Ly6a"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("tdTomato"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Fabp4"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Pdgfra"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Dpt"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Pecam1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Cdh5"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Acta2"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Notch3"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Crabp1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Col11a1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Lrig1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Prdm1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
s
dev.off()

VlnPlot(UI_tdTom, features.plot = c("Ly6a"))

jpeg(file = "UI_tdTom_tsne.jpeg", width = 20, height = 15, units = "cm", res = 500)
TSNEPlot(UI_tdTom, do.label = F, pt.size = 0.85)
dev.off()


library(plotrix)
slices <- c(684, 1371, 356,  133, 96, 89, 94) 
lbls <- c("Papillary_Dermal_Fibros", "Reticular_Dermal_Fibros",
          "Endothelial", "Dermal_Sheath", "Dermal_Papilla", "Adipocytes",
          "Pericytes")

jpeg(file = "UI_tdTom_tsne_pi_chart.jpeg", width = 20, height = 15, units = "cm", res = 500)
pie3D(slices,labels=lbls,explode=0.07,height=0.15,theta=pi/3,start=0,shade=0.6,
      main="UI_tdT_Distribution",
      col=c("#F35E5A","#46A803","#B88612","#18B683",
            "#18A6E4", "#9370FF", "#FA3FCB"))

dev.off()


pie3D(slices,labels=NULL,explode=0.1,height=0.15,theta=pi/3,start=0,shade=0.6,
      main="UI_tdT_Distribution",
      col=c("#F35E5A","#46A803","#B88612","#18B683",
            "#18A6E4", "#9370FF", "#FA3FCB"))











