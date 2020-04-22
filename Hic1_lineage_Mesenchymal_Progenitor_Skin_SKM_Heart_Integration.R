# Goal: Integration of scRNA-seq datasets capturing resident Hic1-lineage progenitors within undamaged heart 
# (GSM4216418), skeletal muscle (GSM2976778), and skin (GSM2910020, this study) using Mutual Nearest Neighbors 
# (MNNs)-based anchoring strategy implemented in Seurat v3 (Stuart et al., 2019).

########################################################################
########## Underhill Lab's skeletal muscle MP for integration ##########
########################################################################
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
Underhill_MP <- Read10X(data.dir = "./Single cell RNA-seq tdTomato MPs quiescent - GEO Download/")
UnderhillMP <- CreateSeuratObject(counts = Underhill_MP, project = "Underhill_MP", min.cells = 3, min.features = 200)
UnderhillMP[["percent.mt"]] <- PercentageFeatureSet(UnderhillMP, pattern = "^mt-")
VlnPlot(UnderhillMP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(UnderhillMP, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(UnderhillMP, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
UnderhillMP <- subset(UnderhillMP, subset = nFeature_RNA > 100 & nFeature_RNA < 3000 & percent.mt < 10) 
#cell filteration parameters based on visual inspection of the QC metrics
UnderhillMP <- NormalizeData(UnderhillMP)
UnderhillMP <- FindVariableFeatures(UnderhillMP, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(UnderhillMP), 10)
plot1 <- VariableFeaturePlot(UnderhillMP)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(UnderhillMP)
UnderhillMP <- ScaleData(UnderhillMP, features = all.genes)
UnderhillMP <- RunPCA(UnderhillMP, features = VariableFeatures(object = UnderhillMP))
save(UnderhillMP, file = "UnderhillMP.Robj")
load(file = "UnderhillMP.Robj")
VizDimLoadings(UnderhillMP, dims = 1:2, reduction = "pca")
DimPlot(UnderhillMP, reduction = "pca")
ElbowPlot(UnderhillMP)
UnderhillMP <- FindNeighbors(UnderhillMP, dims = 1:12)
UnderhillMP <- FindClusters(UnderhillMP, resolution = 0.5)
UnderhillMP <- RunUMAP(UnderhillMP, dims = 1:12)
UnderhillMP <- RunTSNE(UnderhillMP, dims = 1:12)
DimPlot(UnderhillMP, reduction = "umap")
DimPlot(UnderhillMP, reduction = "tsne")
FeaturePlot(UnderhillMP, features = c("Cd34", "Ly6a", "Pdgfra"),reduction = "tsne") #PAN
FeaturePlot(UnderhillMP, features = c("Mki67", "Top2a"),reduction = "tsne") #Polif
FeaturePlot(UnderhillMP, features = c("Col4a1", "Col4a2", "Col6a1", "Col6a2", "Col6a3", "Col15a1"),reduction = "tsne") #FAP1
FeaturePlot(UnderhillMP, features = c("Lum", "Sparcl1", "Podn", "Smoc2", "Mgp", "Bgn"),reduction = "tsne") #FAP1
FeaturePlot(UnderhillMP, features = c("Sfrp4", "Igfbp5", "Sema3c", "Dpp4", "Tgfrb2", "Wnt2"),reduction = "tsne") #FAP2
FeaturePlot(UnderhillMP, features = c("Scx", "Mkx", "Tnmd"),reduction = "tsne") #TenogenicMP
FeaturePlot(UnderhillMP, features = c("Rgs5", "Kcnj8", "Mcam", "Pdgfrb"),reduction = "tsne") #PericyticMP
MP.combined = merge(UI_tdTom, y = UnderhillMP, add.cell.ids = c("SKIN", "SKM"), project = "MP")
View(MP.combined@meta.data)
MP.list <- SplitObject(MP.combined, split.by = "orig.ident")
MP.list <- MP.list[c("Uninjured", "Underhill_MP")]
for (i in 1:length(MP.list)) {
  MP.list[[i]] <- NormalizeData(MP.list[[i]], verbose = FALSE)
  MP.list[[i]] <- FindVariableFeatures(MP.list[[i]], selection.method = "vst", 
                                       nfeatures = 2000, verbose = FALSE)}
reference.list <- MP.list[c("Uninjured", "Underhill_MP")]
MP.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:20)
MP.integrated <- IntegrateData(anchorset = MP.anchors, dims = 1:20)
DefaultAssay(MP.integrated) <- "integrated"
# Ran the standard workflow for visualization and clustering
MP.integrated <- ScaleData(MP.integrated, verbose = FALSE)
save(MP.integrated, file = "MP.integrated.Robj")
load(file = "MP.integrated.Robj")
MP.integrated <- RunPCA(MP.integrated, npcs = 50, verbose = FALSE)
MP.integrated <- RunUMAP(MP.integrated, reduction = "pca", dims = 1:20)
DimPlot(MP.integrated, reduction = "umap", group.by = "orig.ident")

pdf("MP.integrated_Unique_Conserved_Cell_Type_markers.pdf", width = 10, height = 8)
DimPlot(MP.integrated, reduction = "umap", group.by = "orig.ident")
FeaturePlot(MP.integrated, features = c("Scx", "Mkx"),reduction = "umap", min.cutoff = "q5") #TenogenicMP
FeaturePlot(MP.integrated, features = c("Mki67", "Top2a"),reduction = "umap", min.cutoff = "q5") #TenogenicMP
FeaturePlot(MP.integrated, features = c("Pecam1", "Cdh5"),reduction = "umap", min.cutoff = "q5") #TenogenicMP
FeaturePlot(MP.integrated, features = c("Wif1","Rspo3"),reduction = "umap", min.cutoff = "q5") #TenogenicMP
FeaturePlot(MP.integrated, features = c("Col4a1", "Col4a2", "Col6a1", "Col6a2", "Col6a3", "Col15a1"),reduction = "umap", min.cutoff = "q5")
FeaturePlot(MP.integrated, features = c("Lum", "Sparcl1", "Podn", "Smoc2", "Mgp", "Bgn"),reduction = "umap", min.cutoff = "q5")
FeaturePlot(MP.integrated, features = c("Sfrp4", "Igfbp5", "Sema3c", "Dpp4", "Tgfrb2", "Wnt2"),reduction = "umap", min.cutoff = "q5")
FeaturePlot(MP.integrated, features = c("Rgs5", "Kcnj8", "Mcam", "Pdgfrb"), reduction = "umap",  min.cutoff = "q5" )
dev.off()

#### Save Progress as RObj:
save(MP.integrated, file = "MP.integrated.Robj")
save(UnderhillMP, file = "UnderhillMP.Robj")
save(UI_tdTom, file = "UI_tdTom.Robj")

p1 <- DimPlot(MP.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(MP.integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
              repel = TRUE) + NoLegend()
plot_grid(p1, p2)

########################################################################
########## Rossi Lab's heart MP for integration ########################
########################################################################
Fabio_MP <- Read10X(data.dir = "./Single cell Heart Hic1 CreERT2 tdTomato+ cells_undamaged/")
FabioMP <- CreateSeuratObject(counts = Fabio_MP, project = "Fabio_MP", min.cells = 3, min.features = 200)
FabioMP[["percent.mt"]] <- PercentageFeatureSet(FabioMP, pattern = "^mt-")
VlnPlot(FabioMP, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(FabioMP, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(FabioMP, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
FabioMP <- subset(FabioMP, subset = nFeature_RNA > 100 & nFeature_RNA < 3000 & percent.mt < 10)
FabioMP <- NormalizeData(FabioMP)
FabioMP <- FindVariableFeatures(FabioMP, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(FabioMP), 10)
plot1 <- VariableFeaturePlot(FabioMP)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(FabioMP)
FabioMP <- ScaleData(FabioMP, features = all.genes)
FabioMP <- RunPCA(FabioMP, features = VariableFeatures(object = FabioMP))
save(FabioMP, file = "FabioMP.Robj")
load(file = "FabioMP.Robj")
VizDimLoadings(FabioMP, dims = 1:2, reduction = "pca")
DimPlot(FabioMP, reduction = "pca")
ElbowPlot(FabioMP)
FabioMP <- FindNeighbors(FabioMP, dims = 1:15)
FabioMP <- FindClusters(FabioMP, resolution = 0.5)
FabioMP <- RunUMAP(FabioMP, dims = 1:15)
FabioMP <- RunTSNE(FabioMP, dims = 1:15)
DimPlot(FabioMP, reduction = "umap")
DimPlot(FabioMP, reduction = "tsne")
FeaturePlot(FabioMP , features = c("Cd34", "Ly6a", "Pdgfra"),reduction = "tsne") #PAN
FeaturePlot(FabioMP , features = c("Mki67", "Top2a"),reduction = "tsne") #Polif
FeaturePlot(FabioMP , features = c("Col4a1", "Col4a2", "Col6a1", "Col6a2", "Col6a3", "Col15a1"),reduction = "tsne") #FAP1
FeaturePlot(FabioMP , features = c("Lum", "Sparcl1", "Podn", "Smoc2", "Mgp", "Bgn"),reduction = "tsne") #FAP1
FeaturePlot(FabioMP , features = c("Sfrp4", "Igfbp5", "Sema3c", "Dpp4", "Tgfrb2", "Wnt2"),reduction = "tsne") #FAP2
FeaturePlot(FabioMP , features = c("Scx", "Mkx", "Tnmd"),reduction = "tsne") #TenogenicMP
FeaturePlot(FabioMP , features = c("Rgs5", "Kcnj8", "Mcam", "Pdgfrb"),reduction = "tsne") #PericyticMP
MP.SKM.SKIN.HEART.combined = merge(UI_tdTom, y = c(UnderhillMP, FabioMP), add.cell.ids = c("SKIN", "SKM", "HEART"), project = "MP")
View(MP.SKM.SKIN.HEART.combined@meta.data)
MP.SKM.SKIN.HEART.list <- SplitObject(MP.SKM.SKIN.HEART.combined, split.by = "orig.ident")
MP.SKM.SKIN.HEART.list <- MP.SKM.SKIN.HEART.list[c("Uninjured", "Underhill_MP", "Fabio_MP")]
for (i in 1:length(MP.SKM.SKIN.HEART.list)) {
  MP.SKM.SKIN.HEART.list[[i]] <- NormalizeData(MP.SKM.SKIN.HEART.list[[i]], verbose = FALSE)
  MP.SKM.SKIN.HEART.list[[i]] <- FindVariableFeatures(MP.SKM.SKIN.HEART.list[[i]], selection.method = "vst", 
                                                      nfeatures = 2000, verbose = FALSE)}
reference.list <- MP.SKM.SKIN.HEART.list[c("Uninjured", "Underhill_MP", "Fabio_MP")]
MP.SKM.SKIN.HEART.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:20)
MP.SKM.SKIN.HEART.integrated <- IntegrateData(anchorset = MP.SKM.SKIN.HEART.anchors, dims = 1:20)
DefaultAssay(MP.SKM.SKIN.HEART.integrated) <- "integrated"
# Ran the standard workflow for visualization and clustering
MP.SKM.SKIN.HEART.integrated <- ScaleData(MP.SKM.SKIN.HEART.integrated, verbose = FALSE)
save(MP.SKM.SKIN.HEART.integrated, file = "MP.SKM.SKIN.HEART.integrated.Robj")
load(file = "MP.SKM.SKIN.HEART.integrated.Robj")
MP.SKM.SKIN.HEART.integrated <- RunPCA(MP.SKM.SKIN.HEART.integrated, npcs = 50, verbose = FALSE)
MP.SKM.SKIN.HEART.integrated <- RunUMAP(MP.SKM.SKIN.HEART.integrated, reduction = "pca", dims = 1:20)
MP.SKM.SKIN.HEART.integrated <- RunTSNE(MP.SKM.SKIN.HEART.integrated, reduction = "pca", dims = 1:20)
DimPlot(MP.SKM.SKIN.HEART.integrated, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident")
View(MP.SKM.SKIN.HEART.integrated@meta.data)
MP.SKM.SKIN.HEART.integrated <- FindNeighbors(MP.SKM.SKIN.HEART.integrated, dims = 1:15)
MP.SKM.SKIN.HEART.integrated <- FindClusters(MP.SKM.SKIN.HEART.integrated, resolution = 0.5)
DimPlot(MP.SKM.SKIN.HEART.integrated, reduction = "umap", group.by = "orig.ident", label = F)
DimPlot(MP.SKM.SKIN.HEART.integrated, reduction = "umap", group.by = "seurat_clusters", label = T)
DimPlot(MP.SKM.SKIN.HEART.integrated, reduction = "tsne", group.by = "seurat_clusters", label = T)

pdf("MP.SKM.SKIN.HEART.integrated_Signatures.pdf", width = 10, height = 8)
FeaturePlot(MP.SKM.SKIN.HEART.integrated , features = c("Cd34", "Ly6a", "Pdgfra"),reduction = "umap", min.cutoff = "q5") #PAN
FeaturePlot(MP.SKM.SKIN.HEART.integrated , features = c("Mki67", "Top2a"),reduction = "umap", min.cutoff = "q5") #Polif
FeaturePlot(MP.SKM.SKIN.HEART.integrated , features = c("Col4a1", "Col4a2", "Col6a1", "Col6a2", "Col6a3", "Col15a1"),reduction = "umap", min.cutoff = "q5") #FAP1
FeaturePlot(MP.SKM.SKIN.HEART.integrated , features = c("Lum", "Sparcl1", "Podn", "Smoc2", "Mgp", "Bgn"),reduction = "umap", min.cutoff = "q5") #FAP1
FeaturePlot(MP.SKM.SKIN.HEART.integrated , features = c("Sfrp4", "Igfbp5", "Sema3c", "Dpp4", "Tgfrb2", "Wnt2"),reduction = "umap", min.cutoff = "q5") #FAP2
FeaturePlot(MP.SKM.SKIN.HEART.integrated , features = c("Scx", "Mkx", "Tnmd"),reduction = "umap", min.cutoff = "q5") #TenogenicMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Rgs5", "Kcnj8", "Mcam", "Pdgfrb"),reduction = "umap", min.cutoff = "q5") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Rgs5", "Kcnj8", "Mcam", "Pdgfrb"),reduction = "umap", min.cutoff = "q5") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Rgs5", "Kcnj8", "Mcam", "Pdgfrb"),reduction = "umap", min.cutoff = "q5") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Pecam1", "Cdh5"),reduction = "umap", min.cutoff = "q5") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Col11a1", "Rspo3", "Wif1"),reduction = "umap", min.cutoff = "q5") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Tnmd", "Scx", "Kera", "Col22a1", "Mkx", "Rspo3"),reduction = "umap", min.cutoff = "q5") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Crabp1"),reduction = "umap", min.cutoff = "q5") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Cxcl2"),reduction = "umap", min.cutoff = "q5") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Mbp", "Fstl1", "Mpz", "Ngfr", "Thbs1"),reduction = "umap", min.cutoff = "q2") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("tdTomato"),reduction = "umap", min.cutoff = "q2") #PericyticMP
dev.off()

View(MP.SKM.SKIN.HEART.integrated@meta.data)
VlnPlot(MP.SKM.SKIN.HEART.integrated, features = c("tdTomato"),group.by = "orig.ident") #PericyticMP
VlnPlot(MP.SKM.SKIN.HEART.integrated, features = c("Crabp1", "Itga8"),group.by = "seurat_clusters") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Crabp1", "Itga8", "Ednrb", "Acta2"), min.cutoff = "q1") #PericyticMP
DimPlot(MP.SKM.SKIN.HEART.integrated, reduction = "umap", group.by = "seurat_clusters", label = T)
DimPlot(MP.SKM.SKIN.HEART.integrated, reduction = "umap", group.by = "orig.ident", label = T)
DimPlot(MP.SKM.SKIN.HEART.integrated, reduction = "umap", group.by = "", label = T)
View(MP.SKM.SKIN.HEART.integrated@meta.data)

FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Pecam1", "Flt1", "Lyve1"),reduction = "umap", min.cutoff = "q1") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("tdTomato"),reduction = "umap", min.cutoff = "q1") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Hic1"),reduction = "umap", min.cutoff = "q1") #PericyticMP
VlnPlot(MP.SKM.SKIN.HEART.integrated, features = c("tdTomato"), group.by = "orig.ident")
MP.SKM.SKIN.HEART.integrated_cluster14.markers <- FindMarkers(MP.SKM.SKIN.HEART.integrated, ident.1 = 14, min.pct = 0.25)
write.csv(MP.SKM.SKIN.HEART.integrated_cluster14.markers , "MP.SKM.SKIN.HEART.integrated_cluster14.markers.csv")
DimPlot(MP.SKM.SKIN.HEART.integrated, reduction = "umap", group.by = "orig.ident", label = T)
DimPlot(MP.SKM.SKIN.HEART.integrated, reduction = "umap", group.by = "seurat_clusters", label = T)
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Apod"),reduction = "umap", min.cutoff = "q2") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Gpx3", "Pdgfra"),reduction = "umap", min.cutoff = "q2") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Gpx3", "Pdgfra"),reduction = "umap", min.cutoff = "q2") #PericyticMP

### Trying to figure out Cluster 14 Ident:
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Lgals7"),reduction = "umap", min.cutoff = "q2") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Cxcl2"),reduction = "umap", min.cutoff = "q2") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Apod"),reduction = "umap", min.cutoff = "q2") #PericyticMP
load(file = "MP.SKM.SKIN.HEART.integrated.Robj")
DimPlot(MP.SKM.SKIN.HEART.integrated, reduction = "umap", group.by = "seurat_clusters", 
        label = T)
View(MP.SKM.SKIN.HEART.integrated@meta.data)
Idents(MP.SKM.SKIN.HEART.integrated) <- "seurat_clusters"
DimPlot(MP.SKM.SKIN.HEART.integrated, reduction = "umap", group.by = "orig.ident", 
        label = T)
DimPlot(MP.SKM.SKIN.HEART.integrated, reduction = "umap", 
        label = T)
MP.SKM.SKIN.HEART.integrated$annotation_1 = MP.SKM.SKIN.HEART.integrated@active.ident
new.cluster.ids <- c("FAP2_MP2",
                     "FAP1_MP1",
                     "FAP2_MP2",
                     "FAP1_MP1",
                     "FAP1_MP1",
                     "FAP1_MP1",
                     "HFassosiated_MP",
                     "Pericytic_MP",
                     "Endothelial_Cells",
                     "FAP1_MP1",
                     "Pericytic_MP",
                     "Polif_MP",
                     "Endothelial_Cells",
                     "HFassosiated_MP",
                     "UK")
names(new.cluster.ids) <- levels(MP.SKM.SKIN.HEART.integrated)
MP.SKM.SKIN.HEART.integrated <- RenameIdents(MP.SKM.SKIN.HEART.integrated, new.cluster.ids)
View(MP.SKM.SKIN.HEART.integrated@meta.data)
save(MP.SKM.SKIN.HEART.integrated, file = "MP.SKM.SKIN.HEART.integrated.Robj")
write.csv(MP.SKM.SKIN.HEART.integrated@meta.data, file = "MP.SKM.SKIN.HEART.integrated@meta.data.csv")

jpeg(file = "MP.SKM.SKIN.HEART.integrated_umap_nonLabel_origIdent.jpeg", width = 20, height = 15, units = "cm", res = 500)
DimPlot(MP.SKM.SKIN.HEART.integrated, reduction = "umap", label = F, group.by = "orig.ident")
dev.off()

jpeg(file = "MP.SKM.SKIN.HEART.integrated_umap_nonLabel_splitbyOrigIdent.jpeg", width = 45, height = 15, units = "cm", res = 500)
DimPlot(MP.SKM.SKIN.HEART.integrated, reduction = "umap", label = F, split.by = "orig.ident", pt.size = 1)
dev.off()

VlnPlot(MP.SKM.SKIN.HEART.integrated, features = c("tdTomato"), group.by = "orig.ident")
VlnPlot(MP.SKM.SKIN.HEART.integrated, features = c("Hic1"), group.by = "orig.ident")
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("Hic1"),reduction = "umap", min.cutoff = "q2", split.by = "orig.ident") #PericyticMP
FeaturePlot(MP.SKM.SKIN.HEART.integrated, features = c("tdTomato"),reduction = "umap", min.cutoff = "q2", split.by = "orig.ident") #PericyticMP

load(file = "MP.SKM.SKIN.HEART.integrated.Robj")
jpeg(file = "FeaturePlot_tdT_Podn_Rgs5_Dpp4_Wnt2_qchanged_pecam1_Dlk1.jpeg", width = 20, height = 20, units = "cm", res = 500)
FeaturePlot(MP.SKM.SKIN.HEART.integrated , 
            features = c("tdTomato", "Podn", "Rgs5", "Dpp4", "Dlk1", "Pecam1"),
            reduction = "umap", 
            min.cutoff = "q2", max.cutoff = "q95") #FAP2
dev.off()