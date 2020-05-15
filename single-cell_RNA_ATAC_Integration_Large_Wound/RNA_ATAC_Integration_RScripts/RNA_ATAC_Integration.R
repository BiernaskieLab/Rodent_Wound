rm(list=ls())

library(Seurat)
library(ggplot2)
library(hdf5r)

Read10X_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    stop("Please install hdf5r to read HDF5 files")
  }
  if (!file.exists(filename)) {
    stop("File not found")
  }
  infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
  genomes <- names(x = infile)
  output <- list()
  if (!infile$attr_exists("PYTABLES_FORMAT_VERSION")) {
    # cellranger version 3
    if (use.names) {
      feature_slot <- 'features/name'
    } else {
      feature_slot <- 'features/id'
    }
  } else {
    if (use.names) {
      feature_slot <- 'gene_names'
    } else {
      feature_slot <- 'genes'
    }
  }
  for (genome in genomes) {
    counts <- infile[[paste0(genome, '/data')]]
    indices <- infile[[paste0(genome, '/indices')]]
    indptr <- infile[[paste0(genome, '/indptr')]]
    shp <- infile[[paste0(genome, '/shape')]]
    features <- infile[[paste0(genome, '/', feature_slot)]][]
    barcodes <- infile[[paste0(genome, '/barcodes')]]
    sparse.mat <- sparseMatrix(
      i = indices[] + 1,
      p = indptr[],
      x = as.numeric(x = counts[]),
      dims = shp[],
      giveCsparse = FALSE
    )
    if (unique.features) {
      features <- make.unique(names = features)
    }
    rownames(x = sparse.mat) <- features
    colnames(x = sparse.mat) <- barcodes[]
    sparse.mat <- as(object = sparse.mat, Class = 'dgCMatrix')
    # Split v3 multimodal
    if (infile$exists(name = paste0(genome, '/features'))) {
      types <- infile[[paste0(genome, '/features/feature_type')]][]
      types.unique <- unique(x = types)
      if (length(x = types.unique) > 1) {
        message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
        sparse.mat <- sapply(
          X = types.unique,
          FUN = function(x) {
            return(sparse.mat[which(x = types == x), ])
          },
          simplify = FALSE,
          USE.NAMES = TRUE
        )
      }
    }
    output[[genome]] <- sparse.mat
  }
  infile$close_all()
  if (length(x = output) == 1) {
    return(output[[genome]])
  } else{
    return(output)
  }
}

peaks <- Read10X_h5("/home/ssinha/scATAC-Seq/Lib3_out/filtered_peak_bc_matrix.h5")

peaks <- Read10X_h5("/Users/SarthakSinha/Desktop/scATAC Analysis/Lib3_out/filtered_peak_bc_matrix.h5")

CreateGeneActivityMatrix <- function(
  peak.matrix,
  annotation.file,
  seq.levels = c(1:22, "X", "Y"),
  include.body = TRUE,
  upstream = 2000,
  downstream = 0,
  verbose = TRUE
) {
  if (!PackageCheck('GenomicRanges', error = FALSE)) {
    stop("Please install GenomicRanges from Bioconductor.")
  }
  if (!PackageCheck('rtracklayer', error = FALSE)) {
    stop("Please install rtracklayer from Bioconductor.")
  }
  
  # convert peak matrix to GRanges object
  peak.df <- rownames(x = peak.matrix)
  peak.df <- do.call(what = rbind, args = strsplit(x = gsub(peak.df, pattern = ":", replacement = "-"), split = "-"))
  peak.df <- as.data.frame(x = peak.df)
  colnames(x = peak.df) <- c("chromosome", 'start', 'end')
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)
  
  # get annotation file, select genes
  gtf <- rtracklayer::import(con = annotation.file)
  gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, pruning.mode = 'coarse')
  GenomeInfoDb::seqlevelsStyle(gtf) <- "UCSC"
  gtf.genes <- gtf[gtf$type == 'gene']
  
  # Extend definition up/downstream
  if (include.body) {
    gtf.body_prom <- Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
  } else {
    gtf.body_prom <- SummarizedExperiment::promoters(x = gtf.genes, upstream = upstream, downstream = downstream)
  }
  gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom)
  keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 0]
  peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
  gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]
  peak.ids$gene.name <- gene.ids$gene_name
  peak.ids <- as.data.frame(x = peak.ids)
  peak.ids$peak <- paste0(peak.ids$seqnames, ":", peak.ids$start, "-", peak.ids$end)
  annotations <- peak.ids[, c('peak', 'gene.name')]
  colnames(x = annotations) <- c('feature', 'new_feature')
  
  # collapse into expression matrix
  peak.matrix <- as(object = peak.matrix, Class = 'matrix')
  all.features <- unique(x = annotations$new_feature)
  
  if (nbrOfWorkers() > 1) {
    mysapply <- future_sapply
  } else {
    mysapply <- ifelse(test = verbose, yes = pbsapply, no = sapply)
  }
  newmat <- mysapply(X = 1:length(x = all.features), FUN = function(x){
    features.use <- annotations[annotations$new_feature == all.features[[x]], ]$feature
    submat <- peak.matrix[features.use, ]
    if (length(x = features.use) > 1) {
      return(Matrix::colSums(x = submat))
    } else {
      return(submat)
    }
  })
  newmat <- t(x = newmat)
  rownames(x = newmat) <- all.features
  colnames(x = newmat) <- colnames(x = peak.matrix)
  return(as(object = newmat, Class = 'dgCMatrix'))
}

activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = "Mus_musculus.GRCm38.96.gtf", 
                                            seq.levels = c(1:20, "X", "Y"), upstream = 2000, verbose = TRUE)




####################################
############## ROUGH: ##############
####################################
# create a gene activity matrix from the peak matrix and GTF, using chromosomes 1:22, X, and Y.
# Peaks that fall within gene bodies, or 2kb upstream of a gene, are considered
activity.matrix <- CreateGeneActivityMatrix(peak.matrix = peaks, annotation.file = "../data/Homo_sapiens.GRCh37.82.gtf", 
                                            seq.levels = c(1:22, "X", "Y"), upstream = 2000, verbose = TRUE)



peaks <- Read10X_h5("/Users/SarthakSinha/Desktop/filtered_peak_bc_matrix.h5")
infile <- H5File$new("filtered_peak_bc_matrix.h5")

View(infile)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("cicero", version = "3.8")

library(cicero)

# read in matrix data using the Matrix package
indata <- Matrix::readMM("filtered_peak_bc_matrix/matrix.mtx") 

# format cell info
cellinfo <- read.table("filtered_peak_bc_matrix/barcodes.tsv")
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"

# format peak info
peakinfo <- read.table("filtered_peak_bc_matrix/peaks.bed")
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# make CDS
fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
pd <- methods::new("AnnotatedDataFrame", data = cellinfo)
input_cds <-  suppressWarnings(newCellDataSet(indata,
                                              phenoData = pd,
                                              featureData = fd,
                                              expressionFamily=VGAM::binomialff(),
                                              lowerDetectionLimit=0))
input_cds@expressionFamily@vfamily <- "binomialff"
input_cds <- monocle::detectGenes(input_cds)

#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 


load(file = "input_cds")


######################## START #######################
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
fibro_center_atac <- RunUMAP(fibro_center_atac, reduction = "lsi", dims = 1:50)

load(file = "LWC_D14_fibr_trajec_new.Robj")
fibro_center_rna = UpdateSeuratObject(LWC_D14_fibr_trajec_new)

fibro_center_rna$tech <- "rna"

p1 <- DimPlot(fibro_center_atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
p2 <- DimPlot(fibro_center_rna, label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))

DimPlot(fibro_center_atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")

transfer.anchors <- FindTransferAnchors(reference = fibro_center_rna, query = fibro_center_atac, features = VariableFeatures(object = fibro_center_rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = fibro_center_rna@active.ident, 
                                     weight.reduction = fibro_center_atac[["lsi"]])

fibro_center_atac <- AddMetaData(fibro_center_atac, metadata = celltype.predictions)

hist(fibro_center_atac$prediction.score.max)
abline(v = 0.5, col = "red")

table(fibro_center_atac$prediction.score.max > 0.5)

fibro_center_atac_filtered <- subset(fibro_center_atac, subset = prediction.score.max > 0.3)
fibro_center_atac_filtered$predicted.id <- factor(fibro_center_atac_filtered$predicted.id, levels = levels(fibro_center_rna))  # to make the colors match

p1 <- DimPlot(fibro_center_atac_filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(fibro_center_rna, label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
  NoLegend()
CombinePlots(plots = list(p1, p2))



# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
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
coembed <- RunUMAP(coembed, dims = 1:30)
coembed$celltype <- ifelse(!is.na(coembed@active.ident), coembed@active.ident, coembed$predicted.id)


p1 <- DimPlot(coembed, group.by = "tech")
p2 <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE)
CombinePlots(list(p1, p2))


DimPlot(coembed)
FeaturePlot(coembed, features = "Pdgfra", max.cutoff = 500)
FeaturePlot(coembed, features = "Rspo3", max.cutoff = 500)
FeaturePlot(coembed, features = "Crabp1", max.cutoff = 500)
FeaturePlot(coembed, features = "S100a4", max.cutoff = 500)
FeaturePlot(coembed, features = "chr12:100199214-100200794", max.cutoff = 500)


coembed$blacklist_region_fragments[is.na(coembed$blacklist_region_fragments)] <- 0
FeaturePlot(coembed, features = "blacklist_region_fragments", max.cutoff = 500)



coembed <- FindVariableFeatures(coembed, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(coembed), 10)

plot1 <- VariableFeaturePlot(coembed)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

DimPlot(coembed, reduction = "umap")
crabp1.fibros.markers <- FindMarkers(coembed, ident.1 = "Crabp1+ve Fibros", min.pct = 0.25)
head(crabp1.fibros.markers , n = 30)


View(fibro_center_atac)
DimPlot(fibro_center_atac, group.by = "predicted.id", label = TRUE, repel = TRUE, cols = c("red", "blue", "green", "yellow")) + ggtitle("scATAC-seq cells") 

Idents(object = fibro_center_atac) <- "predicted.id"

Rspo3_condensate_peaks <- FindMarkers(fibro_center_atac, ident.1 = "Rspo3+ve Condensate", min.pct = 0.25)
head(Rspo3_condensate_peaks, n = 30)

Rspo3_condensate_peaks <- FindMarkers(fibro_center_atac, ident.1 = "Rspo3+ve Condensate", min.pct = 0.25)
head(Rspo3_condensate_peaks, n = 30)

Rspo3_condensate_peaks <- FindMarkers(fibro_center_atac, ident.1 = "Rspo3+ve Condensate", min.pct = 0.25)
head(Rspo3_condensate_peaks, n = 30)

save(fibro_center_atac, file = "fibro_center_atac.Robj")
save(fibro_center_rna, file = "fibro_center_rna.Robj")
save(coembed, file = "coembed.Robj")

load(file = "fibro_center_atac.Robj")
load(file = "fibro_center_rna.Robj")
load(file = "coembed.Robj")
load(file = "D14_LWC.Robj")
load(file = "LW_center_rna.Robj")

###### REdoing Analysis from ALL of LWD14 Center:
D14_LWC = UpdateSeuratObject(D14_LWC)
DimPlot(D14_LWC, reduction = "tsne")
FeaturePlot(D14_LWC, features = "Crabp1", max.cutoff = 500)
FeaturePlot(D14_LWC, features = "Rspo3", max.cutoff = 500)

DimPlot(fibro_center_atac, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)

LW_center_rna = D14_LWC

LW_center_rna$tech <- "rna"

p1 <- DimPlot(fibro_center_atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
p2 <- DimPlot(LW_center_rna, label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))

DimPlot(fibro_center_atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")

transfer.anchors <- FindTransferAnchors(reference = LW_center_rna, query = fibro_center_atac, features = VariableFeatures(object = LW_center_rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = LW_center_rna@active.ident, 
                                     weight.reduction = fibro_center_atac[["lsi"]])

fibro_center_atac <- AddMetaData(fibro_center_atac, metadata = celltype.predictions)

hist(fibro_center_atac$prediction.score.max)
abline(v = 0.5, col = "red")

table(fibro_center_atac$prediction.score.max > 0.3)

fibro_center_atac_filtered <- subset(fibro_center_atac, subset = prediction.score.max > 0.3)
fibro_center_atac_filtered$predicted.id <- factor(fibro_center_atac_filtered$predicted.id, levels = levels(LW_center_rna))  # to make the colors match

p1 <- DimPlot(fibro_center_atac, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(LW_center_rna, label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
  NoLegend()
CombinePlots(plots = list(p1, p2))

save(LW_center_rna, file = "LW_center_rna.Robj")


FeaturePlot(fibro_center_atac, features = "chr9:54764920-54765284", max.cutoff = 500) # Crabp1
FeaturePlot(fibro_center_atac, features = "chr9:54765868-54766308", max.cutoff = 500) # Crabp1
FeaturePlot(fibro_center_atac, features = "chr9:54774196-54774640", max.cutoff = 500) # Crabp1
FeaturePlot(fibro_center_atac, features = "chr9:54765868-54766308", max.cutoff = 500) # Crabp1

FeaturePlot(fibro_center_atac, features = "chr2:181671255-181671675", max.cutoff = 500) # Sox18

FeaturePlot(fibro_center_atac, features = "chr2:181671255-181671675", max.cutoff = 500) # Sox18

FeaturePlot(fibro_center_atac, features = "chr6:30732584-30732952", max.cutoff = 500) # Mest
FeaturePlot(fibro_center_atac, features = "chr6:30732584-30732952", max.cutoff = 500) # Mest
FeaturePlot(fibro_center_atac, features = "", max.cutoff = 500) # Krt14


FeaturePlot(fibro_center_atac, features = "chr6:30732584-30732952", max.cutoff = 500) # Mest


save(fibro_center_atac, file = "fibro_center_atac_clustered.Robj")
load(file = "fibro_center_atac_clustered.Robj")

# PDGFRA:
FeaturePlot(fibro_center_atac, features = "chr5:75155756-75155987", max.cutoff = 500) # Crabp1

# Dpt:
FeaturePlot(fibro_center_atac, features = "chr1:164807542-164807906", max.cutoff = 500) # Dpt



#### Reclustering:

fibro_center_atac <- AddMetaData(fibro_center_atac, metadata = meta)
fibro_center_atac <- subset(fibro_center_atac, subset = nCount_ATAC > 500)
fibro_center_atac$tech <- "atac"
DefaultAssay(fibro_center_atac) <- "ACTIVITY"
fibro_center_atac <- FindVariableFeatures(fibro_center_atac)
fibro_center_atac <- NormalizeData(fibro_center_atac)
fibro_center_atac <- ScaleData(fibro_center_atac)
DefaultAssay(fibro_center_atac) <- "ATAC"
VariableFeatures(fibro_center_atac) <- names(which(Matrix::rowSums(fibro_center_atac) > 100))
fibro_center_atac <- RunLSI(fibro_center_atac, n = 50, scale.max = NULL)
fibro_center_atac <- RunUMAP(fibro_center_atac, reduction = "lsi", dims = 1:18)
fibro_center_atac <- RunTSNE(fibro_center_atac, reduction = "lsi", dims = 1:18)

all.genes <- rownames(fibro_center_atac)
fibro_center_atac <- ScaleData(fibro_center_atac, features = all.genes)
fibro_center_atac <- RunPCA(fibro_center_atac, features = VariableFeatures(object = fibro_center_atac))
fibro_center_atac <- FindNeighbors(fibro_center_atac, dims = 1:10)
fibro_center_atac <- FindClusters(fibro_center_atac, resolution = 0.5)
ElbowPlot(fibro_center_atac)
fibro_center_atac <- RunUMAP(fibro_center_atac, reduction = "lsi", dims = 1:10)
fibro_center_atac <- RunTSNE(fibro_center_atac, reduction = "lsi", dims = 1:10)
DimPlot(fibro_center_atac, reduction = "umap", label = T)
DimPlot(fibro_center_atac, reduction = "tsne")


FeaturePlot(fibro_center_atac, features = "chr1:164807542-164807906", max.cutoff = 500) # Dpt
FeaturePlot(fibro_center_atac, features = "chr1:164807542-164807906", max.cutoff = 500) # Dpt




DimPlot(fibro_center_atac, reduction = "tsne")

FeaturePlot(fibro_center_atac, features = "chr1:164807542-164807906", max.cutoff = 500) # Dpt


####### July 1st - Reclustering after purification of ATAC Matrix:

load(file = "fibro_center_atac.Robj")
load(file = "fibro_center_rna.Robj")


fibro_center_atac <- AddMetaData(fibro_center_atac, metadata = meta)
fibro_center_atac <- subset(fibro_center_atac, subset = nCount_ATAC > 500)
fibro_center_atac$tech <- "atac"
DefaultAssay(fibro_center_atac) <- "ACTIVITY"
fibro_center_atac <- FindVariableFeatures(fibro_center_atac)
fibro_center_atac <- NormalizeData(fibro_center_atac)
fibro_center_atac <- ScaleData(fibro_center_atac)
DefaultAssay(fibro_center_atac) <- "ATAC"
VariableFeatures(fibro_center_atac) <- names(which(Matrix::rowSums(fibro_center_atac) > 100))
fibro_center_atac <- RunLSI(fibro_center_atac, n = 50, scale.max = NULL)
fibro_center_atac <- RunUMAP(fibro_center_atac, reduction = "lsi", dims = 1:18)
fibro_center_atac <- RunTSNE(fibro_center_atac, reduction = "lsi", dims = 1:18)

all.genes <- rownames(fibro_center_atac)
fibro_center_atac <- ScaleData(fibro_center_atac, features = all.genes)
fibro_center_atac <- RunPCA(fibro_center_atac, features = VariableFeatures(object = fibro_center_atac))
fibro_center_atac <- FindNeighbors(fibro_center_atac, dims = 1:10)
fibro_center_atac <- FindClusters(fibro_center_atac, resolution = 0.5)
ElbowPlot(fibro_center_atac)
fibro_center_atac <- RunUMAP(fibro_center_atac, reduction = "lsi", dims = 1:10)
fibro_center_atac <- RunTSNE(fibro_center_atac, reduction = "lsi", dims = 1:10)
DimPlot(fibro_center_atac, reduction = "umap", label = T)
DimPlot(fibro_center_atac, reduction = "tsne")


FeaturePlot(fibro_center_atac, features = "chr1:164807542-164807906", max.cutoff = 500) # Dpt
FeaturePlot(fibro_center_atac, features = "chr9:54764920-54765284", max.cutoff = 500) # Crabp1
FeaturePlot(fibro_center_atac, features = "chr9:54765868-54766308", max.cutoff = 500) # Crabp1
FeaturePlot(fibro_center_atac, features = "chr9:54774196-54774640", max.cutoff = 500) # Crabp1
FeaturePlot(fibro_center_atac, features = "chr9:54765868-54766308", max.cutoff = 500) # Crabp1
FeaturePlot(fibro_center_atac, features = "chr5:75155756-75155987", max.cutoff = 500) # Pdgfra
FeaturePlot(fibro_center_atac, features = "chr2:181671255-181671675", max.cutoff = 500) # Sox18



sub.fibro_atac <- names(fibro_center_atac@active.ident[fibro_center_atac@active.ident %in% c(3,6,2,5,4)])
sub.fibro_atac <- subset(fibro_center_atac, idents = c(3,6,2,5,4), invert = FALSE)

DimPlot(sub.fibro_atac, reduction = "tsne")
DimPlot(sub.fibro_atac, reduction = "umap")


sub.fibro_atac <- ScaleData(sub.fibro_atac, features = all.genes)
sub.fibro_atac <- RunPCA(sub.fibro_atac, features = VariableFeatures(object = sub.fibro_atac))
sub.fibro_atac <- FindNeighbors(sub.fibro_atac, dims = 1:10)
sub.fibro_atac <- FindClusters(sub.fibro_atac, resolution = 0.5)
ElbowPlot(sub.fibro_atac)
sub.fibro_atac <- RunUMAP(sub.fibro_atac, reduction = "lsi", dims = 1:10)
sub.fibro_atac <- RunTSNE(sub.fibro_atac, reduction = "lsi", dims = 1:10)
DimPlot(sub.fibro_atac, reduction = "umap")
DimPlot(sub.fibro_atac, reduction = "tsne")

# PDGFRA:
FeaturePlot(sub.fibro_atac, features = "chr5:75155756-75155987", max.cutoff = 500) # Pdgfra
# FibroMakers
FeaturePlot(sub.fibro_atac, features = "chr1:164807542-164807906", max.cutoff = 500) # Dpt
FeaturePlot(sub.fibro_atac, features = "chr9:54764920-54765284", max.cutoff = 500) # Crabp1
FeaturePlot(sub.fibro_atac, features = "chr9:54765868-54766308", max.cutoff = 500) # Crabp1
FeaturePlot(sub.fibro_atac, features = "chr9:54774196-54774640", max.cutoff = 500) # Crabp1
FeaturePlot(sub.fibro_atac, features = "chr9:54765868-54766308", max.cutoff = 500) # Crabp1

FeaturePlot(sub.fibro_atac, features = "chr2:181671255-181671675", max.cutoff = 500) # Sox18


fibro_center_rna$tech <- "rna"

p1 <- DimPlot(sub.fibro_atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
p2 <- DimPlot(fibro_center_rna, label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))

DimPlot(sub.fibro_atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")

transfer.anchors <- FindTransferAnchors(reference = fibro_center_rna, query = sub.fibro_atac, features = VariableFeatures(object = fibro_center_rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")


fibro_center_rna$celltype = fibro_center_rna@active.ident

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = fibro_center_rna$celltype, 
                                     weight.reduction = sub.fibro_atac[["lsi"]])

sub.fibro_atac <- AddMetaData(sub.fibro_atac, metadata = celltype.predictions)

hist(sub.fibro_atac$prediction.score.max)
abline(v = 0.5, col = "red")

table(sub.fibro_atac$prediction.score.max > 0.5)

sub.fibro_atac_filtered <- subset(sub.fibro_atac, subset = prediction.score.max > 0.5)
sub.fibro_atac_filtered$predicted.id <- factor(sub.fibro_atac_filtered$predicted.id, levels = levels(fibro_center_rna))  # to make the colors match

p1 <- DimPlot(sub.fibro_atac_filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(fibro_center_rna, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
  NoLegend()
CombinePlots(plots = list(p1, p2))



Idents(object = sub.fibro_atac_filtered) <- "predicted.id"
DimPlot(sub.fibro_atac_filtered, label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)



# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(fibro_center_rna)
refdata <- GetAssayData(fibro_center_rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = sub.fibro_atac[["lsi"]])

# this line adds the imputed data matrix to the pbmc.atac object
sub.fibro_atac[["RNA"]] <- imputation
coembed <- merge(x = fibro_center_rna, y = sub.fibro_atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
ElbowPlot(coembed)
coembed <- RunUMAP(coembed, dims = 1:20)
coembed <- RunTSNE(coembed, dims = 1:20)
coembed$celltype <- ifelse(!is.na(coembed$celltype), coembed$celltype, coembed$predicted.id)

DimPlot(coembed, reduction = "umap")
save(coembed, file = "coembed_updated.Robj")


FeaturePlot(coembed, features = "chr1:164807542-164807906", max.cutoff = 500) # Dpt
FeaturePlot(coembed, features = "chr9:54764920-54765284", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "chr9:54765868-54766308", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "chr9:54774196-54774640", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "chr9:54765868-54766308", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "Crabp1", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "chr8:105636056-105638037", max.cutoff = 500) # Crabp1


p1 <- DimPlot(coembed, group.by = "tech")
p2 <- DimPlot(coembed, label = TRUE, repel = TRUE)
CombinePlots(list(p1, p2))

View(coembed)
save(coembed, file = "coembed_updated.Robj")

coembed <- FindNeighbors(coembed, dims = 1:10)
coembed <- FindClusters(coembed, resolution = 0.3)

p1 <- DimPlot(coembed, group.by = "tech", reduction = "umap")
p2 <- DimPlot(coembed, label = TRUE, repel = TRUE)
CombinePlots(list(p1, p2))

DimPlot(coembed, split.by = "tech", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()




### CURRENT!
DimPlot(coembed, label = TRUE, repel = TRUE, reduction = "umap")
FeaturePlot(coembed, features = "Crabp1", max.cutoff = 500, reduction = "umap") # Crabp1
FeaturePlot(coembed, features = "Rspo3", max.cutoff = 500, reduction = "umap") # Crabp1
FeaturePlot(coembed, features = "S100a4", max.cutoff = 500, reduction = "umap") # Crabp1
FeaturePlot(coembed, features = "Cd200", max.cutoff = 500, reduction = "umap") # Crabp1
FeaturePlot(coembed, features = "Sox18", max.cutoff = 500, reduction = "umap") # Crabp1
FeaturePlot(coembed, features = "Rspo3", max.cutoff = 500, reduction = "umap") # Crabp1

FeaturePlot(coembed, features = "Neat1", max.cutoff = 500, reduction = "umap") # Crabp1
FeaturePlot(coembed, features = "Wipi1", max.cutoff = 500, reduction = "umap") # Crabp1


Crabp1_markers <- FindMarkers(coembed, ident.1 = c(0,2), min.pct = 0.25)
FeaturePlot(coembed, features = "Ptn", max.cutoff = 500, reduction = "umap") # Crabp1
FeaturePlot(coembed, features = "chr19:32755639-32758735", reduction = "umap", max.cutoff = "q70") # Crabp1
FeaturePlot(coembed, features = "chr7:25278302-25279802", reduction = "umap", max.cutoff = "q70") # Crabp1
FeaturePlot(coembed, features = "chr7:25278302-25279802", reduction = "umap", max.cutoff = "q70") # Crabp1
FeaturePlot(coembed, features = "chr2:181671255-181671675", reduction = "umap", max.cutoff = "q70") # Crabp1
FeaturePlot(coembed, features = "chr2:181693423-181694196", reduction = "umap", max.cutoff = "q70") # Crabp1





write.csv(Crabp1_markers, file = "Crabp1_markers.csv")


View(fibro_center_atac@assays$ATAC@counts@Dimnames[[1]])









# Peaks
FeaturePlot(coembed, features = "chr17:13760496-13761153", max.cutoff = 500, reduction = "umap") # Crabp1
FeaturePlot(coembed, features = "chr11:51855270-51859253", max.cutoff = 500, reduction = "umap") # Crabp1
FeaturePlot(coembed, features = "chr17:39842932-39848963", max.cutoff = 500, reduction = "umap") # Crabp1
FeaturePlot(coembed, features = "chr5:146260783-146261568", max.cutoff = 500, reduction = "umap") # Crabp1
FeaturePlot(coembed, features = "chr7:73540350-73542306", max.cutoff = 500, reduction = "umap") # Crabp1




write.csv(sub.fibro_atac_filtered.markers_condensate, file = "sub.fibro_atac_filtered.markers_condensate.csv")

Rspo3_condensate <- FindMarkers(coembed, ident.1 = 5, min.pct = 0.25)
head(Rspo3_condensate, n = 30)



########### Redoing integration with all cells in LWD14 Matrix:
load(file = "fibro_center_atac.Robj")
fibro_center_atac # atac matrix
LW_center_rna # rna matrix

LW_center_rna$tech <- "rna"
LW_center_rna$celltype = LW_center_rna@active.ident

p1 <- DimPlot(fibro_center_atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
p2 <- DimPlot(LW_center_rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))

transfer.anchors <- FindTransferAnchors(reference = LW_center_rna, query = fibro_center_atac, features = VariableFeatures(object = LW_center_rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")


celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = LW_center_rna$celltype, 
                                     weight.reduction = pbmc.atac[["lsi"]])

fibro_center_atac <- AddMetaData(fibro_center_atac, metadata = celltype.predictions)


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



FeaturePlot(coembed, features = "chr2:181671255-181671675", max.cutoff = 500) # Sox18



coembed


coembed_filtered = subset(coembed, idents = c(0,1,2,3,4,5))





p1 <- DimPlot(coembed_filtered, group.by = "tech")
p2 <- DimPlot(coembed_filtered, label = TRUE, repel = TRUE)
CombinePlots(list(p1, p2))






coembed_filtered <- ScaleData(coembed_filtered, features = genes.use, do.scale = FALSE)
coembed_filtered <- RunPCA(coembed_filtered, features = genes.use, verbose = FALSE)
ElbowPlot(coembed_filtered)
coembed_filtered <- RunUMAP(coembed_filtered, dims = 1:20)
coembed_filtered$celltype <- ifelse(!is.na(coembed_filtered@active.ident), coembed_filtered@active.ident, coembed_filtered$predicted.id)


p1 <- DimPlot(coembed_filtered, group.by = "tech")
p2 <- DimPlot(coembed_filtered, group.by = "celltype", label = TRUE, repel = TRUE)
CombinePlots(list(p1, p2))


DimPlot(coembed)
DimPlot(coembed, group.by = "tech")
FeaturePlot(coembed_filtered, features = "Pdgfra", max.cutoff = 500)
FeaturePlot(coembed, features = "Rspo3", max.cutoff = 500)
FeaturePlot(coembed, features = "Crabp1", max.cutoff = 500)
FeaturePlot(coembed, features = "S100a4", max.cutoff = 500)
FeaturePlot(coembed, features = "chr12:100199214-100200794", max.cutoff = 500)

coembed_filtered = subset(coembed, idents = c("Cd200+ve DS", "Crabp1+ve Fibros",
                                              "Rspo3+ve Condensate", "S100a4 DSCs"), 
                          invert = FALSE)




DimPlot(coembed_filtered)

# coembed_filtered reculstering 
coembed_filtered <- ScaleData(coembed_filtered)
#coembed_filtered <- RunPCA(coembed_filtered, features = VariableFeatures(object = coembed_filtered))
coembed_filtered <- FindNeighbors(coembed_filtered, dims = 1:10)
coembed_filtered <- FindClusters(coembed_filtered, resolution = 0.5)
ElbowPlot(sub.fibro_atac)
coembed_filtered <- RunUMAP(coembed_filtered, reduction = "lsi", dims = 1:10)
coembed_filtered <- RunTSNE(coembed_filtered, reduction = "lsi", dims = 1:10)
DimPlot(coembed_filtered, reduction = "umap")
DimPlot(coembed_filtered, reduction = "tsne")

DimPlot(sub.fibro_atac, reduction = "umap")


Idents(object = sub.fibro_atac_filtered) <- "predicted.id"
sub.fibro_atac_filtered.markers <- FindMarkers(sub.fibro_atac_filtered, ident.1 = "Crabp1+ve Fibros", min.pct = 0.25)
head(sub.fibro_atac_filtered.markers, n = 100)

sub.fibro_atac_filtered.markers_condensate <- FindMarkers(sub.fibro_atac_filtered, ident.1 = "Rspo3+ve Condensate", min.pct = 0.25)
head(sub.fibro_atac_filtered.markers_condensate, n = 10)



DimPlot(coembed_filtered, reduction = "tsne")
FeaturePlot(coembed_filtered, features = "chr1:164807542-164807906", max.cutoff = 500) # Dpt
FeaturePlot(coembed, features = "chr9:54764920-54765284", max.cutoff = 500) # Crabp1
FeaturePlot(coembed_filtered, features = "chr9:54765868-54766308", max.cutoff = 500) # Crabp1
FeaturePlot(coembed_filtered, features = "chr9:54774196-54774640", max.cutoff = 500) # Crabp1
FeaturePlot(coembed_filtered, features = "chr9:54765868-54766308", max.cutoff = 500) # Crabp1


coembed_filtered <- RunDiffusion(coembed_filtered,genes.use = coembed_filtered@var.genes)






AverageExpression(object, assays = NULL, features = NULL,
                  return.seurat = FALSE, add.ident = NULL, slot = "data",
                  use.scale = FALSE, use.counts = FALSE, verbose = TRUE, ...)





##############################
sub.fibro_atac # atac matrix
fibro_center_rna # rna matrix

p1 <- DimPlot(sub.fibro_atac, reduction = "umap") + NoLegend() + ggtitle("scATAC-seq")
p2 <- DimPlot(fibro_center_rna, group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend() + ggtitle("scRNA-seq")
CombinePlots(plots = list(p1, p2))


transfer.anchors <- FindTransferAnchors(reference = fibro_center_rna, query = sub.fibro_atac, features = VariableFeatures(object = fibro_center_rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = fibro_center_rna@active.ident, 
                                     weight.reduction = sub.fibro_atac[["lsi"]])
sub.fibro_atac <- AddMetaData(sub.fibro_atac, metadata = celltype.predictions)
hist(sub.fibro_atac$prediction.score.max)
abline(v = 0.5, col = "red")

table(sub.fibro_atac$prediction.score.max > 0.5)

sub.fibro_atac.filtered <- subset(sub.fibro_atac, subset = prediction.score.max > 0.5)
sub.fibro_atac.filtered$predicted.id <- factor(sub.fibro_atac.filtered$predicted.id, levels = levels(fibro_center_rna))  # to make the colors match

p1 <- DimPlot(sub.fibro_atac.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
  NoLegend() + scale_colour_hue(drop = FALSE)
p2 <- DimPlot(fibro_center_rna, label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
  NoLegend()
CombinePlots(plots = list(p1, p2))


genes.use <- VariableFeatures(fibro_center_rna)
refdata <- GetAssayData(fibro_center_rna, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = sub.fibro_atac[["lsi"]])

# this line adds the imputed data matrix to the pbmc.atac object
sub.fibro_atac[["RNA"]] <- imputation
coembed <- merge(x = fibro_center_rna, y = sub.fibro_atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:20)
coembed$celltype <- ifelse(!is.na(coembed@active.ident), coembed@active.ident, coembed$predicted.id)


FeaturePlot(coembed, features = "chr1:164807542-164807906", max.cutoff = 500) # Dpt
FeaturePlot(coembed, features = "chr9:54764920-54765284", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "chr9:54765868-54766308", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "chr9:54774196-54774640", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "chr9:54765868-54766308", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "Crabp1", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "chr8:105636056-105638037", max.cutoff = 500) # Crabp1


p1 <- DimPlot(coembed, split.by = "tech", group.by = "celltype")
p2 <- DimPlot(coembed, group.by = "celltype", label = TRUE, repel = TRUE)
CombinePlots(list(p1, p2))


jpeg(file = "DimPlot_split_tech_group_celltype.jpeg", width = 40, height = 15, units = "cm", res = 500)
DimPlot(coembed, split.by = "tech", group.by = "celltype", pt.size = 1.5)
dev.off()


FeaturePlot(coembed, features = "Kcnc4", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "Slc6a17", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "Aspg", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "Kif26a", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "Trpv3", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "Ttll10", max.cutoff = 500) # Crabp1

FeaturePlot(coembed, features = "Trpv3", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "Cdcp1", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "Urah", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "Dsc3", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "Foxn1", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "Lamb3", max.cutoff = 500) # Crabp1

FeaturePlot(coembed, features = "Rxra", max.cutoff = 500) # Crabp1
FeaturePlot(coembed, features = "Pkp3", max.cutoff = 500) # Crabp1




