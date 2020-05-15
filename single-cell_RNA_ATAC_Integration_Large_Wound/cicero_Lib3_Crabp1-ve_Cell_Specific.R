indata <- Matrix::readMM("/home/ssinha/scATAC-Seq/Lib3_out/filtered_peak_bc_matrix/matrix.mtx")
indata@x[indata@x > 0] <- 1
cellinfo <- read.table("/home/ssinha/scATAC-Seq/Lib3_out/filtered_peak_bc_matrix/barcodes.tsv")
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"
peakinfo <- read.table("/home/ssinha/scATAC-Seq/Lib3_out/filtered_peak_bc_matrix/peaks.bed")
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name
row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)
fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
pd <- methods::new("AnnotatedDataFrame", data = cellinfo)
input_cds <-  suppressWarnings(newCellDataSet(indata,
                                              phenoData = pd,
                                              featureData = fd,
                                              expressionFamily=VGAM::binomialff(),
                                              lowerDetectionLimit=0))
input_cds@expressionFamily@vfamily <- "binomialff"
input_cds <- monocle::detectGenes(input_cds)
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 

## Running Cicero:
#Dimensionality reduction and creation of a Cicero CDS
set.seed(2017)
input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)
input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                             reduction_method = 'tSNE', norm_method = "none")


crabp1_id <- row.names(subset(fData(input_cds), site_name == "chr9_54765653_54766062"))
pdgfra_id <- row.names(subset(fData(input_cds), site_name == "chr5_75149771_75157509"))


cth_2 <- newCellTypeHierarchy()
cth_2 <- addCellType(cth_2, "Pdgfra+Crabp1-fibros", classify_func = function(x) { x[pdgfra_id,] == 1 & x[crabp1_id,] == 0} )
input_cds_2 <- classifyCells(input_cds, cth_2, 0.1)
## Filter out just the ZBTB7b+PBMC cells:
input_cds_2 <- input_cds_2[,pData(input_cds_2)$CellType == "Pdgfra+Crabp1-fibros"]


#length((fData(input_cds_2))
#length(pData(input_cds_2))

agg_cds_2 <- aggregate_nearby_peaks(input_cds_2, distance = 10000)
agg_cds_2 <- detectGenes(agg_cds_2)
agg_cds_2 <- estimateSizeFactors(agg_cds_2)
agg_cds_2 <- estimateDispersions(agg_cds_2) ### Warning - some peaks (aka featrues) might be eliminated as outliners
diff_test_res_2 <- differentialGeneTest(agg_cds_2, fullModelFormulaStr = "~CellType", cores = 20)
ordering_sites_2 <- row.names(subset(diff_test_res_2, qval < 1))

#jpeg(file = "plot_pc_variance_explained_agg_cds.jpeg", width = 40, height = 10, units = "cm", res = 500)
#plot_pc_variance_explained(agg_cds, return_all = F)
#dev.off()

agg_cds_2 <- reduceDimension(agg_cds_2,
                           max_components = 2,
                           norm_method = 'log',
                           num_dim = 4, #### Choose # of dimentions appropriately according to plot
                           reduction_method = 'tSNE',
                           verbose = T,
                           perplexity = 10)

agg_cds_2 <- clusterCells(agg_cds_2, verbose = F)

jpeg(file = "plot_cell_clusters_PDGFRA+_Crabp1-_1.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_cell_clusters(agg_cds_2, color_by = 'as.factor(Cluster)')
dev.off()

clustering_DA_sites_2 <- differentialGeneTest(agg_cds_2, #Takes a few minutes
                                            fullModelFormulaStr = '~CellType')

ordering_sites_2 <- row.names(clustering_DA_sites_2)[order(clustering_DA_sites_2$qval)][1:1000]
agg_cds_2 <- setOrderingFilter(agg_cds_2, ordering_sites_2)

agg_cds_2 <- reduceDimension(agg_cds_2, max_components = 2,
                           residualModelFormulaStr="~as.numeric(num_genes_expressed)",
                           reduction_method = 'DDRTree')
agg_cds_2 <- orderCells(agg_cds_2)


####### SKIP
jpeg(file = "plot_cell_trajectory_state_1PDGFRA+_Crabp1-.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_cell_trajectory(agg_cds_2, color_by = "State")
dev.off()
jpeg(file = "plot_cell_trajectory_Pseudotime_1PDGFRA+_Crabp1-.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_cell_trajectory(agg_cds_2, color_by = "Pseudotime")
dev.off()
####### SKIP

pData(input_cds_2)$Pseudotime <- pData(agg_cds_2)[colnames(input_cds_2),]$Pseudotime
pData(input_cds_2)$State <- pData(agg_cds_2)[colnames(input_cds_2),]$State

input_cds_lin_2 <- input_cds_2[,row.names(subset(pData(input_cds_2), State  != 5))]




####### PDGFRA+ CRABP1 + Fibros
jpeg(file = "plot_accessibility_in_pseudotime_pdgfra_1 - PDGFRA+ CRABP1 + Fibros.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin[c("chr5_75149771_75157509")])
dev.off()

jpeg(file = "plot_accessibility_in_pseudotime_crabp1_accessibility_PDGFRA+_Crabp1+ Fibro.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin[c("chr9_54765653_54766062")])
dev.off()

jpeg(file = "plot_accessibility_in_pseudotime_rspo3promoter_PDGFRA+_Crabp1+ Fibro.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin[c("chr10_29535452_29536158")])
dev.off()

jpeg(file = "plot_accessibility_in_pseudotime_rspo3putativeEnhancer_PDGFRA+_Crabp1+ Fibro.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin[c("chr10_29312723_29313879")])
dev.off()

jpeg(file = "plot_accessibility_in_pseudotime_Cd24a_PDGFRA+_Crabp1+.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin[c("chr10_43577658_43580162")])
dev.off()

jpeg(file = "plot_accessibility_in_pseudotime_Btbd10_promoter_PDGFRA+_Crabp1+.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin[c("chr7_113377514_113378736")])
dev.off()

jpeg(file = "plot_accessibility_in_pseudotime_sox18_PDGFRA+_Crabp1+.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin[c("chr2_181670381_181671868")])
dev.off()



####### PDGFRA+ CRABP1 - Fibros
jpeg(file = "plot_accessibility_in_pseudotime_PDGFRA_accessibility_PDGFRA+_Crabp1-.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin_2[c("chr5_75149771_75157509")])
dev.off()

jpeg(file = "plot_accessibility_in_pseudotime_crabp1_accessibility_PDGFRA+_Crabp1-.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin_2[c("chr9_54765653_54766062")])
dev.off()

jpeg(file = "plot_accessibility_in_pseudotime_rspo3promoter_PDGFRA+_Crabp1-.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin_2[c("chr10_29535452_29536158")])
dev.off()

jpeg(file = "plot_accessibility_in_pseudotime_rspo3putativeEnhancer_PDGFRA+_Crabp1-.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin_2[c("chr10_29312723_29313879")])
dev.off()

jpeg(file = "plot_accessibility_in_pseudotime_Cd24a_PDGFRA+_Crabp1-.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin_2[c("chr10_43577658_43580162")])
dev.off()

jpeg(file = "plot_accessibility_in_pseudotime_Btbd10_promoter_PDGFRA+_Crabp1-.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin_2[c("chr7_113377514_113378736")])
dev.off()

jpeg(file = "plot_accessibility_in_pseudotime_sox18_PDGFRA+_Crabp1-.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin_2[c("chr2_181670381_181671868")])
dev.off()


