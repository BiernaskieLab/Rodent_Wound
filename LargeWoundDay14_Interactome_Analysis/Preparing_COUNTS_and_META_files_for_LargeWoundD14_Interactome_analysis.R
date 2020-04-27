### Generating COUNTS file and METADATA file to perform interactome analysis between 
### Small Wound Day 8 Cells using CellPhoneDB:

# 1. In SEURAT (using Seurat V2):
library(Seurat)
D14_LWC <- FilterCells(wound, cells.use = c(sample_1, sample_2), subset.names = "nGene", low.thresholds = 400, high.thresholds = 5500)
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


# 2. META DATA PROCESSING:
# Briefly, name the two columns cells, cell_names.

# 3.DATA PROCESSING through gProfiler:
# Open the @data spreadsheet, copy paste the list onto gProfiler Ortho search.
# Download resulting excel.
# From the downloaded excel, select the gene row and human ensemblID rows, and copy it into a new spreadsheet.
# Select both rows, delete duplicates, by only selecting the gene row.
# Copy/paste this list onto @data spreadsheet.

# 4. Upload at: https://www.cellphonedb.org/explore-sc-rna-seq