### Generating COUNTS file and METADATA file to perform interactome analysis between 
### Small Wound Day 8 Cells using CellPhoneDB:

# 1. In SEURAT (using Seurat V2):
library(Seurat)
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


# 2. META DATA PROCESSING:
# Briefly, name the two columns cells, cell_names.

# 3.DATA PROCESSING through gProfiler:
# Open the @data spreadsheet, copy paste the list onto gProfiler Ortho search.
# Download resulting excel.
# From the downloaded excel, select the gene row and human ensemblID rows, and copy it into a new spreadsheet.
# Select both rows, delete duplicates, by only selecting the gene row.
# Copy/paste this list onto @data spreadsheet.

# 4. Upload at: https://www.cellphonedb.org/explore-sc-rna-seq