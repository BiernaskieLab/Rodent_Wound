### Generating COUNTS file and METADATA file to perform interactome analysis between 
### Dermal Condensate and Epithelial Placode Cells using CellPhoneDB:

# 1. In SEURAT (using Seurat V2):
library(Seurat)
TSNEPlot(LW14_plac_condensate, do.label = TRUE, pt.size = 0.8)
LW14_plac_condensate = SubsetData(LW14_plac_condensate, subset.raw = TRUE)
LW14_plac_condensate <- StashIdent(LW14_plac_condensate, save.name = "Cell_names")

CPD = as.data.frame.array(LW14_plac_condensate@data)
CPD <- as.matrix(CPD)
write.csv(CPD, file = "CPD.csv")

CPD_2 = as.data.frame.array(LW14_plac_condensate@meta.data)
CPD_2 <- as.matrix(CPD_2)
write.csv(CPD_2, file = "CPD_2.csv")


# 2. META DATA PROCESSING:
# Briefly, name the two columns cells, cell_names.

# 3.DATA PROCESSING through gProfiler:
# Open the @data spreadsheet, copy paste the list onto gProfiler Ortho search.
# Download resulting excel.
# From the downloaded excel, select the gene row and human ensemblID rows, and copy it into a new spreadsheet.
# Select both rows, delete duplicates, by only selecting the gene row.
# Copy/paste this list onto @data spreadsheet.

# 4. Upload at: https://www.cellphonedb.org/explore-sc-rna-seq