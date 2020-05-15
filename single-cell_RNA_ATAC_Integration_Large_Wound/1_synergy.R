library(data.table)
library(SummarizedExperiment)
library(GenomicRanges)
library(diffloop)
library(dplyr)

peaks <- diffloop::bedToGRanges("peaks.bed") # taken from the outs folder

############ Function to get counts ############
getCountsFromFrags <- function(frag_gz_file,
                               peaks_gr,
                               barcodes){
  
  # Make GRanges of fragments that are solid for the cells that we care about
  frags_valid <- data.table::fread(paste0("zcat < ", frag_gz_file)) %>% 
    data.frame() %>% filter(V4 %in% barcodes) %>%  # filter for barcodes in our search set
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)
  
  # Get a denominator, per cell
  denom <- table(GenomicRanges::mcols(frags_valid)$V4)
  barcodes_found <- names(denom)
  
  # Get the overlaps with peaks
  ovPEAK <- GenomicRanges::findOverlaps(peaks_gr, frags_valid)
  
  # Establish a numeric index for the barcodes for sparse matrix purposes
  id <- factor(as.character(GenomicRanges::mcols(frags_valid)$V4), levels = barcodes_found)
  
  # Make sparse matrix with counts with peaks by  unique barcode
  countdf <- data.frame(peaks = S4Vectors::queryHits(ovPEAK),
                        sample = as.numeric(id)[S4Vectors::subjectHits(ovPEAK)]) %>%
    dplyr::group_by(peaks,sample) %>% dplyr::summarise(count = n()) %>% data.matrix()
  
  m <- Matrix::sparseMatrix(i = c(countdf[,1], length(peaks_gr)),
                            j = c(countdf[,2], length(barcodes_found)),
                            x = c(countdf[,3],0))
  colnames(m) <- barcodes_found
  
  # Make a polished colData
  colData <- data.frame(
    sample = barcodes_found,
    depth = as.numeric(denom),
    FRIP = Matrix::colSums(m)/as.numeric(denom)
  )
  # Make sure that the SE can be correctly constructed
  stopifnot(all(colData$sample == colnames(m)))
  
  # Make summarized Experiment
  SE <- SummarizedExperiment::SummarizedExperiment(
    rowRanges = peaks_gr,
    assays = list(counts = m),
    colData = colData
  )
  return(SE)
}
############ END Function to get counts ############

# Import PBMC to a summarized experiment
bc_pub_pbmc <- as.character(read.table("/home/ssinha/scATAC-Seq/PBMC_5K_try1/outs/filtered_peak_bc_matrix/barcodes.tsv")[,1])
pbmc_SE <- getCountsFromFrags("/home/ssinha/scATAC-Seq/PBMC_5K_try1/outs/fragments.tsv.gz", peaks, bc_pub_pbmc)

setwd("/home/ssinha/scATAC-Seq/PBMC_5K_try1/outs")

bamfiles <- c("possorted_bam.bam")
fragment_counts <- getCounts(bamfiles, peaks, 
                             paired =  TRUE, 
                             by_rg = TRUE, 
                             format = "bam", 
                             colData = DataFrame(celltype = c("pbmc")))


# Now compatible with chromVAR... run the pipeline
# filterPeaks; etc.
SE <- addGCBias(se, genome = BSgenome.Hsapiens.UCSC.hg19)

# Match and compute kmers
KmerMatch <- matchKmers(7, se, BSgenome.Hsapiens.UCSC.hg19)
kmerdev <- computeDeviations(se, KmerMatch)
kmerdevTable <- assays(kmerdev)[["deviations"]]

write.table(kmerdevTable, file = "kmer_Zscores.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)