load(file = "input_cds.Robj")
load(file = "annotation_sample.Robj")
input_cds <- annotate_cds_by_site(input_cds, annotation_sample)

A = fData(input_cds)$symbol
B = fData(input_cds)$feature
C = rownames(fData(input_cds))
rownames(fData(input_cds)) = paste(A,B,C, sep = '_')
rownames(input_cds) = rownames(fData(input_cds))



to_be_tested <- row.names(subset(fData(input_cds)))
q_subset <- input_cds[to_be_tested,]
diff_pt <- differentialGeneTest(q_subset, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 16)
diff_pt[,c("site_name", "pval", "qval", "symbol", "feature")]
#save(diff_pt, file = "diff_pt_gene_labelled.Robj")
#load(file = "diff_pt_gene_labelled.Robj")

sig_gene_names <- row.names(subset(diff_pt, qval < 0.5))
jpeg(file = "heatmap_sig_gene_names_300_2_1.jpeg", width = 40, height = 300, units = "cm", res = 50)
plot_pseudotime_heatmap(input_cds[sig_gene_names,],
                        num_clusters = 3,
                        cores = 2,
                        show_rownames = T)
dev.off()

write.csv(diff_pt, file = "diff_pt.csv")


jpeg(file = "heatmap_sig_gene_names_300_2.jpeg", width = 35, height = 360, units = "cm", res = 500)
plot_pseudotime_heatmap(input_cds[sig_gene_names,],
                        num_clusters = 3,
                        cores = 2,
                        show_rownames = T)
dev.off()
jpeg(file = "plot_accessibility_in_pseudotime_rspo3promoter_1.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin[c("chr10_29535452_29536158")])
dev.off()
jpeg(file = "plot_accessibility_in_pseudotime_rspo3putativeEnhancer_1.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin[c("chr10_29312723_29313879")])
dev.off()

jpeg(file = "plot_accessibility_in_pseudotime_sox18.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin[c("chr2_181670381_181671868")])
dev.off()

jpeg(file = "plot_accessibility_in_pseudotime_Hey1.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin[c("chr3_8664429_8665173")])
dev.off()

####### chr10_29535452_29536158 = Rspo3

jpeg(file = "plot_connections_Rspo3_conetworks_1.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_connections(conns, "chr10", 26005452, 31996158, 
                 gene_model = annotation_sample, 
                 coaccess_cutoff = .10, 
                 connection_width = .5, 
                 collapseTranscripts = "longest")
dev.off()

chr2_181670381_181671868

jpeg(file = "plot_connections_Sox18_conetworks_3.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_connections(conns, "chr2", 181000000, 182000000, 
                 gene_model = annotation_sample, 
                 coaccess_cutoff = .20,
                 connection_width = .5,
                 collapseTranscripts = "longest")
dev.off()


jpeg(file = "plot_accessibility_in_pseudotime_Zbtb46.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin[c("chr2_181410434_181411497")])
dev.off()
jpeg(file = "plot_accessibility_in_pseudotime_Zgpat.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin[c("chr2_181381158_181381418")])
dev.off()
jpeg(file = "plot_accessibility_in_pseudotime_Rara.jpeg", width = 20, height = 10, units = "cm", res = 500)
plot_accessibility_in_pseudotime(input_cds_lin[c("chr11_98935890_98945891")])
dev.off()


####### Differential Accessibility Over Time:
library(plotrix)

slices <- c(21127, 72886) 
lbls <- c("Protein Coding", "Non-Protein Coding")

jpeg(file = "Coding_regions_pi_chart.jpeg", width = 20, height = 15, units = "cm", res = 500)
pie3D(slices,labels=NULL,explode=0.07,height=0.15,theta=pi/5,start=1.0,shade=0.3,
      main="Coding_regions_pi_chart_1",
      col=c("Black","Grey"))
dev.off()



slices <- c(1016, 1601,656,110,176,471,428) 
lbls <- c("Antisense", "lincRNA", "TEC", "Bidirectional Promoter lncRNA", "Processed Pseudogene", "Processed Transcript", "Other")

jpeg(file = "NonCoding_regions_pi_chart.jpeg", width = 20, height = 15, units = "cm", res = 500)
pie3D(slices,labels=NULL,explode=0.03,height=0.15,theta=pi/5,start=1.0,shade=0.3,
      main="Coding_regions_pi_chart_1",
      col=c("#F35E5A","#46A803","#B88612","#18B683",
            "#18A6E4", "#9370FF", "#FA3FCB"))
dev.off()



####### Unaltered Genome:
genome_u = read.table(file="/home/ssinha/scATAC-Seq/unfiltered_genome/mm10.chrom.sizes", header = F)
conns_u <- run_cicero(cicero_cds, genome_u)








