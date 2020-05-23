
load(file = "UI_tdTom.Robj")
TSNEPlot(UI_tdTom, do.label = TRUE, pt.size = 0.8)

### Dpp4 (CD26) = Papillary Dermal Fibroblast Marker
FeaturePlot(UI_tdTom, features.plot = c("Dpp4"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Pecam1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Dlk1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Ly6a"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Pparg"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Thy1"), reduction.use = "tsne", cols.use = c("grey", "blue"))
FeaturePlot(UI_tdTom, features.plot = c("Ncam2"), reduction.use = "tsne", cols.use = c("grey", "blue"))


View(UI_tdTom@meta.data)
UI_tdTom <- Seurat::StashIdent(UI_tdTom, save.name = "detailed.cluster.ident")
UI_tdTom <- SetAllIdent(UI_tdTom, id = "res.0.1")
TSNEPlot(UI_tdTom, do.label = TRUE, pt.size = 0.8)

table(UI_tdTom@meta.data$detailed.cluster.ident, UI_tdTom@meta.data$orig.ident)

####### Synergy Analysis:
library(Seurat)
load(file = "wound_tdT.Robj")

jpeg(file = "TSNEPlot_wound_tdT.jpeg", width = 20, height = 15, units = "cm", res = 500)
TSNEPlot(wound_tdT, do.label = TRUE, pt.size = 0.4)
dev.off()

jpeg(file = "FeaturePlot_wound_tdT_Dpt.jpeg", width = 20, height = 15, units = "cm", res = 500)
FeaturePlot(wound_tdT, features.plot = c("Dpt"), reduction.use = "tsne", cols.use = c("yellow", "red"),min.cutoff = "q10", max.cutoff = "q95", pt.size = 0.7)
dev.off()

jpeg(file = "FeaturePlot_wound_tdT_Pdgfra.jpeg", width = 20, height = 15, units = "cm", res = 500)
FeaturePlot(wound_tdT, features.plot = c("Pdgfra"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10", max.cutoff = "q95", pt.size = 0.7)
dev.off()

##### Female Cells
jpeg(file = "FeaturePlot_wound_tdT_Xist.jpeg", width = 20, height = 15, units = "cm", res = 500)
FeaturePlot(wound_tdT, features.plot = c("Xist"), reduction.use = "tsne", cols.use = c("yellow", "red"), max.cutoff = "q95", pt.size = 0.7)
dev.off()
##### Male Cells
jpeg(file = "FeaturePlot_wound_tdT_Eif2s3y.jpeg", width = 20, height = 15, units = "cm", res = 500)
FeaturePlot(wound_tdT, features.plot = c("Eif1s3y"), reduction.use = "tsne", cols.use = c("yellow", "red"), pt.size = 0.7)
dev.off()


jpeg(file = "FeaturePlot_wound_tdT_tdTomato.jpeg", width = 20, height = 15, units = "cm", res = 500)
FeaturePlot(wound_tdT, features.plot = c("tdTomato"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10", max.cutoff = "q95", pt.size = 0.7)
dev.off()

jpeg(file = "FeaturePlot_wound_tdT_Gli1.jpeg", width = 20, height = 15, units = "cm", res = 500)
FeaturePlot(wound_tdT, features.plot = c("Gli1"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10")
dev.off()

jpeg(file = "FeaturePlot_wound_tdT_Pecam1.jpeg", width = 20, height = 15, units = "cm", res = 500)
FeaturePlot(wound_tdT, features.plot = c("Pecam1"), reduction.use = "tsne", cols.use = c("yellow", "red"),min.cutoff = "q10", max.cutoff = "q95", pt.size = 0.7)
dev.off()

jpeg(file = "FeaturePlot_wound_tdT_Pecam1.jpeg", width = 20, height = 15, units = "cm", res = 500)
FeaturePlot(wound_tdT, features.plot = c("Pecam1"), reduction.use = "tsne", cols.use = c("yellow", "red"), min.cutoff = "q10")
dev.off()

jpeg(file = "FeaturePlot_wound_tdT_Rgs5.jpeg", width = 20, height = 15, units = "cm", res = 500)
FeaturePlot(wound_tdT, features.plot = c("Rgs5"), reduction.use = "tsne", cols.use = c("yellow", "red"),min.cutoff = "q10", max.cutoff = "q95", pt.size = 0.7)
dev.off()

cluster9.markers <- FindMarkers(wound_tdT, ident.1 = c(9), min.pct = 0.25, test.use = "negbinom")
cluster9.markers <- cluster9.markers[order(cluster9.markers$avg_logFC),]
write.csv(cluster9.markers,file = "cluster9.markers.csv")




###### Displaying the number of cells in all Seurat Clusters as part of the plot legend on a tSNE plot ######
# Calculate number of cells per cluster from object@ident
cell.num <- table(LWD14fibr_subset@ident)

# Add cell number per cluster to cluster labels
ClusterLabels = paste("Cluster", names(cell.num), paste0("(n = ", cell.num, ")"))

# Order legend labels in plot in the same order as 'ClusterLabels'
ClusterBreaks = names(cell.num)

# Plot tSNE with new legend labels for clusters

jpeg(file = "TSNEPlot_LWD14fibr_subset_with_number_of_cells_in_a_cluster.jpeg", width = 20, height = 15, units = "cm", res = 500)
TSNEPlot(LWD14fibr_subset, do.return = T, do.label = TRUE, pt.size = 0.4) +
  scale_colour_discrete(breaks = ClusterBreaks, 
                        labels = ClusterLabels) +
  labs(x = "t-SNE 1",
       y = "t-SNE 2")

dev.off()
