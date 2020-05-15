
############### FINAL:
save(coembed, file = "coembed_updated.Robj")
load(file = "coembed_updated.Robj")

p1 <- DimPlot(coembed, group.by = "tech", reduction = "umap")
p2 <- DimPlot(coembed, label = TRUE, repel = TRUE)
CombinePlots(list(p1, p2))


DimPlot(coembed, split.by = "tech", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()

View(coembed@meta.data)




trial = FindMarkers(coembed,
            ident.1 = c(0),
            min.pct = 0.25)


View(peaks)






View(trial)

## Extra:









## Extra:
sub.fibro_atac <- FindNeighbors(sub.fibro_atac, dims = 1:10)
sub.fibro_atac <- FindClusters(sub.fibro_atac, resolution = 0.2)
sub.fibro_atac <- RunUMAP(sub.fibro_atac, reduction = "umap", dims = 1:10)
DimPlot(sub.fibro_atac, reduction = "umap")

transfer.anchors <- FindTransferAnchors(reference = fibro_center_rna, query = sub.fibro_atac, features = VariableFeatures(object = fibro_center_rna), 
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = fibro_center_rna$celltype, 
                                     weight.reduction = sub.fibro_atac[["lsi"]])

sub.fibro_atac <- AddMetaData(sub.fibro_atac, metadata = celltype.predictions)


DimPlot(sub.fibro_atac, reduction = "umap", group.by = "predicted.id")
