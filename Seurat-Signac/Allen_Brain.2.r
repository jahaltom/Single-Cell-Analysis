library(Seurat)
library(ggplot2)
library(readr)
library(Signac)
library(scProportionTest)




merged_KA<-readRDS("06122024_combined_juvenile_integrated.rds")


Idents(merged_KA)<-"seurat_clusters"

pdf("DimPlot.06122024_combined_juvenile_integrated.clusters.pdf",width=12,height=13)
DimPlot(merged_KA, reduction = "umap_pca", split.by = "condition", label = TRUE, repel = TRUE)+NoLegend()
dev.off()

prop_test <- sc_utils(merged_KA)

prop_test_treatment <- permutation_test(
    prop_test, cluster_identity = "seurat_clusters",
    sample_1 = "normoxia", sample_2 = "hypoxia",
    sample_identity = "condition"
)

pdf("permutation_plot.06122024_combined_juvenile_integrated.pdf",width=12,height=13)
permutation_plot(prop_test_treatment)
dev.off()



allen_brain_ACA<-readRDS("ACA_10Xv2_seurat.rds")






transfer.anchors <- FindTransferAnchors(
  reference = allen_brain_ACA,
  query = merged_KA,
  query.assay = "RNA",
  reduction = 'pcaproject', reference.reduction = "pca",
  dims = 1:20
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen_brain_ACA$subclass,
  weight.reduction = merged_KA[['harmony']],
  dims = 1:20
)


colnames(predicted.labels)[which(names(predicted.labels) == "predicted.id")]  <- "HPF_subclass"

merged_KA <- AddMetaData(object = merged_KA, metadata = predicted.labels)

pdf("DimPlot.Tranfer.Subclass.pdf",width=12,height=13)
DimPlot(merged_KA, reduction = "umap_pca", split.by = "condition",group.by = "HPF_subclass", label = TRUE, repel = TRUE)+NoLegend()
dev.off()


prop_test <- sc_utils(merged_KA)

prop_test_treatment <- permutation_test(
    prop_test, cluster_identity = "HPF_subclass",
    sample_1 = "normoxia", sample_2 = "hypoxia",
    sample_identity = "condition"
)

pdf("permutation_plot.Transfer.subclass.pdf",width=12,height=13)
permutation_plot(prop_test_treatment)
dev.off()


Idents(merged_KA)<-"HPF_subclass"

#predicted_id_markers<-FindAllMarkers(merged_KA, only.pos = TRUE)

#write.csv(predicted_id_markers, "predicted_subclassACA_markers.csv")




predicted.labels.class <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen_brain_ACA$class,
  weight.reduction = merged_KA[['harmony']],
  dims = 1:20
)

colnames(predicted.labels.class)[which(names(predicted.labels.class) == "predicted.id")]  <- "HPF_class"

merged_KA <- AddMetaData(object = merged_KA, metadata = predicted.labels.class)
pdf("DimPlot.Tranfer.Class.pdf",width=12,height=13)
DimPlot(merged_KA, reduction = "umap_pca", split.by = "condition",group.by = "HPF_class", label = TRUE, repel = TRUE)+NoLegend()
dev.off()






prop_test <- sc_utils(merged_KA)

prop_test_treatment <- permutation_test(
    prop_test, cluster_identity = "HPF_class",
    sample_1 = "normoxia", sample_2 = "hypoxia",
    sample_identity = "condition"
)

pdf("permutation_plot.Transfer.class.pdf",width=12,height=13)
permutation_plot(prop_test_treatment)
dev.off()




Idents(merged_KA)<-"HPF_class"

#predicted_id_markers<-FindAllMarkers(merged_KA, only.pos = TRUE)

#write.csv(predicted_id_markers, "predicted_classACA_id_markers.csv")





saveRDS(merged_KA,"done.rds" )
