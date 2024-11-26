# Quality control

```r
library(Signac)
library(Seurat)
library(dplyr)

juvenile <- readRDS("06122024_combined_juvenile_integrated.allen_brain_projection.rds")
DefaultAssay(juvenile)<-"ATAC"
```

Calulate QC metrics that will be used for filtering.

Mononucleosomal / nucleosome-free ratio.  Low NucleosomeSignal indicates a higher proportion of fragments in nucleosome-free regions (open chromatin, typically active regulatory regions).
```r
juvenile <- NucleosomeSignal(juvenile)
```
Compute TSS enrichment score per cell. THis represents the ratio of fragments centered at the TSS to fragments in TSS-flanking regions. Higher is better. 
```r
juvenile <- TSSEnrichment(object = juvenile)
```


#Total number of fragments in peaks: A measure of cellular sequencing depth / complexity. Cells with very few reads may need to be excluded due to low sequencing depth. Cells with extremely high levels may represent doublets, nuclei clumps, or other artefacts.
#add fraction of reads in peaks. Represents the fraction of all fragments that fall within ATAC-seq peaks. Cells with low values (i.e. <15-20%) often represent low-quality cells or technical artifacts that should be removed. Note that this value can be sensitive to the set of peaks used.
#juvenile$pct_reads_in_peaks <- juvenile$peak_region_fragments / juvenile$passed_filters * 100

Add blacklist ratio. Representing reads which are often associated with technical artifacts.
```r
juvenile$blacklist_ratio <- FractionCountsInRegion(
  object = juvenile, 
  assay = 'ATAC',
  regions = blacklist_hg38_unified
)
```


Note that the last three metrics can be obtained from the output of CellRanger

Visualize

nCount_peaks vs TSS enrichment
```r

png("DensityScatter_QC_TSS_vs_nCount_peaks.png",width=15,height=15,units="in",res=300)
DensityScatter(juvenile, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
```
![alt text](https://github.com/jahaltom/Single-Cell-Analysis/blob/main/Seurat-Signac/scATAC/images/DensityScatter_QC_TSS_vs_nCount_peaks.png)

Distribution of each QC metric
```r
#'pct_reads_in_peaks'
png("VlnPlot_QC.png",width=15,height=15,units="in",res=300)
VlnPlot(
  object = juvenile,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()
```
![alt text](https://github.com/jahaltom/Single-Cell-Analysis/blob/main/Seurat-Signac/scATAC/images/VlnPlot_QC.png)

Subseting Seurat object
```r
juvenile <- subset(
  x = juvenile,
  subset = nCount_peaks >= 1000 &
    nCount_peaks <= 25000 &
   # pct_reads_in_peaks > 40 &
    blacklist_ratio <= 0.05 &
    nucleosome_signal < 10 &
    TSS.enrichment >= 2
)

juvenile
png("DensityScatter_QC_TSS_vs_nCount_peaks.2.png",width=15,height=15,units="in",res=300)
DensityScatter(juvenile, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()



#'pct_reads_in_peaks'
png("VlnPlot_QC.2.png",width=15,height=15,units="in",res=300)
VlnPlot(
  object = juvenile,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

saveRDS(juvenile, "06122024_combined_juvenile_integrated.allen_brain_projection.filtered.rds")

```
# Normalization and linear dimensional reduction
Signac performs term frequency-inverse document frequency (TF-IDF) normalization. This is a two-step normalization procedure, that both normalizes across cells to 
correct for differences in cellular sequencing depth, and across peaks to give higher values to more rare peaks.

```r
juvenile <- RunTFIDF(juvenile)
juvenile <- FindTopFeatures(juvenile, min.cutoff = 'q0')
```
ingular value decomposition (SVD) on the TD-IDF matrix, using the features (peaks) selected above
```r
juvenile <- RunSVD(juvenile)
```

correlation between each LSI component and sequencing depth
```r
png(""LSIvsSeqDepth.png",width=15,height=15,units="in",res=300)
DepthCor(juvenile)
dev.off()
```






#Non-linear dimension reduction and clustering
```r
juvenile <- RunUMAP(object = juvenile, reduction = 'lsi', dims = 2:30)
juvenile <- FindNeighbors(object = juvenile, reduction = 'lsi', dims = 2:30)
juvenile <- FindClusters(object = juvenile, verbose = FALSE, algorithm = 3)

png("Cluster.png",width=15,height=15,units="in",res=300)
DimPlot(object = juvenile, label = TRUE) + NoLegend()
dev.off()
```

# Gene activites

# Differential accessibility and Motif Analysis

```r
df=read.csv("DA/32_Immune.da.csv",row.names=1)

open_hypoxia <- rownames(df[df$avg_log2FC > 3, ])
open_normoxia <- rownames(df[df$avg_log2FC < -3, ])

pdf("VlnDA.pdf",width=25,height=15)
plot1 <- print(VlnPlot(
  object = srt,
  features = open_hypoxia[1],
  pt.size = 0.5,
  idents = c("32_Immune_hypoxia","32_Immune_normoxia")
))
dev.off()

```

```r
pdf("FeatureDA.pdf",width=25,height=15)
plot2 <- print(FeaturePlot(
  object = srt,
  features = open_hypoxia[1],
  pt.size = 0.5
))
dev.off()

```



```r
srt <- SortIdents(srt)

# find DA peaks overlapping gene of interest
regions_highlight <- subsetByOverlaps(StringToGRanges(open_hypoxia), LookupGeneCoords(srt, "Gad1"))
pdf("CoveragePlotDA.GAD1.pdf",width=25,height=15)
CoveragePlot(
  object = srt,
  region = "Gad1",
  region.highlight = regions_highlight,
  extend.upstream = 1000,
  extend.downstream = 1000
)
dev.off()
```
