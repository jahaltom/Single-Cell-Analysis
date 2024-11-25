# Quality control

```r
library(Signac)
library(Seurat)
library(dplyr)

juvenile <- readRDS("06122024_combined_juvenile_integrated.allen_brain_projection.rds")
DefaultAssay(juvenile)<-"ATAC"
```

Mononucleosomal / nucleosome-free ratio.  Low NucleosomeSignal indicates a higher proportion of fragments in nucleosome-free regions (open chromatin, typically active regulatory regions).
```r
juvenile <- NucleosomeSignal(juvenile)
```
Compute TSS enrichment score per cell. ratio of fragments centered at the TSS to fragments in TSS-flanking regions. Higher the better. 
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

nCount_peaks vs TSS enrichment
```r
pdf("DensityScatter_QC_TSS_vs_nCount_peaks.pdf",width=25,height=15)
DensityScatter(juvenile, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()
```

Distribution of each QC metric
```r
#'pct_reads_in_peaks'
pdf("VlnPlot_QC.pdf",width=25,height=15)
VlnPlot(
  object = juvenile,
  features = c('nCount_peaks', 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()
```

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
saveRDS(juvenile, "06122024_combined_juvenile_integrated.allen_brain_projection.filtered.rds")
```
# Normalization and linear dimensional reduction