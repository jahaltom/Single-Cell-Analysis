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

png("DensityScatter_QC_TSS_vs_nCount_peaks.png",width=15,height=15,units="in",res=300)
DensityScatter(juvenile, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()



#'pct_reads_in_peaks'
png("VlnPlot_QC.png",width=15,height=15,units="in",res=300)
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
