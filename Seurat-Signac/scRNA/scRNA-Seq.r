library(Seurat)
library("pheatmap")

#List of mito and inflamation genes
GeneList=read.csv("23_9.19_Genelist.txt",sep="\t")
geneList=GeneList$HumanGene
geneList=geneList[!is.na(geneList)]
geneList=t(t(geneList)) 
colnames(geneList)=c("HumanGene")
geneList=unique(geneList)



GeneList= GeneList[!duplicated(GeneList$HumanGene),]
ann=GeneList[,c("Major.Pathways")]
ann=t(t(ann))
rownames(ann)=GeneList$HumanGene
colnames(ann)=c("Major.Pathways")
ann=data.frame(ann)



newCols <- colorRampPalette(colorspace::qualitative_hcl(palette = "set3",length(unique(GeneList$Major.Pathways))))
mycolors <- newCols(length(unique(GeneList$Major.Pathways)))
names(mycolors) <- unique(GeneList$Major.Pathways)

ann_color = list("Major.Pathways"=mycolors)






#Read in seurat rds object
srt=readRDS("RNA_AllAssays_Final_Subset.rds")
#Normalize
srtNorm=NormalizeData(srt,assay = 'RNA', normalization.method = 'LogNormalize',scale.factor = median(srt$nCount_RNA))
  
  
  
#Create metadata that has cell type and condition
srtNorm$celltype.stim <- paste(srtNorm$MarkerAnnotations, srtNorm$Status, sep = "_")
#Lable Idents with this metadata
Idents(srtNorm) <- "celltype.stim"

#All cell types in this study
cells=c("B", "CD4", "CD8", "DC", "HSPC", "NK", "PC", "pDC","CD14 Monocytes", "CD16 Monocytes")

#####################Perform DGE
for (c in cells){
      
resultsSFH=FindMarkers(srtNorm, ident.1 = paste(c,"_Severe_Followup",sep=""), ident.2 = paste(c,"_Healthy",sep=""))
resultsSFH$Contrast="Severe_Followup_vs_Healthy"
#Inf and mito genes merge
resultsSFH$HumanGene=rownames(resultsSFH)
infMitoSFH=merge(resultsSFH, geneList, by = 'HumanGene')


resultsSH=FindMarkers(srtNorm, ident.1 = paste(c,"_Severe",sep=""), ident.2 = paste(c,"_Healthy",sep=""))
resultsSH$Contrast="Severe_vs_Healthy"
#Inf and mito genes merge
resultsSH$HumanGene=rownames(resultsSH)
infMitoSH=merge(resultsSH, geneList, by = 'HumanGene')


resultsMFH=FindMarkers(srtNorm, ident.1 = paste(c,"_Mild_Followup",sep=""), ident.2 = paste(c,"_Healthy",sep=""))
resultsMFH$Contrast="Mild_Followup_vs_Healthy"
#Inf and mito genes merge
resultsMFH$HumanGene=rownames(resultsMFH)
infMitoMFH=merge(resultsMFH, geneList, by = 'HumanGene')


resultsMH=FindMarkers(srtNorm, ident.1 = paste(c,"_Mild",sep=""), ident.2 = paste(c,"_Healthy",sep=""))
resultsMH$Contrast="Mild_vs_Healthy"
#Inf and mito genes merge
resultsMH$HumanGene=rownames(resultsMH)
infMitoMH=merge(resultsMH, geneList, by = 'HumanGene')

#All DEGs
results=rbind(resultsSFH,resultsSH,resultsMFH,resultsMH)
write.table(results,paste(c,"_DGE.tsv",sep=""),sep = '\t',row.names = FALSE)
#All Inf/Mito DEGs
results=rbind(infMitoSFH,infMitoSH,infMitoMFH,infMitoMH)
write.table(results,paste(c,"InfMito_DGE.tsv",sep=""),sep = '\t',row.names = FALSE)


#####################Get log2 avg norm counts and make plot(s)
#Lable Idents with cell type
Idents(srtNorm) <- srtNorm$MarkerAnnotations
#Subset to specific cell type
srtNormSub=subset(x = srtNorm, idents = "CD14 Monocytes")

#Extract norm counts
srtNormSub_counts=GetAssayData(srtNormSub, assay = "RNA", slot = "data")


#Transpose to have gene name as columns and sample ID as row. 
srtNormSub_counts=t(as.matrix(srtNormSub_counts))
#Gather metadata for sample ID and condition
md=as.matrix(srtNormSub$Status)

#Merge by sample ID
srtM=merge(srtNormSub_counts, md, by = 'row.names', all = TRUE)
#Drop unnecessary column 
srtM=srtM[ , -which(names(srtM) %in% c("Row.names"))]
#Average by condition 
srtM=aggregate(. ~ V1, srtM, FUN = function(x) mean(as.numeric(as.character(x))))
#Transpose once more
srtM=t(srtM)
colnames(srtM)=(srtM)["V1",]
#Drop unnecessary row
srtM=srtM[!(row.names(srtM) %in% c("V1")),]
#Convert to df and add gene name as column
srtM=as.data.frame(srtM)
srtM$HumanGene=rownames(srtM)


#Severe_Followup_vs_Healthy. Gather mito and inflamation genes that are DE. 
degs = infMitoSFH[infMitoSFH$p_val_adj <=0.05, ]
degs=degs$HumanGene
degs=t(t(degs)) 
colnames(degs)=c("HumanGene")
#Merge with norm counts. 
srtPlot=merge(srtM, degs, by = 'HumanGene')
rownames(srtPlot)=srtPlot$HumanGene

arr=split(srtPlot, ceiling(seq_along(srtPlot$HumanGene)/80))
for (i in 1:length(arr)){
        dfArr=arr[[i]]
        dfArr=dfArr[ , -which(names(dfArr) %in% c("HumanGene"))]
        #Log2 transform
        x=rownames(dfArr)
        dfArr=log2(as.data.frame(sapply(dfArr, as.numeric)))
        rownames(dfArr)=x


        png(paste(c,i,"InfMitoDEGs.png",sep=""),width=10,height=15,units="in",res=300)
        pheatmap(dfArr,cluster_rows=FALSE,cluster_cols=FALSE,cellwidth = 20, cellheight = 10,annotation_row=ann,annotation_colors = ann_color,main=paste("Severe_Followup_vs_Healthy",c,sep=""))
        dev.off()  
           
}}
            
            
            
            
            
             
            
            
            
            
            
            
             
