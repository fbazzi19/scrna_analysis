#remember to set the working directory accordingly!

#----library calls----
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggpubr)

#----Colors----
#color pallete to use for figure generation throughout
colorpallete <- c("#DA9C7E", "#66B386", "#E6D939", 
                  "#EDC2C2", "#B1C2FF", "#BCB0CE",
                  "#EB7B70", "#B06CC8", "#626CB2")

#----load in and initialize the data----
decidua.data <- load("./SRA782908_SRS3815606.sparse.RData")
decidua.data<- sm
rm(sm)
#remove ensemblID from rownames (eg, everything after '_ENSG...')
splitnames <- strsplit(row.names(decidua.data), "_ENSG")
splitnames <- matrix(unlist(splitnames),ncol=2,byrow=T)
#make the names unique
names <- make.unique(splitnames[,1])
row.names(decidua.data)<-names
rm(names, splitnames)
#initialize seurat object with raw non normalized data
decidua <- CreateSeuratObject(counts = decidua.data, project = "decidua", min.cells = 3, min.features = 200)
decidua
#look at some genes
decidua.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
#column names
head(colnames(decidua))

#----Quality Control----
#select rows with name starting with MT- (from mitochondrial genome)
grep("^MT-",rownames(decidua),value = TRUE)
#percentage of counts originating from mitochondrial genome
#high mt indicates Poor sample quality, leading to a 
#high fraction of apoptotic or lysing cells
decidua[["percent.mt"]] <- PercentageFeatureSet(decidua, pattern = "^MT-")
#ribosomal protein genes
grep("^RP[LS]",rownames(decidua),value = TRUE)
#percentage of counts originating from ribosomal protein genes
decidua[["percent.rbp"]] <- PercentageFeatureSet(decidua, pattern = "^RP[LS]")
#meta data
head(decidua@meta.data)

# Visualize QC metrics as violin plots - also adding the RPL genes
VlnPlot(decidua, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), ncol = 4)
#without dots:
qc_plotnodots <- VlnPlot(decidua, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), 
                         ncol = 4, pt.size=0)

qc_plotnodots[[1]] <- qc_plotnodots[[1]] +
  scale_fill_manual(values = colorpallete[2])+
  ggtitle("Number of Genes") +
  theme(plot.title = element_text(hjust = 0.5))
qc_plotnodots[[2]] <- qc_plotnodots[[2]] +
  scale_fill_manual(values = colorpallete[4])+
  ggtitle("Number of Reads") +
  theme(plot.title = element_text(hjust = 0.5))
qc_plotnodots[[3]] <- qc_plotnodots[[3]] +
  scale_fill_manual(values = colorpallete[5]) +
  ggtitle("Mitochondrial %") +
  theme(plot.title = element_text(hjust = 0.5))
qc_plotnodots[[4]] <- qc_plotnodots[[4]] +
  scale_fill_manual(values = colorpallete[8])+
  ggtitle("Ribosomal Protein %") +
  theme(plot.title = element_text(hjust = 0.5))

qc_plotnodots &
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.background = element_rect(fill = "#F9FCF3"))

#check if the different parameters are correlated with each other
plot1 <- FeatureScatter(decidua, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(decidua, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot3 <- FeatureScatter(decidua, feature1 = "nCount_RNA", feature2 = "percent.rbp")
plot3

#thresholds for quality control
decidua <- subset(decidua, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)
decidua

#----Normalize Data----
#log of the normalized counts (counts per 10,000 reads)
decidua <- NormalizeData(decidua, normalization.method = "LogNormalize", scale.factor = 10000)
#raw counts
decidua[["RNA"]]$counts
#normalized counts 
decidua[["RNA"]]$data

#genes that have the highest mean expression across cells
apply(decidua[["RNA"]]$data,1,mean) -> gene.expression
sort(gene.expression, decreasing = TRUE) -> gene.expression
head(gene.expression, n=50)
#distribution of one of these genes compared to a housekeeping gene
#a lot with no expression due to dropout from biological/technical reasons
VlnPlot(decidua, features = c("TMSB4X","PGK1"))

#----Cell cycle----
#cell cycle specific genes
cc.genes.updated.2019
#guess which cell cycle phase each cell is in
CellCycleScoring(decidua, s.features = cc.genes.updated.2019$s.genes, 
                 g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE) -> decidua

decidua[[]]

#----Variable genes----
#compute the mean-variance relationship of each gene and chooses the 
#2000 genes with highest variance. 
decidua <- FindVariableFeatures(decidua, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(decidua), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(decidua)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
x11()
plot1 + plot2

#----Scaling----
#scale data so mean expression across cells is zero and variance is 1
all.genes <- rownames(decidua)
decidua <- ScaleData(decidua, features = all.genes)
#look at scaled values
decidua@assays$RNA
decidua[["RNA"]]$scale.data["CCL21",]

#----Dimension Reduction----
#PCA performed on the variable features only
decidua <- RunPCA(decidua, features = VariableFeatures(object = decidua))
# Examine and visualize PCA results a few different ways
print(decidua[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(decidua, dims = 1:2, reduction = "pca")
DimPlot(decidua, reduction = "pca")
#with ndims we can choose how many PC to plot
ElbowPlot(decidua, ndims=20)+ #somewhere between 10 and 15
  theme(plot.background = element_rect(fill = "#F9FCF3"))
#different method to choose pc
pc.touse <- (decidua$pca@stdev)^2
pc.touse <- pc.touse/sum(pc.touse)
pc.touse <- cumsum(pc.touse)[1:50]
pc.touse <- min(which(pc.touse>=0.75))
pc.touse #12

#----Clustering----
#construct knn graph based on euclidean distance in PCA space
decidua12 <- FindNeighbors(decidua, dims = 1:12) #dims decided in prev step
decidua15 <- FindNeighbors(decidua, dims = 1:15) #dims decided in prev step
#iteratively group cells together into clusters
#larger resolution better for larger datasets
decidua12_5 <- FindClusters(decidua12, resolution = 0.5) #can be adjusted
decidua12_25 <- FindClusters(decidua12, resolution = 0.25) #can be adjusted
decidua15_25 <- FindClusters(decidua15, resolution = 0.25) #can be adjusted

# Look at cluster IDs of the first 5 cells
head(Idents(decidua12_25), 5)
head(Idents(decidua15_25), 5)
head(decidua12_1[[]],5)

#clustered plotted in space of first two PCA components
DimPlot(decidua12_5, reduction = "pca")
DimPlot(decidua12_25, reduction = "pca")
DimPlot(decidua15_25, reduction = "pca")
#plot using t_sne
decidua12_5 <- RunTSNE(decidua12_5, dims=1:12) #same dims used for clustering
DimPlot(decidua12_5, reduction = "tsne")
decidua12_25 <- RunTSNE(decidua12_25, dims=1:12) #same dims used for clustering
DimPlot(decidua12_25, reduction = "tsne")
decidua15_25 <- RunTSNE(decidua15_25, dims=1:15) #same dims used for clustering
DimPlot(decidua15_25, reduction = "tsne")
#plot with umap
decidua12_5 <- RunUMAP(decidua12_5, dims = 1:12)
DimPlot(decidua12_5, reduction = "umap")
decidua12_25 <- RunUMAP(decidua12_25, dims = 1:12)
DimPlot(decidua12_25, reduction = "umap")+
  theme(plot.background = element_rect(fill = "#F9FCF3"))
decidua15_25 <- RunUMAP(decidua15_25, dims = 1:15)
DimPlot(decidua15_25, reduction = "umap")+
  theme(plot.background = element_rect(fill = "#F9FCF3"))

#from here on, I continue with the 12 dimension 0.25 resolution one,
#as I decided it would be better since the clusters are more defined

#check whether quality parameters influenced clustering
#keep in mind that if a cell type is still able to be assigned
#it is not a big problem
vclust1<-VlnPlot(decidua12_25,features="nCount_RNA")+
  ggtitle("Number of Reads")+
  theme(plot.background = element_rect(fill = "#F9FCF3"))+
  scale_fill_manual(values = c(colorpallete, colorpallete[1:3]))
vclust2<-VlnPlot(decidua12_25,features="nFeature_RNA")+
  ggtitle("Number of Genes")+
  theme(plot.background = element_rect(fill = "#F9FCF3"))+
  scale_fill_manual(values = c(colorpallete, colorpallete[1:3]))
vclust3<-VlnPlot(decidua12_25,features="percent.mt")+
  ggtitle("Mitochondrial %")+
  theme(plot.background = element_rect(fill = "#F9FCF3"))+
  scale_fill_manual(values = c(colorpallete, colorpallete[1:3]))
vclust4<-VlnPlot(decidua12_25,features="percent.rbp")+
  ggtitle("Ribosomal Protein %")+
  theme(plot.background = element_rect(fill = "#F9FCF3"))+
  scale_fill_manual(values = c(colorpallete, colorpallete[1:3]))
ggarrange(vclust1, vclust2, vclust3, vclust4,
          ncol = 2, nrow = 2)

vclust1<-VlnPlot(decidua15_25,features="nCount_RNA")+
  ggtitle("Number of Reads")+
  theme(plot.background = element_rect(fill = "#F9FCF3"))+
  scale_fill_manual(values = c(colorpallete, colorpallete[1:5]))
vclust2<-VlnPlot(decidua15_25,features="nFeature_RNA")+
  ggtitle("Number of Genes")+
  theme(plot.background = element_rect(fill = "#F9FCF3"))+
  scale_fill_manual(values = c(colorpallete, colorpallete[1:5]))
vclust3<-VlnPlot(decidua15_25,features="percent.mt")+
  ggtitle("Mitochondrial %")+
  theme(plot.background = element_rect(fill = "#F9FCF3"))+
  scale_fill_manual(values = c(colorpallete, colorpallete[1:5]))
vclust4<-VlnPlot(decidua15_25,features="percent.rbp")+
  ggtitle("Ribosomal Protein %")+
  theme(plot.background = element_rect(fill = "#F9FCF3"))+
  scale_fill_manual(values = c(colorpallete, colorpallete[1:5]))
ggarrange(vclust1, vclust2, vclust3, vclust4,
          ncol = 2, nrow = 2)

#check whether cell cycle influenced clustering
decidua12_25@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")

#----Assigning types to clusters----
#one vs all comparison
#we return only genes "over expressed", found in at least 25% of the cells, and with a logFC threshold of at least 0.25
decidua12_25.markers <- FindAllMarkers(decidua12_25, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#output top 5 genes for each cluster, sorted by logfc
top5all_12_25 <- decidua12_25.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)


#plot markers with heatmap
FeaturePlot(decidua12_25, features = "NKG7",
            cols = c(colorpallete[6], colorpallete[7])) +
  theme(plot.background = element_rect(fill = "#F9FCF3"))
FeaturePlot(decidua12_25, features = "TIMP3",
            cols = c(colorpallete[6], colorpallete[7])) +
  theme(plot.background = element_rect(fill = "#F9FCF3"))
FeaturePlot(decidua12_25, features = c("PRL", "IGFBP1", "NKG7", "CD3D"))

decidua12_25.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10_12_25
DoHeatmap(decidua12_25, features = top10_12_25$gene) +
  scale_fill_gradientn(colors=c(colorpallete[7], colorpallete[3]))+
  NoLegend()+
  theme(plot.background = element_rect(fill = "#F9FCF3"))

#here we see come clusters are very similar
#overlaps:
#----3, 6, 11----
cluster3AND611.markers <- FindMarkers(decidua12_25, ident.1 = c(3,6,11), min.pct = 0.25, test.use = "wilcox")
cluster3AND611.markers <- cluster3AND611.markers[order(-cluster3AND611.markers$avg_log2FC),]
head(cluster3AND611.markers, n = 10)
#NK Cells (KLRC1, GNLY), Gamma Delta T Cells (GNLY)

#3 v 6,11
cluster3611.markers <- FindMarkers(decidua12_25, ident.1 = 3, ident.2 = c(6,11), min.pct = 0.25, test.use = "wilcox")
#overexpressed by cluster 3 but not 6 or 11
cluster3611.markers <- cluster3611.markers[order(-cluster3611.markers$avg_log2FC),]
head(cluster3611.markers, n = 15)
#Monocytes (ZFP36L2 nc)

#6 v 3,11
cluster6311.markers <- FindMarkers(decidua12_25, ident.1 = 6, ident.2 = c(3,11), min.pct = 0.25, test.use = "wilcox")
#overexpressed by cluster 6 but not 3 or 11
cluster6311.markers <- cluster6311.markers[order(-cluster6311.markers$avg_log2FC),]
head(cluster6311.markers, n = 15)
#no relevant markers

#11 v 3,6
cluster1136.markers <- FindMarkers(decidua12_25, ident.1 = 11, ident.2 = c(3,6), min.pct = 0.25, test.use = "wilcox")
#overexpressed by cluster 11 but not 3 or 6
cluster1136.markers <- cluster1136.markers[order(-cluster1136.markers$avg_log2FC),]
head(cluster1136.markers, n = 15)
#Gamma Delta T Cells: marker TOP2A, TROAP

#overexpressed by cluster 0,6 but not 12
cluster1136.markers <- cluster1136.markers[order(cluster1136.markers$avg_log2FC),]
head(cluster1136.markers, n = 10)
#NK Cells: markers KLRB1

#Final Labels:
#3:NK Cells, 6:NK Cells, 11:Gamma Delta T Cells

#----0, 1, 2----
cluster0.markers <- FindMarkers(decidua12_25, ident.1 = 0, min.pct = 0.25, test.use = "wilcox")
cluster0.markers <- cluster0.markers[order(-cluster0.markers$avg_log2FC),]
head(cluster0.markers, n = 10)
#No relevant Markers

cluster1.markers <- FindMarkers(decidua12_25, ident.1 = 1, min.pct = 0.25, test.use = "wilcox")
cluster1.markers <- cluster1.markers[order(-cluster1.markers$avg_log2FC),]
head(cluster1.markers, n = 10)
#Fibroblasts (ADAMTS5)

cluster2.markers <- FindMarkers(decidua12_25, ident.1 = 2, min.pct = 0.25, test.use = "wilcox")
cluster2.markers <- cluster2.markers[order(-cluster2.markers$avg_log2FC),]
head(cluster2.markers, n = 15)
#Fibroblasts (ADAM33, nc)

#0 and 1 and 2
cluster0AND12.markers <- FindMarkers(decidua12_25, ident.1 = c(0,1,2), min.pct = 0.25, test.use = "wilcox")
cluster0AND12.markers <- cluster0AND12.markers[order(-cluster0AND12.markers$avg_log2FC),]
head(cluster0AND12.markers, n = 10)
#Fibroblasts(SCARA5)

#0 v 1,2
cluster012.markers <- FindMarkers(decidua12_25, ident.1 = 0, ident.2 = c(1,2), min.pct = 0.25, test.use = "wilcox")
#overexpressed by cluster 0 but not 1,2
cluster012.markers <- cluster012.markers[order(-cluster012.markers$avg_log2FC),]
head(cluster012.markers, n = 10)
#Ductal Cells (SLPI*)

#1 v 0,2
cluster102.markers <- FindMarkers(decidua12_25, ident.1 = 1, ident.2 = c(0,2), min.pct = 0.25, test.use = "wilcox")
#overexpressed by cluster 1 but not 0,2
cluster102.markers <- cluster102.markers[order(-cluster102.markers$avg_log2FC),]
head(cluster102.markers, n = 10)
#Fibroblasts (ADAMTS5)

#2 v 0,1
cluster201.markers <- FindMarkers(decidua12_25, ident.1 = 2, ident.2 = c(0,1), min.pct = 0.25, test.use = "wilcox")
#overexpressed by cluster 2 but not 0,1
cluster201.markers <- cluster201.markers[order(-cluster201.markers$avg_log2FC),]
head(cluster201.markers, n = 10)
#Erythroid-like and erythroid precursor cells (THBS1)
#Fibroblasts (ADAM33 nc)

#overexpressed by cluster 0,1 but not 2
cluster201.markers <- cluster201.markers[order(cluster201.markers$avg_log2FC),]
head(cluster201.markers, n = 10)

#Final Labels
#0:Decidual Stromal Cells, 1:Fibroblasts , 2:Fibroblasts


#----4----
cluster4.markers <- FindMarkers(decidua12_25, ident.1 = 4, min.pct = 0.25, test.use = "wilcox")
cluster4.markers <- cluster4.markers[order(-cluster4.markers$avg_log2FC),]
head(cluster4.markers, n = 15)
#Dendritic Cells (TREM2)

#Final Label: Dendritic Cells

#----5----
cluster5.markers <- FindMarkers(decidua12_25, ident.1 = 5, min.pct = 0.25, test.use = "wilcox")
cluster5.markers <- cluster5.markers[order(-cluster5.markers$avg_log2FC),]
head(cluster5.markers, n = 15)
#T Cells (IL7R)

#Final Label: T Cells

#----7----
cluster7.markers <- FindMarkers(decidua12_25, ident.1 = 7, min.pct = 0.25, test.use = "wilcox")
cluster7.markers <- cluster7.markers[order(-cluster7.markers$avg_log2FC),]
head(cluster7.markers, n = 15)
#Cholangiocytes (ELF3, LCN2), Ductal Cells (PERP, WFDC2)

#Final Label: Ductal Cells

#----8, 10----
cluster8.markers <- FindMarkers(decidua12_25, ident.1 = 8, min.pct = 0.25, test.use = "wilcox")
cluster8.markers <- cluster8.markers[order(-cluster8.markers$avg_log2FC),]
head(cluster8.markers, n = 10)
#Endothelial Cells (MMRN2)

cluster10.markers <- FindMarkers(decidua12_25, ident.1 = 10, min.pct = 0.25, test.use = "wilcox")
cluster10.markers <- cluster10.markers[order(-cluster10.markers$avg_log2FC),]
head(cluster10.markers, n = 10)
#Endothelial Cells (STAB2)

#8 and 10
cluster8AND10.markers <- FindMarkers(decidua12_25, ident.1 = c(8,10), min.pct = 0.25, test.use = "wilcox")
cluster8AND10.markers <- cluster8AND10.markers[order(-cluster8AND10.markers$avg_log2FC),]
head(cluster8AND10.markers, n = 10)
#Endothelial Cells (ERG)

#8 v 10
cluster810.markers <- FindMarkers(decidua12_25, ident.1 = 8, ident.2 = 10, min.pct = 0.25, test.use = "wilcox")
#overexpressed by cluster 8 but not 10
cluster810.markers <- cluster810.markers[order(-cluster810.markers$avg_log2FC),]
head(cluster810.markers, n = 10)
#Endothelial Cells (ACKR1)

#overexpressed by cluster 10 but not 8
cluster810.markers <- cluster810.markers[order(cluster810.markers$avg_log2FC),]
head(cluster810.markers, n = 10)
#Endothelial Cells (TMEM100)

#Final Labels:
#8:Endothelial Cells, 10: Endothelial Cells

#----9----
cluster9.markers <- FindMarkers(decidua12_25, ident.1 = 9, min.pct = 0.25, test.use = "wilcox")
cluster9.markers <- cluster9.markers[order(-cluster9.markers$avg_log2FC),]
head(cluster9.markers, n = 15)
#Endothelial Cells/Pancreatic Stellate Cells (RGS5), Fibroblasts (FRZB)
#Pancreatic Stellate Cells (NDUFA4L2)

#Final Labels
#9:Pancreatic Stellate Cells





#----Cell Types and Markers----
#NK Cells: cluster 3,6 (Marker: KLRB1)
#Gamma Delta T Cells: cluster 11 (Marker: TOP2A)
#Fibroblasts: cluster 1,2 (Marker: SCARA5)
#Decidual Stomal Cells: cluster 0 (Marker: PRL)
#Endothelial Cells: cluster 8,10 (Marker: ERG)
#"Pancreatic Stellate cells: cluster 9 (Marker: NDUFA4L2, RGS5)
#Dendritic Cells: cluster 4 (Marker: TREM2)
#Ductal Cells: Cluster 7 (Marker: ELF3)
#T Cells: Cluster 5 (Marker: IL7R)

marker_list <- c("NKG7", "TOP2A", "SCARA5", "PRL", "ERG", "NDUFA4L2",
                 "TREM2", "ELF3", "CD3D")

#plotting of marker genes
FeaturePlot(decidua12_25, features = marker_list, cols = c(colorpallete[6], colorpallete[7])) &
  theme(plot.background = element_rect(fill = "#F9FCF3"),
        panel.background = element_rect(fill = "#F9FCF3"))
VlnPlot(decidua12_25, features = marker_list)
DotPlot(decidua12_25, features = marker_list, cols=c("lightgrey", colorpallete[2]))+
  theme(plot.background = element_rect(fill = "#F9FCF3"))


#add new labels to seurat obj
new.cluster.ids <- c("Decidualized Stromal", "Fibroblast", "Fibroblast",
                     "NK", "Dendritic", "T",
                     "NK", "Ductal", "Endothelial",
                     "Pancreatic Stellate", "Endothelial", "Gamma Delta T")
names(new.cluster.ids) <- levels(decidua12_25)
decidua12_25 <- RenameIdents(decidua12_25, new.cluster.ids)
DimPlot(decidua12_25, reduction = "umap", label = TRUE, pt.size = 0.5) +
  theme(plot.background = element_rect(fill = "#F9FCF3"))+
  NoLegend()

DimPlot(decidua12_25, reduction = "tsne", label = TRUE, pt.size = 0.5) +
  theme(plot.background = element_rect(fill = "#F9FCF3"))+
  NoLegend()
