{
library(Seurat)
library(ggplot2)
library(plotly)
library(dplyr)
library(Matrix)
}

#Set working directory
setwd('')

#Expression matrix file name
path_data=''
raw_counts=read.table(path_data,header=T,sep=',',row.names=1)

#Once the Seurat object is created, the data processing can be stopped at any
#time and the Seurat object exported in a rds file. This file can then be loaded
#to resume the processing at the step it was left to

#Export the seurat object 
saveRDS(mydata, '')

#Import the seurat object
mydata=readRDS('')


########################### Data processing ######################################

#Load data in a Seurat object and free RAM by deleting the expression matrix
mydata=CreateSeuratObject(raw_counts,min.cells=3,min.genes=200,project = 'Seurat pipeline')
ngenes=nrow(raw_counts)
rm(raw_counts)

#Data normalization
mydata <- NormalizeData(mydata, normalization.method = 'LogNormalize', scale.factor = 10000)
mydata<- FindVariableFeatures(mydata, selection.method = 'vst', nfeatures =ngenes)
all.genes <- rownames(mydata)
mydata <- ScaleData(mydata, features = all.genes,do.scale=TRUE,do.center=TRUE)

#PCA
mydata <- RunPCA(mydata, npcs=100)

#Some visualization
VizDimLoadings(mydata, dims = 1:2, reduction = 'pca')
DimPlot(mydata, reduction = 'pca')
DimHeatmap(mydata, dims = 1, cells = 500, balanced = TRUE)

#UMAP/Louvain clustering
dims=100
mydata=FindNeighbors(mydata,reduction='pca',dims=1:dims)
mydata=FindClusters(mydata,resolution=1)
mydata=RunUMAP(mydata, dims = 1:10,n.neighbors=20,metric='cosine',min.dist=0.1)
DimPlot(mydata, reduction = 'umap')



###Find markers genes in the clusters created, parameters are adjustable
#only.pos(TRUE or FALSE): Only genes overexpressed in the cluster or not  
#min.pct: % of barcodes in the cluster in which the gene must be expressed
#logfc.threshold: Minimum threshold of the gene differential expression with other clusters
cluster.markers=FindAllMarkers(mydata,only.pos=T, min.pct = 0.25, logfc.threshold = 1)
cluster.markers=cluster.markers[cluster.markers$p_val_adj<0.05,]

#Write table of markers genes per cluster
out_markers=''
write.csv(cluster.markers,out_markers)



#######################Analysis of specific genes expression#########################
genes=c() #1 or more genes can be entered

#Gene expression in the UMAP plot
FeaturePlot(mydata,features=genes

#Heatmap of genes expression
DoHeatmap(subset(mydata, downsample = 1000), features = features, size = 2)