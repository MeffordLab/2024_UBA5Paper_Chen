# Analyzing HChen's scRNAseq data (10X) from human cortical spheroids
# November 11, 2024
# Christy LaFlamme
##########################################################################################
# LOAD LIBRARIES
##########################################################################################
library(Seurat) # version 5.1.0
library(dplyr)
library(Matrix)
library(ggplot2)
library(data.table)
library(patchwork)
library(SeuratObject)
library(SeuratData)
library(multtest)
library(metap)
library(stringr)
##########################################################################################
# FUNCTIONS
##########################################################################################
# function to remove lowly expressed gene from raw Seurat object
removeLowlyExpGenes <- function(raw, number = 1) {
  counts <- GetAssayData(raw, assay = "RNA") # get the count data
  total <- as.data.frame(rowSums(counts)) # sum the rows and convert to dataframe
  cutoff <- as.numeric(as.numeric(number)+1) # get the cutoff number (user input for max + 1)
  number <- length(total[which(total$`rowSums(counts)` < cutoff),]) # get the number of genes to remove
  print(paste("Removing", number, "genes that have only 0 or 1 count(s) across all the cells.", sep = " ")) # print statement
  
  # subset the features in the raw Seurat object
  indexes <- which(total$`rowSums(counts)` < cutoff)
  genes <- row.names(total)[indexes]
  counts <- counts[-(which(rownames(counts) %in% genes)),]
  raw <- subset(raw, features = rownames(counts))
  
  # return the raw Seurat object
  return(raw)
}
##########################################################################################
# function to remove genes starting with ^ENSG and ^LINC from the raw count matrix Seurat object
removeENSGandLINCGenes <- function(raw) {
  
  # subset the features of the raw object
  features = NULL
  features <- Features(raw)
  print(paste0("The raw object has a total of ", length(features), " features before filtering out ^ENSG and ^LINC genes."))
  
  # determine the number of features starting with ENSG and LINC
  num.ensg <- table(str_extract(features, "^ENSG"))
  num.linc <- table(str_extract(features, "^LINC"))
  
  # remove genes starting with ^ENSG
  not.ensg <- which(is.na(str_extract(features, "^ENSG"))) # these are the indexes of those that are na's (in other words those that aren't ENSG)
  genes.keep <- features[not.ensg] # get the gene names to keep
  num.genes.remain <- length(genes.keep) # store the number of genes remaining
  print(paste0(num.ensg, " genes removed starting with 'ENSG'. There are ", num.genes.remain, " genes remaining."))
  
  raw.ensg.removed <- subset(raw, features = genes.keep) # subset out the genes
  
  # remove genes starting with ^LINC
  features <- Features(raw.ensg.removed)
  not.linc <- which(is.na(str_extract(features, "^LINC"))) # these are the indexes of those that are na's (in other words those that aren't LINC)
  genes.keep <- features[not.linc] # get the gene names to keep
  num.genes.remain <- length(genes.keep) # store the number of genes remaining
  print(paste0(num.linc, " genes removed starting with 'LINC'. There are ", num.genes.remain, " genes remaining."))
  
  raw.ensg.linc.removed <- subset(raw.ensg.removed, features = genes.keep) # subset out the genes
  
  # return the raw Seurat object
  return(raw.ensg.linc.removed)
}
##################################################
# create Seurat Object function to load 10X data into a raw seurat object and perform QC

createSO <- function(data, project) {
  
  raw <- CreateSeuratObject(counts = data, project = project, min.cells = 3, min.features = 200)
  
  # MITOCHONDRIA
  # get info for mitochondrial DNA percentage
  raw[["percent.mt"]] <- PercentageFeatureSet(raw, pattern = "^MT-")
  
  # Visualize mitochondrial QC metrics as a violin plot
  pdf("./raw_vlnplot.pdf")
  vln <- VlnPlot(raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(vln)
  dev.off()
  
  # FEATURE SCATTER
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  pdf("./raw_featurescatter.pdf")
  plot1 <- FeatureScatter(raw, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2)
  dev.off()
  
  # FILTERING
  # from Helen's manuscript (Chen et al., 2024 BioRxiv)
  # 1) cells with fewer than 1,000 or more than 20,000 unique molecular
  # identifiers were removed
  # 2) cells expressing less than 500 or larger than 5,000 unique genes were
  # considered outliers and discarded
  # 3) cells with a mitochondrial transcript proportion higher than 20%
  # were filtered out
  raw <- subset(raw, subset = nCount_RNA > 1000 & nCount_RNA < 20000 & nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20)
  
  # REMOVE LOWLY EXPRESSING GENES
  # remove genes that only have 1 count or less across all cells (lowly expressed genes)
  # input: raw Seurat object, the maximum number of genes to filter cells (i.e. 1 means that cells with 0 or 1 count will be removed)
  raw <- removeLowlyExpGenes(raw, number = 1)
  
  #REMOVE ENSG and LINC genes (RNA genes, pseudogenes)
  # Added 09-12-24
  raw <- removeENSGandLINCGenes(raw)
  
  # return the raw seurat object
  return(raw)
}
##################################################
# normSO to normalize Seurat objects

normSO <- function(raw){
  
  # normalize the data (default parameters)
  norm <- NormalizeData(raw, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # feature selection
  norm <- FindVariableFeatures(norm, selection.method = "vst", nfeatures = 2000)
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(norm), 10)
  
  # plot variable features with and without labels
  pdf("./norm_variablefeatures_plot.pdf", width = 10, height = 5)
  plot1 <- VariableFeaturePlot(norm)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  print(plot1 + plot2)
  dev.off()
  
  # scaling the data
  all.genes <- rownames(norm)
  norm <- ScaleData(norm, features = all.genes)
  
  # run PCA
  norm <- RunPCA(norm, features = VariableFeatures(object = norm))
  
  # examine PCA
  print(norm[["pca"]], dims = 1:5, nfeatures = 5)
  
  # vizualize the dimensions
  
  pdf("./norm_vizdimloadings.pdf")
  vizdim <- VizDimLoadings(norm, dims = 1:2, reduction = "pca")
  print(vizdim)
  dev.off()
  
  pdf("./norm_dimplot.pdf")
  dimplot <- DimPlot(norm, reduction = "pca") + NoLegend()
  print(dimplot)
  dev.off()
  
  pdf("./norm_dimheatmap_dim1.pdf")
  dimheat <- DimHeatmap(norm, dims = 1, cells = 500, balanced = TRUE)
  print(dimheat)
  dev.off()
  
  pdf("./norm_dimheatmap_dim1-30")
  dimheat30 <- DimHeatmap(norm, dims = 1:30, cells = 500, balanced = TRUE)
  print(dimheat30)
  dev.off()
  
  pdf("./norm_elbowplot.pdf")
  elbow <- ElbowPlot(norm, ndims = 40)
  print(elbow)
  dev.off()
  
  # return the norm object
  return(norm)
}

# function for quick within sample umap clustering
clusterSO <- function(norm, dims = 30){
  
  # run single sample clustering for QC purposes
  norm <- FindNeighbors(norm, dims = 1:as.numeric(dims))
  norm <- FindClusters(norm, resolution = 0.5)
  head(Idents(norm), 5)
  norm <- RunUMAP(norm, dims = 1:as.numeric(dims))
  
  # save the object
  saveRDS(norm, file = "./norm_umap.rds")
  
  pdf("./umap_dimplot.pdf")
  umap.pl <- DimPlot(norm, reduction = "umap")
  print(umap.pl)
  dev.off()
  
  # finding cluster biomarkers
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  markers <- FindAllMarkers(norm, only.pos = TRUE)
  write.table(markers, file = "./ManuscriptFilters_cluster_markers.txt", sep = "\t", quote = F, row.names = F)
  
  # get top 40 marker genes
  top40 <- topXgenes(markers, log2FC.column = 2, cluster.column = 6, LFCcutoff = 0.1, padjcutoff = 0.05, numbergenes = 40) 
  write.table(top40, file = "./ManuscriptFilters_cluster_markers_top40.txt", sep = "\t", quote = F, row.names = F)
  
  return(norm)
}
##################################################
# function for selecting the top number of genes from marker gene list
# select the number of genes to be selected for each cluster on which to perform cell type identity
# for clusters starting with zero!
# topXgenes <- function(markers, log2FC.column = 2, cluster.column = 6, LFCcutoff = 1, padjcutoff = 0.05, numbergenes = 40) {
#   
#   # first filter and rank the marker genes across the clusters
#   markers.filtered <- markers %>% filter(avg_log2FC > as.numeric(LFCcutoff)) %>% filter(p_val_adj < as.numeric(padjcutoff))
#   markers.filtered <- markers.filtered[order(markers.filtered[,as.numeric(cluster.column)], -markers.filtered[,as.numeric(log2FC.column)]),]
#   
#   topX <- NULL
#   topX = as.data.frame(markers.filtered[markers.filtered$cluster == 0, "gene"][1:as.numeric(numbergenes)])
#   colnames(topX) <- "cluster0" #assumption that first cluster is cluster0
#   
#   # loop over the clusters and grab the first X number of genes to create a nice and neat dataframe
#   for (i in 1:as.numeric(length(table(markers.filtered$cluster))-1)) {
#     topX[,paste0("cluster",i)] <- markers.filtered[markers.filtered$cluster == i, "gene"][1:as.numeric(numbergenes)]
#   }
#   
#   return(topX)
#   
# }
##################################################
# for clusters starting with one!!!
topXgenes <- function(markers, log2FC.column = 2, cluster.column = 6, LFCcutoff = 1, padjcutoff = 0.05, numbergenes = 40) {
  
  # first filter and rank the marker genes across the clusters
  markers.filtered <- markers %>% filter(avg_log2FC > as.numeric(LFCcutoff)) %>% filter(p_val_adj < as.numeric(padjcutoff))
  markers.filtered <- markers.filtered[order(markers.filtered[,as.numeric(cluster.column)], -markers.filtered[,as.numeric(log2FC.column)]),]
  
  topX <- NULL
  topX = as.data.frame(markers.filtered[markers.filtered$cluster == 1, "gene"][1:as.numeric(numbergenes)])
  colnames(topX) <- "cluster1" #assumption that first cluster is cluster0
  
  # loop over the clusters and grab the first X number of genes to create a nice and neat dataframe
  for (i in 1:as.numeric(length(table(markers.filtered$cluster)))) {
    topX[,paste0("cluster",i)] <- markers.filtered[markers.filtered$cluster == i, "gene"][1:as.numeric(numbergenes)]
  }
  
  return(topX)
  
}
##########################################################################################################################################################
# Load the data from 10X
##################################################

# set working directory for data upload - mom (e8278-2) from family #1 (SRM301059)
setwd("/research_jude/rgs01_jude/groups/meffogrp/projects/meffogrp_auto/common/cab/10XSCRNASEQ/MEFFO-301059-10XSCRNASEQ/")
data <- Read10X(data.dir = "./2611941_HC020/196759465/filtered_feature_bc_matrix/")
data1 <- Read10X(data.dir = "./2611942_HC021/196759467/filtered_feature_bc_matrix/")

# set working directory for data upload - proband (e4559-1) from family #2 (836687
setwd("/research_jude/rgs01_jude/groups/meffogrp/projects/meffogrp_auto/common/cab/10XSCRNASEQ/MEFFO-836687-10XSCRNASEQ/")
data2 <- Read10X(data.dir = "./3209730_HC066/239872970/filtered_feature_bc_matrix/")
data3 <- Read10X(data.dir = "./3209731_HC067/239872972/filtered_feature_bc_matrix/")
data4 <- Read10X(data.dir = "./3209734_HC070/239872978/filtered_feature_bc_matrix/")

##################################################
# MOTHER FROM FAMILY #1

# e8278-2 CO D100 cl1 mom rep #1
setwd("/research_jude/rgs01_jude/groups/meffogrp/projects/EpilepsyOmics/common/10X/MEFFO-301059-10XSCRNASEQ_and_MEFFO-836687-10XSCRNASEQ/2611941_HC020/")
raw <- createSO(data, "HC020")
# saveRDS(raw, "./2611941_HC020_raw.RDS")
# raw <- readRDS("./2611941_HC020_raw.RDS")
norm <- normSO(raw)
# saveRDS(norm, "./2611941_HC020_norm.RDS")
norm <- readRDS("./2611941_HC020_norm.RDS")
norm.umap <- clusterSO(norm, dims = 30)
# saveRDS(norm.umap, "./2611941_HC020_norm.umap.RDS")

# e8278-2 CO D100 cl1 mom rep #2
setwd("/research_jude/rgs01_jude/groups/meffogrp/projects/EpilepsyOmics/common/10X/MEFFO-301059-10XSCRNASEQ_and_MEFFO-836687-10XSCRNASEQ/2611942_HC021/")
raw1 <- createSO(data1, "HC021")
# saveRDS(raw1, "./2611942_HC021_raw.RDS")
# raw <- readRDS("./2611942_HC021_raw.RDS")
norm1 <- normSO(raw1)
# saveRDS(norm1, "./2611942_HC021_norm.RDS")
norm1 <- readRDS("./2611942_HC021_norm.RDS")
norm.umap1 <- clusterSO(norm1, dims = 30)
# saveRDS(norm.umap1, "./2611942_HC021_norm.umap.RDS")

##################################################
# PROBAND FROM FAMILY #2

# e4559-1 CO D100 proband rep#1
setwd("/research_jude/rgs01_jude/groups/meffogrp/projects/EpilepsyOmics/common/10X/MEFFO-301059-10XSCRNASEQ_and_MEFFO-836687-10XSCRNASEQ/3209730_HC066/")
raw2 <- createSO(data2, "HC066")
# saveRDS(raw2, "./3209730_HC066_raw.RDS")
# raw2 <- readRDS("./3209730_HC066_raw.RDS")
norm2 <- normSO(raw2)
# saveRDS(norm2, "./3209730_HC066_norm.RDS")
norm2 <- readRDS("./3209730_HC066_norm.RDS")
norm.umap2 <- clusterSO(norm2, dims = 30)
# saveRDS(norm.umap2, "./3209730_HC066_norm.umap.RDS")

# e4559-1 CO D100 proband rep#2
setwd("/research_jude/rgs01_jude/groups/meffogrp/projects/EpilepsyOmics/common/10X/MEFFO-301059-10XSCRNASEQ_and_MEFFO-836687-10XSCRNASEQ/3209731_HC067/")
raw3 <- createSO(data3, "HC067")
# saveRDS(raw3, "./3209731_HC067_raw.RDS")
raw3 <- readRDS("./3209731_HC067_raw.RDS")
norm3 <- normSO(raw3)
# saveRDS(norm3, "./3209731_HC067_norm.RDS")
norm3 <- readRDS("./3209731_HC067_norm.RDS")
norm.umap3 <- clusterSO(norm3, dims = 30)
# saveRDS(norm.umap3, "./3209731_HC067_norm.umap.RDS")

# e4559-1 CO D100 proband rep#3
setwd("/research_jude/rgs01_jude/groups/meffogrp/projects/EpilepsyOmics/common/10X/MEFFO-301059-10XSCRNASEQ_and_MEFFO-836687-10XSCRNASEQ/3209734_HC070/")
raw4 <- createSO(data4, "HC070")
# saveRDS(raw4, "./3209734_HC070_raw.RDS")
raw4 <- readRDS("./3209734_HC070_raw.RDS")
norm4 <- normSO(raw4)
# saveRDS(norm4, "./3209734_HC070_norm.RDS")
# norm4 <- readRDS("./3209734_HC070_norm.RDS")
norm.umap4 <- clusterSO(norm4, dims = 30)
# saveRDS(norm.umap4, "./3209734_HC070_norm.umap.RDS")

setwd("/research_jude/rgs01_jude/groups/meffogrp/projects/EpilepsyOmics/common/10X/MEFFO-301059-10XSCRNASEQ_and_MEFFO-836687-10XSCRNASEQ/HC067_and_HC070/")
# merge HC067 and HC070 before merging
raw34 <- merge(raw3, y = raw4, add.cell.ids = c("HC067", "HC070"))
raw34$orig.ident <- "HC067_and_HC070"
norm34 <- normSO(raw34)
saveRDS(norm34, "./HC067_and_HC070_norm.RDS")
# norm34 <- readRDS("./HC067_and_HC070_norm.RDS")
norm.umap34 <- clusterSO(norm34, dims = 30)
saveRDS(norm.umap34, "./HC067_and_HC070_norm.umap.RDS")

############################################
# filter based on FOXG1 > 0
############################################
setwd("/research_jude/rgs01_jude/groups/meffogrp/projects/EpilepsyOmics/common/10X/MEFFO-301059-10XSCRNASEQ_and_MEFFO-836687-10XSCRNASEQ/HC067_HC070_merged_FOXG1_filtered_output/")
norm.foxg1 <- subset(x = norm, subset = FOXG1 > 0)
# saveRDS(norm.foxg1, "./HC020_raw_FOXG1_filtered.RDS")
norm.foxg1 <- readRDS("./HC020_raw_FOXG1_filtered.RDS")

norm1.foxg1 <- subset(x = norm1, subset = FOXG1 > 0)
# saveRDS(norm1.foxg1, "./HC021_raw_FOXG1_filtered.RDS")
norm1.foxg1 <- readRDS("./HC021_raw_FOXG1_filtered.RDS")

norm2.foxg1 <- subset(x = norm2, subset = FOXG1 > 0)
# saveRDS(norm2.foxg1, "./HC066_raw_FOXG1_filtered.RDS")
norm2.foxg1 <- readRDS("./HC066_raw_FOXG1_filtered.RDS")

norm34.foxg1 <- subset(x = norm34, subset = FOXG1 > 0)
# saveRDS(norm34.foxg1, "./HC067_and_HC070_raw_FOXG1_filtered.RDS")
norm34.foxg1 <- readRDS("./HC067_and_HC070_raw_FOXG1_filtered.RDS")

############################################
# integration and analysis
############################################
setwd("/research_jude/rgs01_jude/groups/meffogrp/projects/EpilepsyOmics/common/10X/MEFFO-301059-10XSCRNASEQ_and_MEFFO-836687-10XSCRNASEQ/HC067_HC070_merged_FOXG1_filtered_output_redo/")
# integrate the data across all the samples so that cells can be clustered and assessed at the same time; also for DGE between conditions across same cell types
######################
# trying by creating a list of separate seurat objects
# followed this tutorial: https://satijalab.org/seurat/archive/v4.3/integration_introduction
norm.list <- list(HC020 = norm.foxg1, HC021 = norm1.foxg1, HC066 = norm2.foxg1, HC067_and_HC070 = norm34.foxg1)

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = norm.list)
brain.anchors <- FindIntegrationAnchors(object.list = norm.list, anchor.features = features)

# this command creates an 'integrated' data assay
norm.combined <- IntegrateData(anchorset = brain.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(norm.combined) <- "integrated"

norm.combined <- ScaleData(norm.combined, verbose = FALSE)
norm.combined <- RunPCA(norm.combined, npcs = 30, verbose = FALSE)
norm.combined <- RunUMAP(norm.combined, reduction = "pca", dims = 1:30)
norm.combined <- FindNeighbors(norm.combined, reduction = "pca", dims = 1:30)
norm.combined <- FindClusters(norm.combined, resolution = 0.5)

# start the cluster number from 1 instead of 0
norm.combined$seurat_clusters <- as.factor(as.numeric(as.character(norm.combined$seurat_clusters)) + 1)

# superimpose 067 and 070 plots; did this before I went back and fully merged them before integration
# table(norm.combined$orig.ident)
# norm.combined$orig.ident <- gsub("HC067","HC1_6770_combined",norm.combined$orig.ident)
# norm.combined$orig.ident <- gsub("HC070","HC1_6770_combined",norm.combined$orig.ident)

# Try specifying which metadata column to use explicitly so that clusters start from 1
Idents(norm.combined) <- "seurat_clusters"

# Visualization
pdf("./norm_combined_UMAP_HC020_HC021_HC066_HC067_HC070_6770_combined_clusters_renumb_FOXG1.pdf", width = 10, height = 5)
p1 <- DimPlot(norm.combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(norm.combined, reduction = "umap", label = TRUE, repel = TRUE)
print(p1 + p2)
dev.off()

pdf("./norm_combined_UMAP_SPLIT_HC020_HC021_HC066_HC067_HC070_6770_combined_clusters_renumb_FOXG1_nolabels.pdf", width = 15, height = 5)
p3 <- DimPlot(norm.combined, reduction = "umap", split.by = "orig.ident", group.by="seurat_clusters") # label = TRUE
print(p3)
dev.off()

# saving data
# saveRDS(norm.combined, file = "./norm_combined_HC020_HC021_HC066_HC067_HC070.RDS")
norm.combined <- readRDS("./norm_combined_HC020_HC021_HC066_HC067_HC070.RDS")
######################
# For performing differential expression after integration, we switch back to the original
DefaultAssay(norm.combined) <- "RNA"

# Once integrative analysis is complete, you can rejoin the layers - which collapses the individual datasets together and recreates the original counts and data layers. 
# You will need to do this before performing any differential expression analysis. However, you can always resplit the layers in case you would like to reperform integrative analysis.
norm.combined <- JoinLayers(norm.combined)

# find markers for every cluster compared to all remaining cells, report only the positive ones
for (i in 1:13) {
  markers <- NULL
  markers <- FindConservedMarkers(norm.combined, ident.1 = as.numeric(i), grouping.var = "orig.ident")
  head(markers)
  write.table(markers, paste0("./norm_combined_FindConservedMarkers_cluster",toString(i),".txt"), quote = F, sep = "\t", row.names = T)
}

# find all markers
markers <- FindAllMarkers(norm.combined, only.pos = TRUE)

write.table(markers, paste0("./norm_combined_FindAllMarkers_ALL.txt"), quote = F, sep = "\t", row.names = T)
# markers <- as.data.frame(read.table("./norm_combined_FindAllMarkers_ALL.txt"))
# markers <- as.data.frame(fread("./norm_combined_FindAllMarkers_ALL.txt", drop = 1))

# filter by p adj < 0.05 and LFC > 1; and then rank by control LFC to determien top 40 marker genes of cell clusters
markers.filtered <- markers %>% filter(avg_log2FC > 1) %>% filter(p_val_adj < 0.05)
markers.filtered <- markers.filtered[order(markers.filtered[,6], -markers.filtered[,2]),]
write.table(markers.filtered, paste0("./norm_combined_FindAllMarkers_filtered_LFC1_padj0.05.txt"), quote = F, sep = "\t", row.names = T)

######################
# use function to get the top genes

# testing out this function and it works! -Christy LaFlamme on March 18, 2024
top40 <- topXgenes(markers, LFCcutoff = 1, padjcutoff = 0.05, numbergenes = 40)
write.table(top40, paste0("./norm_combined_top40genes_per13clusters.txt"), quote = F, sep = "\t", row.names = T)

top100 <- topXgenes(markers, LFCcutoff = 1, padjcutoff = 0.05, numbergenes = 100)
write.table(top100, paste0("./norm_combined_top100genes_per13clusters.txt"), quote = F, sep = "\t", row.names = T)

######################
# now, plot expression of some genes across the umap

# pdf("./norm.combined_FeaturePlot_FOXG1_NKX2-2_PAX6_SHH_GRIN2B.pdf", width = 25, height = 20)
# genexp <- FeaturePlot(norm.combined, features = c("FOXG1", "NKX2-2", "PAX6", "SHH", "GRIN2B"), min.cutoff = "q9", split.by = "orig.ident") & theme(legend.position = "right")
# print(genexp)
# dev.off()
# 
# pdf("./norm.combined_FeaturePlot_Excitatory_Neurons.pdf", width = 25, height = 20)
# genexp <- FeaturePlot(norm.combined, features = c("NRN1", "MEF2C", "GRIN2B", "SLC17A7", "SLC17A6"), min.cutoff = "q9", split.by = "orig.ident") & theme(legend.position = "right")
# print(genexp)
# dev.off()
# 
# pdf("./norm.combined_FeaturePlot_GABA_Interneurons.pdf", width = 25, height = 20)
# genexp <- FeaturePlot(norm.combined, features = c("DLX1", "DLX2", "DLX5", "SST", "NPY", "SP8"), min.cutoff = "q9", split.by = "orig.ident") & theme(legend.position = "right")
# print(genexp)
# dev.off()
# 
# pdf("./norm.combined_FeaturePlot_MainFig2.pdf", width = 25, height = 20)
# genexp <- FeaturePlot(norm.combined, features = c("GAD1", "GAD2", "CALB1", "CALB2", "SCGN", "TBR1", "SATB2"), min.cutoff = "q9", split.by = "orig.ident") & theme(legend.position = "right")
# print(genexp)
# dev.off()

# plots for helen's paper after super-imposing both replicates
norm.combined$orig.ident <- gsub("HC020","e8278-2_combined",norm.combined$orig.ident)
norm.combined$orig.ident <- gsub("HC021","e8278-2_combined",norm.combined$orig.ident)
norm.combined$orig.ident <- gsub("HC066","e4559-1_combined",norm.combined$orig.ident)
norm.combined$orig.ident <- gsub("HC067_and_HC070","e4559-1_combined",norm.combined$orig.ident)

pdf("./norm.combined_FeaturePlot_GAD1_GAD2_CR_SCGN_TBR1_SATB2_combinedsamples.pdf", width = 10, height = 24)
genexp <- FeaturePlot(norm.combined, features = c("GAD1", "GAD2", "CALB2", "SCGN", "TBR1", "SATB2"), min.cutoff = "q9", split.by = "orig.ident") & theme(legend.position = "right")
print(genexp)
dev.off()

pdf("./norm.combined_FeaturePlot_DLX1_DLX2_DLX5_SST_combinedsamples.pdf", width = 10, height = 16)
genexp <- FeaturePlot(norm.combined, features = c("DLX1", "DLX2", "DLX5", "SST"), min.cutoff = "q9", split.by = "orig.ident") & theme(legend.position = "right")
print(genexp)
dev.off()

pdf("./norm.combined_FeaturePlot_EOMES_NEUROG1_combinedsamples.pdf", width = 10, height = 8)
genexp <- FeaturePlot(norm.combined, features = c("EOMES", "NEUROG1"), min.cutoff = "q9", split.by = "orig.ident") & theme(legend.position = "right")
print(genexp)
dev.off()

pdf("./norm.combined_FeaturePlot_NPY_SP8_SP9_VIP_CB_combinedsamples.pdf", width = 10, height = 20)
genexp <- FeaturePlot(norm.combined, features = c("OXT", "SP8", "SP9", "AVP", "CALB1"), min.cutoff = "q9", split.by = "orig.ident") & theme(legend.position = "right")
print(genexp)
dev.off()

pdf("./norm.combined_FeaturePlot_HOPX_combinedsamples.pdf", width = 10, height = 4)
genexp <- FeaturePlot(norm.combined, features = c("HOPX"), min.cutoff = "q9", split.by = "orig.ident") & theme(legend.position = "right")
print(genexp)
dev.off()

pdf("./norm.combined_FeaturePlot_Excitatory_Markers_combinedsamples.pdf", width = 10, height = 20)
genexp <- FeaturePlot(norm.combined, features = c("NRN1", "MEF2C", "GRIN2B", "SLC17A7", "SLC17A6"), min.cutoff = "q9", split.by = "orig.ident") & theme(legend.position = "right")
print(genexp)
dev.off()

pdf("./norm.combined_FeaturePlot_UBA5_combinedsamples.pdf", width = 10, height = 4)
genexp <- FeaturePlot(norm.combined, features = c("UBA5"), min.cutoff = "q9", split.by = "orig.ident") & theme(legend.position = "right")
print(genexp)
dev.off()

pdf("./norm.combined_FeaturePlot_NPY_CALB2_combinedsamples.pdf", width = 10, height = 8)
genexp <- FeaturePlot(norm.combined, features = c("NPY",  "CALB2"), min.cutoff = "q9", split.by = "orig.ident") & theme(legend.position = "right")
print(genexp)
dev.off()

######################
# get average expression value across all clusters for all samples
avgexp <- as.data.frame(AverageExpression(object = norm.combined, group.by = c('orig.ident', 'seurat_clusters'))$RNA)
write.table(avgexp, file = "./norm.combined_avgexp_values_across_clusters_bysample.txt", sep = "\t", quote = F, row.names = T)

uba5_avgexp <- t(as.data.frame(avgexp %>% filter(rownames(avgexp) == "UBA5")))
write.table(uba5_avgexp, file = "./norm.combined_avgexp_values_across_clusters_bysample_UBA5_values.txt", sep = "\t", quote = F, row.names = T)

######################
# now, want to get proportion of cells for each cell type within each orig ident (to compare proportions across HC samples)

cellcounts.bycluster <- as.data.frame(table(paste0(norm.combined@meta.data$orig.ident, norm.combined@meta.data$seurat_clusters)))
write.table(cellcounts.bycluster, file = "./norm.combined_cellcounts.bycluster.txt", sep = "\t", quote = F, row.names=F)

######################
# foxg1
VlnPlot(norm.combined, features = c("FOXG1"), split.by = "orig.ident")
uba5 <- VlnPlot(norm.combined, features = c("UBA5"), split.by = "orig.ident")

# foxg1 data
foxg1 <- FetchData(norm.combined, vars = "FOXG1") 
hist(foxg1$FOXG1)

# foxg1 data
uba5_data <- FetchData(norm.combined, vars = "UBA5") 
hist(uba5_data$rna_UBA5)

# individual level

# HC020
foxg1.0 <- FetchData(norm, vars = "FOXG1") 
norm #6169
length(which(foxg1.0$FOXG1 > 0)) # 3980

# HC021
foxg1.1 <- FetchData(norm1, vars = "FOXG1") 
norm1 #6381
length(which(foxg1.1$FOXG1 > 0)) # 4260

# HC066
foxg1.2 <- FetchData(norm2, vars = "FOXG1") 
norm2 # 6946
length(which(foxg1.2$FOXG1 > 0)) # 3891

# HC067 and HC070
foxg1.34 <- FetchData(norm34, vars = "FOXG1") 
norm34 # 5540
length(which(foxg1.34$FOXG1 > 0)) # 3009








