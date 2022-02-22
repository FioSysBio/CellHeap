#Packages
library(DropletUtils)
library(Seurat)
library(dplyr)
library(patchwork)
library(scales)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(MAST)
library(SingleCellExperiment)
library(Matrix)
library(reshape2)




#Loading the data
#severe
#balf
patient_1 <- read.table("/home/Script_SD/CellRanger_output/PRJNA661032/SAMN15967311_v6/SoupX/features_not_infected_severe_SAMN15967311.tsv", header=TRUE,row.names=1)
patient_2 <- read.table("/home/Script_SD/CellRanger_output/PRJNA661032/SAMN15967304_v6/SoupX/features_not_infected_severe_SAMN15967304.tsv", header=TRUE,row.names=1)
patient_3 <- read.table("/home/Script_SD/CellRanger_output/PRJNA661032/SAMN15967297_v6/SoupX/features_not_infected_severe_SAMN15967297.tsv", header=TRUE,row.names=1)
patient_5 <- read.table("/home/Script_SD/CellRanger_output/PRJNA661032/SAMN15967318_v6/SoupX/features_not_infected_severe_SAMN15967318.tsv", header=TRUE,row.names=1)

#Create Seurat object
patient_1 <- CreateSeuratObject(counts = patient_1, project = "patient1")
patient_2 <- CreateSeuratObject(counts = patient_2, project = "patient2")
patient_3 <- CreateSeuratObject(counts = patient_3, project = "patient3")
patient_5 <- CreateSeuratObject(counts = patient_5, project = "patient5")

patient_1$patient1 <- "patient_1"
patient_2$patient1 <- "patient_2"
patient_3$patient1 <- "patient_3"
patient_5$patient1 <- "patient_5"

#Normalization and findVariableFeatures
patient1 <-NormalizeData(patient_1,normalization.method = "LogNormalize", scale.factor = 10000)
patient1 <- FindVariableFeatures(patient1, selection.method = "vst", nfeatures = 2000)
patient2  <-NormalizeData(patient_2,normalization.method = "LogNormalize", scale.factor = 10000)
patient2  <- FindVariableFeatures(patient2, selection.method = "vst", nfeatures = 2000)
patient3 <-NormalizeData(patient_3,normalization.method = "LogNormalize", scale.factor = 10000)
patient3 <- FindVariableFeatures(patient3, selection.method = "vst", nfeatures = 2000)
patient5 <-NormalizeData(patient_5,normalization.method = "LogNormalize", scale.factor = 10000)
patient5 <- FindVariableFeatures(patient5, selection.method = "vst", nfeatures = 2000)


#Integration of two samples
#The parameter dims is based in paper Liao et al,2020.
immune.anchors <- FindIntegrationAnchors(object.list = list(patient1, patient2,patient3,patient5 ),  dims= 1:50)


#Create the integrated assays
immune.combined <- IntegrateData(anchorset = immune.anchors, dims= 1:30)
DefaultAssay(immune.combined) <- "integrated"

#Scaling the dataset
immune.combined <- ScaleData(immune.combined)

#Run PCA clustering and visualization
immune.combined <- RunPCA(immune.combined, npcs = 50)

#Run UMAP clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:50)

#Calculate k.param nearest neighbors for a dataset.
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:50)

#Identify clusters of cells by SNN modularity.
#The parameter resolution is based in paper Liao et al 2020
immune.combined <- FindClusters(immune.combined, resolution = 1.2)


#Find markers for every cluster compared to all remaining cells
immune.combined@misc$markers <- FindAllMarkers(immune.combined, assay = 'RNA',logfc.threshold = 0.25, only.pos = TRUE, test.use = 'MAST')

#Save txt file
write.table(immune.combined@misc$markers,file='marker_MAST.txt',row.names = FALSE,quote = FALSE,sep = '\t')

#Identify differential expressed genes (DEGs)
immune.combined$celltype.orig.ident <- paste(Idents(immune.combined), immune.combined$patient1, sep = "_")

#Create a column in the meta.data slot to hold both the clusters
immune.combined$seurat_clusters <- Idents(immune.combined)

Idents(immune.combined) <- "celltype.orig.ident"

#Rename of clusters (cell subpopulation) after output Enrich from PanglaoDB Augmented 2021 database about markers contained in marker_MAST.txt file
#macrophages subpopulation genes dataset: pSTAT1, TNF, IL6, IL1B, CXCL10, CXCL9, IDO1, IRF5.
immune.combined <- RenameIdents(immune.combined, `10_patient_1` = "Macrophage",`10_patient_2` = "Macrophage", `10_patient_3` = "Macrophage",`10_patient_5` = "Macrophage")


DefaultAssay(immune.combined) <- "RNA"

#Violin plots of the genes of interest of each sample with integrated data
dpi = 300
png(file="violin_fcgr1a_gene.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
VlnPlot(immune.combined, features = "human----FCGR1A", idents = "Macrophage", group.by = "orig.ident")
dev.off()

png(file="violin_fcer1g_gene.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
VlnPlot(immune.combined, features = "human----FCER1G",idents = "Macrophage", group.by = "orig.ident")
dev.off()

png(file="violin_hck_gene.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
VlnPlot(immune.combined, features = "human----HCK", idents = "Macrophage", group.by = "orig.ident")
dev.off()

png(file="violin_syk_gene.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
VlnPlot(immune.combined, features = "human----SYK", idents = "Macrophage", group.by = "orig.ident")
dev.off()

png(file="violin_raf1_gene.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
VlnPlot(immune.combined, features = "human----RAF1", idents = "Macrophage", group.by = "orig.ident")
dev.off()

png(file="violin_mapk8_gene.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
VlnPlot(immune.combined, features = "human----MAPK8", idents = "Macrophage", group.by = "orig.ident")
dev.off()

png(file="violin_mapk11_gene.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
VlnPlot(immune.combined, features = "human----MAPK11", idents = "Macrophage", group.by = "orig.ident")
dev.off()

png(file="violin_mapk1_gene.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
VlnPlot(immune.combined, features = "human----MAPK1", idents = "Macrophage", group.by = "orig.ident")
dev.off()





