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
#mild
not_virus_mild <- read.csv("features_not_infected_mild_liao.tsv", header=TRUE,row.names=1)

#Create Seurat object
all_mild <- CreateSeuratObject(counts = not_virus_mild, project = "mild")
all_mild$control <- "mild"

#Print dataset size
table(Idents(all_mild))

#Normalization and findVariableFeatures
all_mild <-NormalizeData(all_mild,normalization.method = "LogNormalize", scale.factor = 10000)
all_mild <- FindVariableFeatures(all_mild, selection.method = "vst", nfeatures = 2000)

#control
not_virus_control <- read.csv("features_not_infected_control_liao.tsv", header=TRUE,row.names=1)

#Create Seurat object
all_control <- CreateSeuratObject(counts = not_virus_control, project = "control")
all_control$control <- "control"

#Print dataset size
table(Idents(all_control))

#Normalization and findVariableFeatures
#The parameters are based in paper Liao et al,2020.
all_control  <-NormalizeData(all_control,normalization.method = "LogNormalize", scale.factor = 10000)
all_control  <- FindVariableFeatures(all_control, selection.method = "vst", nfeatures = 2000)

#Integration of two samples
#The parameter dims is based in paper Liao et al,2020.
immune.anchors <- FindIntegrationAnchors(object.list = list(all_control, all_mild),  dims= 1:50)


#Create the integrated assays
immune.combined <- IntegrateData(anchorset = immune.anchors, dims= 1:50)
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

#Violin plots of the genes of interest of each sample with integrated data
dpi = 300
png(file="violin_fcn1_gene_integ_liao.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
VlnPlot(immune.combined, features = "human----FCN1", group.by = "orig.ident")
dev.off()

png(file="violin_spp1_gene_integ_liao.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
VlnPlot(immune.combined, features = "human----SPP1", group.by = "orig.ident")
dev.off()

png(file="violin_fabp4_gene_integ_liao.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
VlnPlot(immune.combined, features = "human----FABP4", group.by = "orig.ident")
dev.off()

#RNA data
DefaultAssay(immune.combined) <- "RNA"

#Find markers for every cluster compared to all remaining cells
immune.combined@misc$markers <- FindAllMarkers(immune.combined, assay = 'RNA',only.pos = TRUE, test.use = 'MAST')

#Identify differential expressed genes (DEGs) across conditions
immune.combined$celltype.orig.ident <- paste(Idents(immune.combined), immune.combined$control, sep = "_")

#Create a column in the meta.data slot to hold both the clusters
immune.combined$seurat_clusters <- Idents(immune.combined)

Idents(immune.combined) <- "celltype.orig.ident"
 
#Rename of clusters (cell subpopulation) after output Enrich from PanglaoDB Augmented 2021 database about markers contained in marker_MAST.txt file
immune.combined <- RenameIdents(immune.combined, `3_mild` = "Macrophage_mild", `3_control` = "Macrophage_control")


#Find the genes that are different between conditions
interferon.response <- FindMarkers(immune.combined, ident.1 = "Macrophage_mild", ident.2 = "Macrophage_control", logfc.threshold = 0.25, verbose = FALSE, test.use = "MAST", only.pos = FALSE)

#Violin plots of the genes of interest of each sample - RNA data
dpi = 300
png(file="violin_fcn1_gene_liao.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
VlnPlot(immune.combined, features = "human----FCN1", group.by = "orig.ident")
dev.off()

png(file="violin_spp1_gene_liao.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
VlnPlot(immune.combined, features = "human----SPP1", group.by = "orig.ident")
dev.off()

png(file="violin_fabp4_gene_liao.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
VlnPlot(immune.combined, features = "human----FABP4", group.by = "orig.ident")
dev.off()

