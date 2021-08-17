#packages
library(DropletUtils)
library(Seurat)
library(dplyr)
library(patchwork)
library(scales)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(gplots)

#Loading the data 
#mild
#Output cellranger
sp54.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/test54/new4_0/54_hv/outs/filtered_feature_bc_matrix")
sp55.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/test_55/new4_0/55_hv/outs/filtered_feature_bc_matrix")
sp57.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/test_57/new4_0/57_hv/outs/filtered_feature_bc_matrix")

#create Seurat object
sp54 <- CreateSeuratObject(counts = out, project = "sp54.data",min.cells = 1, min.genes = 150)
sp55 <- CreateSeuratObject(counts = out, project = "sp55.data",min.cells = 1, min.genes = 150)
sp57 <- CreateSeuratObject(counts = out, project = "sp57.data",min.cells = 1, min.genes = 150)

#merge
alldata <- merge(sp54, c(sp55,sp57), add.cell.ids=c("sp54.data","sp55.data","sp57.data"))

#Calculates the percentage of mitochondrial RNA reads
alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^MT-")

#Plot before metrics QC 
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
VlnPlot(alldata, features = feats, pt.size = 0.1, ncol = 3) + 
     NoLegend()


#Apply filters
#cell filtering parameters for selecting cells based on number of genes/cell, UMI counts/cell, and percent mitochondrial genes according to Wauters et al, 2021. 
alldata.filtered.nCount_RNA <- subset(sp54, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 20 & nCount_RNA > 301)


#Plot after metrics QC
VlnPlot(alldata.filtered.nCount_RNA, features = feats, pt.size = 0.1, ncol = 3) + 
   NoLegend()

#Export .tsv file
counts <- as.matrix(alldata.filtered.nCount_RNA@assays$RNA@counts)
write.csv(counts, file="counts.tsv")