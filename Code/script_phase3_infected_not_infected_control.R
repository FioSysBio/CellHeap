#packages
library(Seurat)
library(dplyr)
library(patchwork)
library(scales)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(gplots)

#Loading the data 
#control
ctrl1.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/control/control1_hv/outs/filtered_feature_bc_matrix")
ctrl2.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/control/control2/control2_hv/outs/filtered_feature_bc_matrix")
ctrl3.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/control/control3/control3_hv/outs/filtered_feature_bc_matrix")

#create Seurat object
#control
ctrl1 <- CreateSeuratObject(counts = ctrl1.data, project = "control",min.cells = 1, min.genes = 150)
ctrl2 <- CreateSeuratObject(counts = ctrl2.data, project = "control",min.cells = 1, min.genes = 150)
ctrl3 <- CreateSeuratObject(counts = ctrl3.data, project = "control",min.cells = 1, min.genes = 150)

#merge files
alldata <- merge(ctrl1, c(ctrl2,ctrl3), add.cell.ids=c("ctrl1.data","ctrl2.data","ctrl3.data"))

#Calculates the percentage of mitochondrial RNA reads
alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^MT-")

#metrics QC
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

#Apply filters
# cell filtering parameters for selecting cells based on number of genes/cell, UMI counts/cell, and percent mitochondrial genes according to Wauters et al, 2021. 
alldata.filtered.nCount_RNA <- subset(alldata, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 20 & nCount_RNA > 301)

#Export .tsv and .csv files
counts <- as.matrix(alldata.filtered.nCount_RNA@assays$RNA@counts)
write.csv(counts, file="features_control.tsv")
write.csv(alldata.filtered.nCount_RNA@meta.data, file="metadata_control.csv")

#import of seurat filter file
sample <- read.table("features_control.tsv", sep=",", header=T, row.names = 1)

#select rows containing sarscov2
sample_sarscov2 <- sample[grep("^sarscov2", row.names(sample)),]

#transpose data frame
sample_sarscov2 <- t(sample_sarscov2)

#print rows where all columns are zero
not_infected <- rownames(sample_sarscov2[ rowSums(sample_sarscov2)==0, ] )

#print rows where all columns are different of zero
infected <- rownames(sample_sarscov2[ rowSums(sample_sarscov2)>0, ] )

#dataframe not infected
write.table(sample[,not_infected], file= "features_not_infected_control.tsv", quote=FALSE, sep=',', col.names = TRUE)

#dataframe  infected
write.table(sample[,infected], file= "features_infected_control.tsv", quote=FALSE, sep=',', col.names = TRUE)

#import dataframe not infected
not_virus <- read.csv("features_not_infected_control.tsv", header=TRUE,row.names=1)

#create Seurat object only not infected
not_infected_sample <- CreateSeuratObject(counts = not_virus)
