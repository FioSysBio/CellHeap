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
#mild
#Output cellranger
sp54.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/test54/new4_0/54_hv/outs/filtered_feature_bc_matrix")
sp55.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/test_55/new4_0/55_hv/outs/filtered_feature_bc_matrix")
sp57.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/test_57/new4_0/57_hv/outs/filtered_feature_bc_matrix")

#severe
sp56.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/test_56/new4_0/56_hv/outs/filtered_feature_bc_matrix")
sp58.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/test_58/new4_0/58_hv/outs/filtered_feature_bc_matrix")
sp49.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/test_49/new4_0/49_hv/outs/filtered_feature_bc_matrix")
sp50.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/test50/new4_0/50_hv/outs/filtered_feature_bc_matrix")
sp51.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/test_51/new4_0/51_hv/outs/filtered_feature_bc_matrix")

#control
ctrl1.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/control/control1_hv/outs/filtered_feature_bc_matrix")
ctrl2.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/control/control2/control2_hv/outs/filtered_feature_bc_matrix")
ctrl3.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/control/control3/control3_hv/outs/filtered_feature_bc_matrix")

#create Seurat object
#mild
sp54 <- CreateSeuratObject(counts = sp54.data, project = "mild",min.cells = 1, min.genes = 150)
sp55 <- CreateSeuratObject(counts = sp55.data, project = "mild",min.cells = 1, min.genes = 150)
sp57 <- CreateSeuratObject(counts = sp57.data, project = "mild",min.cells = 1, min.genes = 150)

#severe
sp56 <- CreateSeuratObject(counts = sp56.data, project = "severe",min.cells = 1, min.genes = 150)
sp58 <- CreateSeuratObject(counts = sp58.data, project = "severe",min.cells = 1, min.genes = 150)
sp59 <- CreateSeuratObject(counts = sp59.data, project = "severe",min.cells = 1, min.genes = 150)
sp49 <- CreateSeuratObject(counts = sp49.data, project = "severe",min.cells = 1, min.genes = 150)
sp50 <- CreateSeuratObject(counts = sp50.data, project = "severe",min.cells = 1, min.genes = 150)
sp51 <- CreateSeuratObject(counts = sp51.data, project = "severe",min.cells = 1, min.genes = 150)

#control
ctrl1 <- CreateSeuratObject(counts = ctrl1.data, project = "control",min.cells = 1, min.genes = 150)
ctrl2 <- CreateSeuratObject(counts = ctrl2.data, project = "control",min.cells = 1, min.genes = 150)
ctrl3 <- CreateSeuratObject(counts = ctrl3.data, project = "control",min.cells = 1, min.genes = 150)

#merge files
alldata <- merge(sp54, c(sp55, sp57,sp56, sp58, sp59, sp49,sp50,sp51, ctrl1,ctrl2,ctrl3), add.cell.ids=c("sp54.data","sp55.data","sp57.data","sp56.data","sp58.data","sp59.data","sp49.data","sp50.data","sp51.data","ctrl1.data","ctrl2.data","ctrl3.data"))

#Calculates the percentage of mitochondrial RNA reads
alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^MT-")

#metrics QC
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

#Apply filters
# cell filtering parameters for selecting cells based on number of genes/cell, UMI counts/cell, and percent mitochondrial genes according to Wauters et al, 2021. 
alldata.filtered.nCount_RNA <- subset(alldata, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 20 & nCount_RNA > 301)

#Export .tsv file
counts <- as.matrix(alldata.filtered.nCount_RNA@assays$RNA@counts)
write.csv(counts, file="features.tsv")

#import of seurat filter file
sample <- read.table("features.tsv", sep=",", header=T, row.names = 1)

#select rows containing sarscov2
sample_sarscov2 <- sample[grep("^sarscov2", row.names(sample)),]

#transpose data frame
sample_sarscov2 <- t(sample_sarscov2)

#print rows where all columns are zero
not_infected <- rownames(sample_sarscov2[ rowSums(sample_sarscov2)==0, ] )

#print rows where all columns are different of zero
infected <- rownames(sample_sarscov2[ rowSums(sample_sarscov2)>0, ] )

#dataframe not infected
write.table(sample[,not_infected], file= "features_not_infected.tsv", quote=FALSE, sep=',', col.names = TRUE)

#dataframe  infected
write.table(sample[,infected], file= "features_infected.tsv", quote=FALSE, sep=',', col.names = TRUE)

#import dataframe not infected
not_virus <- read.csv("features_not_infected", header=TRUE,row.names=1)

#create Seurat object only not infected
not_infected_sample <- CreateSeuratObject(counts = not_virus)



