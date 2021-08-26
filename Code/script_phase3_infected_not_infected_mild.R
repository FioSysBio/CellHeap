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

#create Seurat object
#mild
sp54 <- CreateSeuratObject(counts = sp54.data, project = "mild",min.cells = 1, min.genes = 150)
sp55 <- CreateSeuratObject(counts = sp55.data, project = "mild",min.cells = 1, min.genes = 150)
sp57 <- CreateSeuratObject(counts = sp57.data, project = "mild",min.cells = 1, min.genes = 150)

#merge files
alldata <- merge(sp54, c(sp55, sp57), add.cell.ids=c("sp54.data","sp55.data","sp57.data"))

#Calculates the percentage of mitochondrial RNA reads
alldata[["percent.mt"]] <- PercentageFeatureSet(alldata, pattern = "^MT-")

#metrics QC
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

#Apply filters
# cell filtering parameters for selecting cells based on number of genes/cell, UMI counts/cell, and percent mitochondrial genes according to Wauters et al, 2021. 
alldata.filtered.nCount_RNA <- subset(alldata, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 20 & nCount_RNA > 301)

#Export .tsv and .csv files
counts <- as.matrix(not_infected_sample@assays$RNA@counts)
write.csv(counts, file="features_mild.tsv")
write.csv(alldata.filtered.nCount_RNA@meta.data, file="metadata_mild.csv")

#import of seurat filter file
sample <- read.table("features_mild.tsv", sep=",", header=T, row.names = 1)

#select rows containing sarscov2
sample_sarscov2 <- sample[grep("^sarscov2", row.names(sample)),]

#transpose data frame
sample_sarscov2 <- t(sample_sarscov2)

#print rows where all columns are zero
not_infected <- rownames(sample_sarscov2[ rowSums(sample_sarscov2)==0, ] )

#print rows where all columns are different of zero
infected <- rownames(sample_sarscov2[ rowSums(sample_sarscov2)>0, ] )

#dataframe not infected
write.table(sample[,not_infected], file= "features_not_infected_mild.tsv", quote=FALSE, sep=',', col.names = TRUE)

#dataframe  infected
write.table(sample[,infected], file= "features_infected_mild.tsv", quote=FALSE, sep=',', col.names = TRUE)
