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

#Violin plot before filters
png("/scratch/inova-covd19/vanessa.silva/control_before_filters_liao.png")
VlnPlot(alldata, features = feats, pt.size = 0.1, ncol = 3) + 
+      NoLegend()
dev.off()

#Apply filters
# cell filtering parameters for selecting cells based on number of genes/cell, UMI counts/cell, and percent mitochondrial genes according to Liao et al, 2020. 
alldata.new <- subset(alldata, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 0.1 & nCount_RNA > 1000)

#Violin plot after filters
png("/scratch/inova-covd19/vanessa.silva/control_after_filters_liao.png")
VlnPlot(alldata.new, features = feats, pt.size = 0.1, ncol = 3) + 
     NoLegend()
dev.off()

#Export .tsv and .csv files
counts <- as.matrix(alldata.new@assays$RNA@counts)
write.csv(counts, file="features_control_liao.tsv")

#import of seurat filter file
sample <- read.table("features_control_liao.tsv", sep=",", header=T, row.names = 1)

#select rows containing sarscov2
sample_sarscov2 <- sample[grep("^sarscov2", row.names(sample)),]

#transpose data frame
sample_sarscov2 <- t(sample_sarscov2)

#print rows where all columns are zero
not_infected <- rownames(sample_sarscov2[ rowSums(sample_sarscov2)==0, ] )

#print rows where all columns are different of zero
infected <- rownames(sample_sarscov2[ rowSums(sample_sarscov2)>0, ] )

#### Export data #####

#dataframe not infected 
write.table(sample[,not_infected], file= "features_not_infected_control_liao.tsv", quote=FALSE, sep=',', col.names = TRUE)
