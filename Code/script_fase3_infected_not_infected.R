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
#Output cellranger
sp56.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/test_56/new4_0/56_hv/outs/filtered_feature_bc_matrix")
sp58.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/test_58/new4_0/58_hv/outs/filtered_feature_bc_matrix")
sp49.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/test_49/new4_0/49_hv/outs/filtered_feature_bc_matrix")
sp50.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/test50/new4_0/50_hv/outs/filtered_feature_bc_matrix")
sp51.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/test_51/new4_0/51_hv/outs/filtered_feature_bc_matrix")
spbal1.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal01/bal01_hv/outs/filtered_feature_bc_matrix")
spbal2.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal02/bal02_hv/outs/filtered_feature_bc_matrix")
spbal3.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal03/bal03_hv/outs/filtered_feature_bc_matrix")
spbal4.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal04/bal04_hv/outs/filtered_feature_bc_matrix")
spbal5.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal05/bal05_hv/outs/filtered_feature_bc_matrix")
spbal6.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal06/bal06_hv/outs/filtered_feature_bc_matrix")
spbal7.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal07/bal07_hv/outs/filtered_feature_bc_matrix")
spbal8.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal08/bal08_hv/outs/filtered_feature_bc_matrix")
spbal9.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal09/bal09_hv/outs/filtered_feature_bc_matrix")
spbal10.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal10/bal10_hv/outs/filtered_feature_bc_matrix")
spbal11.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal11/bal11_hv/outs/filtered_feature_bc_matrix")
spbal14.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal14/bal14_hv/outs/filtered_feature_bc_matrix")
spbal16.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal16/bal16_hv/outs/filtered_feature_bc_matrix")
spbal17.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal17/bal17_hv/outs/filtered_feature_bc_matrix")
spbal18.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal18/bal18_hv/outs/filtered_feature_bc_matrix")
spbal22.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal22/bal22_hv/outs/filtered_feature_bc_matrix")
spbal23.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal23/bal23_hv/outs/filtered_feature_bc_matrix")
spbal24.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal24/bal24_hv/outs/filtered_feature_bc_matrix")
spbal25.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal25/bal25_hv/outs/filtered_feature_bc_matrix")
spbal27.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal27/bal27_hv/outs/filtered_feature_bc_matrix")
spbal29.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/bal29/bal29_hv/outs/filtered_feature_bc_matrix")
spsar01.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar01/sar01_hv/outs/filtered_feature_bc_matrix")
spsar02.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar02/sar02_hv/outs/filtered_feature_bc_matrix")
spsar03.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar03/sar03_hv/outs/filtered_feature_bc_matrix")
spsar04.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar04/sar04_hv/outs/filtered_feature_bc_matrix")
spsar05.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar05/sar05_hv/outs/filtered_feature_bc_matrix")
spsar06.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar06/sar06_hv/outs/filtered_feature_bc_matrix")
spsar07.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar07/sar07_hv/outs/filtered_feature_bc_matrix")
spsar08.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar08/sar08_hv/outs/filtered_feature_bc_matrix")
spsar09.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar09/sar09_hv/outs/filtered_feature_bc_matrix")
spsar10.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar10/sar10_hv/outs/filtered_feature_bc_matrix")
spsar11.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar11/sar11_hv/outs/filtered_feature_bc_matrix")
spsar12.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar12/sar12_hv/outs/filtered_feature_bc_matrix")
spsar13.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar13/sar13_hv/outs/filtered_feature_bc_matrix")
spsar14.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar14/sar14_hv/outs/filtered_feature_bc_matrix")
spsar15.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar15/sar15_hv/outs/filtered_feature_bc_matrix")
spsar16.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar16/sar16_hv/outs/filtered_feature_bc_matrix")
spsar17.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar17/sar17_hv/outs/filtered_feature_bc_matrix")
spsar18.data <- Read10X(data.dir = "/scratch/inova-covd19/vanessa.silva/sar18/sar18_hv/outs/filtered_feature_bc_matrix")

#create Seurat object
sp56 <- CreateSeuratObject(counts = out, project = "sp56.data",min.cells = 1, min.genes = 150)
sp58 <- CreateSeuratObject(counts = sp58.data, project = "sp58.data",min.cells = 1, min.genes = 150)
sp59 <- CreateSeuratObject(counts = sp59.data, project = "sp59.data",min.cells = 1, min.genes = 150)
sp49 <- CreateSeuratObject(counts = sp49.data, project = "sp49.data",min.cells = 1, min.genes = 150)
sp50 <- CreateSeuratObject(counts = sp50.data, project = "sp50.data",min.cells = 1, min.genes = 150)
sp51 <- CreateSeuratObject(counts = sp51.data, project = "sp51.data",min.cells = 1, min.genes = 150)
spbal1 <- CreateSeuratObject(counts = spbal1.data, project = "spbal1.data",min.cells = 1, min.genes = 150)
spbal2 <- CreateSeuratObject(counts = spbal2.data, project = "spbal2.data",min.cells = 1, min.genes = 150)
spbal3 <- CreateSeuratObject(counts = spbal3.data, project = "spbal3.data",min.cells = 1, min.genes = 150)
spbal4 <- CreateSeuratObject(counts = spbal4.data, project = "spbal4.data",min.cells = 1, min.genes = 150)
spbal5 <- CreateSeuratObject(counts = spbal5.data, project = "spbal5.data",min.cells = 1, min.genes = 150)
spbal6 <- CreateSeuratObject(counts = spbal6.data, project = "spbal6.data",min.cells = 1, min.genes = 150)
spbal7 <- CreateSeuratObject(counts = spbal7.data, project = "spbal7.data",min.cells = 1, min.genes = 150)
spbal8 <- CreateSeuratObject(counts = spbal8.data, project = "spbal8.data",min.cells = 1, min.genes = 150)
spbal9 <- CreateSeuratObject(counts = spbal9.data, project = "spbal9.data",min.cells = 1, min.genes = 150)
spbal10 <- CreateSeuratObject(counts = spbal10.data, project = "spbal10.data",min.cells = 1, min.genes = 150)
spbal11 <- CreateSeuratObject(counts = spbal11.data, project = "spbal11.data",min.cells = 1, min.genes = 150)
spbal14 <- CreateSeuratObject(counts = spbal14.data, project = "spbal14.data",min.cells = 1, min.genes = 150)
spbal16 <- CreateSeuratObject(counts = spbal16.data, project = "spbal16.data",min.cells = 1, min.genes = 150)
spbal17 <- CreateSeuratObject(counts = spbal17.data, project = "spbal17.data",min.cells = 1, min.genes = 150)
spbal18 <- CreateSeuratObject(counts = spbal18.data, project = "spbal18.data",min.cells = 1, min.genes = 150)
spbal22 <- CreateSeuratObject(counts = spbal22.data, project = "spbal22.data",min.cells = 1, min.genes = 150)
spbal23 <- CreateSeuratObject(counts = spbal23.data, project = "spbal23.data",min.cells = 1, min.genes = 150)
spbal24 <- CreateSeuratObject(counts = spbal24.data, project = "spbal24.data",min.cells = 1, min.genes = 150)
spbal25 <- CreateSeuratObject(counts = spbal25.data, project = "spbal25.data",min.cells = 1, min.genes = 150)
spbal27 <- CreateSeuratObject(counts = spbal27.data, project = "spbal27.data",min.cells = 1, min.genes = 150)
spbal29 <- CreateSeuratObject(counts = spbal29.data, project = "spbal29.data",min.cells = 1, min.genes = 150)
spsar01 <- CreateSeuratObject(counts = spsar01.data, project = "spsar01.data",min.cells = 1, min.genes = 150)
spsar02 <- CreateSeuratObject(counts = spsar02.data, project = "spsar02.data",min.cells = 1, min.genes = 150)
spsar03 <- CreateSeuratObject(counts = spsar03.data, project = "spsar03.data",min.cells = 1, min.genes = 150)
spsar04 <- CreateSeuratObject(counts = spsar04.data, project = "spsar04.data",min.cells = 1, min.genes = 150)
spsar05 <- CreateSeuratObject(counts = spsar05.data, project = "spsar05.data",min.cells = 1, min.genes = 150)
spsar06 <- CreateSeuratObject(counts = spsar06.data, project = "spsar06.data",min.cells = 1, min.genes = 150)
spsar07 <- CreateSeuratObject(counts = spsar07.data, project = "spsar07.data",min.cells = 1, min.genes = 150)
spsar08 <- CreateSeuratObject(counts = spsar08.data, project = "spsar08.data",min.cells = 1, min.genes = 150)
spsar09 <- CreateSeuratObject(counts = spsar09.data, project = "spsar09.data",min.cells = 1, min.genes = 150)
spsar10 <- CreateSeuratObject(counts = spsar10.data, project = "spsar10.data",min.cells = 1, min.genes = 150)
spsar11 <- CreateSeuratObject(counts = spsar11.data, project = "spsar11.data",min.cells = 1, min.genes = 150)
spsar12 <- CreateSeuratObject(counts = spsar12.data, project = "spsar12.data",min.cells = 1, min.genes = 150)
spsar13 <- CreateSeuratObject(counts = spsar13.data, project = "spsar13.data",min.cells = 1, min.genes = 150)
spsar14 <- CreateSeuratObject(counts = spsar14.data, project = "spsar14.data",min.cells = 1, min.genes = 150)
spsar15 <- CreateSeuratObject(counts = spsar15.data, project = "spsar15.data",min.cells = 1, min.genes = 150)
spsar16 <- CreateSeuratObject(counts = spsar16.data, project = "spsar16.data",min.cells = 1, min.genes = 150)
spsar17 <- CreateSeuratObject(counts = spsar17.data, project = "spsar17.data",min.cells = 1, min.genes = 150)
spsar18 <- CreateSeuratObject(counts = spsar18.data, project = "spsar18.data",min.cells = 1, min.genes = 150)

#merge files
alldata <- merge(sp56, c(sp58, sp59, sp49,sp50,sp51,spbal1,spbal2,spbal3,spbal4,spbal5,spbal6, spbal7,spbal8, spbal9,spbal10,spbal11,spbal14,spbal16,spbal17, spbal18, spbal22, spbal23, spbal24,spbal25,spbal27,spbal29,spsar01,spsar02,spsar03,spsar04,spsar05,spsar06,spsar07,spsar08,spsar09,spsar10,spsar11,spsar12,spsar13,spsar14,spsar15,spsar16,spsar17, spsar18), add.cell.ids=c("sp56.data","sp58.data","sp59.data","sp49.data","sp50.data","sp51.data", "spbal1.data","spbal2.data","spbal3.data","spbal4.data","spbal05","spbal6.data", "spbal7.data", "spbal8.data", "spbal9.data", "spbal10.data", "spbal11.data", "spbal14.data", "spbal16.data", "spbal17.data", "spbal18.data", "spbal22.data", "spbal23.data", "spbal24.data","spbal25.data", "spbal27.data",  "spbal29.data","spsar01.data","spsar02.data", "spsar03.data", "spsar04.data","spsar05.data","spsar06.data","spsar07.data","spsar08.data","spsar09.data", "spsar10.data", "spsar11.data","spsar12.data","spsar13.data", "spsar14.data","spsar15.data","spsar16.data", "spsar17.data", "spsar18.data" ))

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



