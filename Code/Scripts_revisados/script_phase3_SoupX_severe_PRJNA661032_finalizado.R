setwd("/home/Script_SD/CellRanger_output/PRJNA661032/SAMN15967319_v6/SoupX")
#packages
library(Seurat)
library(dplyr)
library(patchwork)
library(scales)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(SoupX)
library(DropletUtils)
library(conflicted)

#SoupX - Modo automático usando os outputs do CellRanger...o SoupX diz que os resultados serão melhores caso tenha alguma informações de clustering básico e diz que "If you are using 10X data mapped with cellranger, the default clustering produced by cellranger is automatically loaded and used"

#Carrega os dados dos outputs diretos do CellRanger, cria o objeto SoupChannel e estima o soup profile

SAMN15967319_severe_SpX <- load10X("/home/Script_SD/CellRanger_output/PRJNA661032/SAMN15967319_v6/outs")

#Estima o parâmetro "rho" que representa a fração de contaminação, onde rho igual a 0 significa nenhuma contaminação e rho igual a 1 significa 100% de UMIs em um droplet são "soup"

SAMN15967319_severe_SpX <- autoEstCont(SAMN15967319_severe_SpX, doPlot=FALSE)

#Limpar os dados, ou seja gerar a matriz corrigida sem contaminação

out_SAMN15967319_severe_SpX <- adjustCounts(SAMN15967319_severe_SpX)

#Caso queira salvar a nova matriz corrigida no formato do CellRanger, deve utilizar o seguinte comando

#DropletUtils:::write10xCounts("./SoupX_SAMN15967319_severe", SAMN15967319_severe_SpX)

#Pode utilizar as matrizes corrigidas provenientes do SoupX diretamente no Seurat através da função CreateSeuratObject

SAMN15967319_severe <- CreateSeuratObject(counts = out_SAMN15967319_severe_SpX, project = "PRJNA661032_Severe",min.cells = 1)
write.csv(SAMN15967319_severe@meta.data, file="metadata_severe_SAMN15967319.csv")

#Discutir qual a melhor abordagem e se seria possível fazer o merge de todas as amostras do mesmo grupo

#Calculates the percentage of mitochondrial RNA reads
SAMN15967319_severe[["percent.mt"]] <- PercentageFeatureSet(SAMN15967319_severe, pattern = "^human----MT-")

#Métricas
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

png("severe_before_filters_SAMN15967319.png")

VlnPlot(SAMN15967319_severe, features = feats, pt.size = 0.1, ncol = 3) +
     NoLegend()
dev.off()

#Etapa dos filtros de controle de qualidade definidas nas análises anteriores baseadas no Wauters et al, 2021
# cell filtering parameters for selecting cells based on number of genes/cell, UMI counts/cell, and percent mitochondrial genes according to Wauters et al, 2021.
SAMN15967319_severe_filtered <- subset(SAMN15967319_severe, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 20 & nCount_RNA > 301)

png("severe_after_filters_SAMN15967319.png")

VlnPlot(SAMN15967319_severe_filtered, features = feats, pt.size = 0.1, ncol = 3) +
     NoLegend()
dev.off()

#As próximas etapas envolvem a busca pelo RNA viral nas amostras, apenas mantive os comandos do script anterior de fase 3
#Export .tsv and .csv files
out_test <- as.matrix(SAMN15967319_severe_filtered@assays$RNA@counts)
write.table(out_test, file="features_severe_SAMN15967319.tsv", quote=FALSE, sep='\t', col.names = TRUE)

#import of seurat filter file
sample <- read.table("features_severe_SAMN15967319.tsv", sep="\t", header=T, row.names = 1)

#select rows containing sarscov2
sample_sarscov2 <- sample[grep("virus-v6", row.names(sample)),,drop=FALSE]

#transpose data frame
sample_sarscov2_transposta <- t(sample_sarscov2)

#print rows where all columns are zero - CORRIGIDO
not_infected <- rownames(sample_sarscov2_transposta)[which(rowSums(sample_sarscov2_transposta)==0)]

#print rows where all columns are different of zero - CORRIGIDO
infected <- rownames(sample_sarscov2_transposta)[which(rowSums(sample_sarscov2_transposta)>0)]

#dataframe not infected - CORRIGIDO
write.table(sample[,not_infected,drop=FALSE], file= "features_not_infected_severe_SAMN15967319.tsv", quote=FALSE, sep='\t', col.names = TRUE)

#dataframe  infected - CORRIGIDO
write.table(sample[,infected,drop=FALSE], file= "features_infected_severe_SAMN15967319.tsv", quote=FALSE, sep='\t', col.names = TRUE)
