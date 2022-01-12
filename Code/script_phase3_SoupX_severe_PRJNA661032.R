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

SAMN15967269_severe_SpX <- load10X("/home/Script_SD/CellRanger_output/PRJNA661032/SAMN15967269_v6/outs")

#Estima o parâmetro "rho" que representa a fração de contaminação, onde rho igual a 0 significa nenhuma contaminação e rho igual a 1 significa 100% de UMIs em um droplet são "soup"

SAMN15967269_severe_SpX <- autoEstCont(SAMN15967269_severe_SpX)

#Limpar os dados, ou seja gerar a matriz corrigida sem contaminação

out_SAMN15967269_severe_SpX <- adjustCounts(SAMN15967269_severe_SpX)

#Caso queira salvar a nova matriz corrigida no formato do CellRanger, deve utilizar o seguinte comando

#DropletUtils:::write10xCounts("./SoupX_SAMN15967269_severe", SAMN15967269_severe_SpX)

#Pode utilizar as matrizes corrigidas provenientes do SoupX diretamente no Seurat através da função CreateSeuratObject

SAMN15967269_severe <- CreateSeuratObject(counts = out_SAMN15967269_severe_SpX, project = "PRJNA661032_Severe",min.cells = 1)

#Discutir qual a melhor abordagem e se seria possível fazer o merge de todas as amostras do mesmo grupo


#Calculates the percentage of mitochondrial RNA reads
SAMN15967269_severe[["percent.mt"]] <- PercentageFeatureSet(SAMN15967269_severe, pattern = "^MT-")

#Métricas
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

#Etapa dos filtros de controle de qualidade definidas nas análises anteriores baseadas no Wauters et al, 2021
# cell filtering parameters for selecting cells based on number of genes/cell, UMI counts/cell, and percent mitochondrial genes according to Wauters et al, 2021.
SAMN15967269_severe.filtered.nCount_RNA <- subset(SAMN15967269_severe, subset = nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 20 & nCount_RNA > 301)

#As próximas etapas envolvem a busca pelo RNA viral nas amostras, apenas mantive os comandos do script anterior de fase 3
#Export .tsv and .csv files
out_test <- as.matrix(SAMN15967269_severe.filtered.nCount_RNA@assays$RNA@counts)
write.csv(out_test, file="features_severe_SAMN15967269.tsv")
write.csv(SAMN15967269_severe.filtered.nCount_RNA@meta.data, file="metadata_severe_SAMN15967269.csv")

#import of seurat filter file
sample <- read.table("features_severe_SAMN15967269.tsv", sep=",", header=T, row.names = 1)

#select rows containing sarscov2
sample_sarscov2 <- sample[grep("^sarscov2", row.names(sample)),]

#transpose data frame
sample_sarscov2 <- t(sample_sarscov2)

#print rows where all columns are zero
not_infected <- rownames(sample_sarscov2[ rowSums(sample_sarscov2)==0, ] )

#print rows where all columns are different of zero
infected <- rownames(sample_sarscov2[ rowSums(sample_sarscov2)>0, ] )

#dataframe not infected
write.table(sample[,not_infected], file= "features_not_infected_severe_SAMN15967269.tsv", quote=FALSE, sep=',', col.names = TRUE)

#dataframe  infected
write.table(sample[,infected], file= "features_infected_severe_SAMN15967269.tsv", quote=FALSE, sep=',', col.names = TRUE)
