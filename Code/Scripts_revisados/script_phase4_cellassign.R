#Libraries.
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
library(cellassign)
#libraries suggested by Gabriel.
#human database library.
#https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html
library(org.Hs.eg.db)
library(edgeR)
library(scater)
library(loomR)

#Samples identified with the ID.
#Samples are BALF only from severe cases of uninfected cells.
#Load the uninfected features file reported in phase 3, indicating the separator that is the tab.
SAMN15967304_sev_balf <- read.table("/home/Script_SD/CellRanger_output/PRJNA661032/SAMN15967304_v6/SoupX/features_not_infected_severe_SAMN15967304.tsv", header=TRUE, sep = "\t")

#Create the Seurat object with the sample ID.
SAMN15967304_sev_balf <- CreateSeuratObject(counts = SAMN15967304_sev_balf, project = "SAMN15967304_sev_balf")
test2.sce <- as.SingleCellExperiment(SAMN15967304_sev_balf)
print("imprimindo_test2.sce")
print(test2.sce)
print("final_impressao_test2.sce")

#Normalization and findVariableFeatures.
#The selection of methods for normalization and features were based on the Seurat vignette.
SAMN15967304_sev_balf <-NormalizeData(SAMN15967304_sev_balf,normalization.method = "LogNormalize", scale.factor = 10000)
SAMN15967304_sev_balf <- FindVariableFeatures(SAMN15967304_sev_balf, selection.method = "vst", nfeatures = 2000)

#ScaleData() defaults to only scaling the highly variable genes identified earlier (2000), which will not affect PCA and clustering results.
#Scaling the dataset
SAMN15967304_sev_balf <- ScaleData(SAMN15967304_sev_balf)

#Run PCA clustering and visualization.
#The parameter npcs = 50 is standard and found in the Seurat vignette.
SAMN15967304_sev_balf <- RunPCA(SAMN15967304_sev_balf, npcs = 50)

#The dimension file was tested and we found that from 16 it is already stable, so we set the dimension to 20.
#The number of replicates was based on the Seurat vignette.
SAMN15967304_sev_balf <- JackStraw(SAMN15967304_sev_balf, num.replicate = 100)
SAMN15967304_sev_balf <- ScoreJackStraw(SAMN15967304_sev_balf, dims = 1:20)

#Dimension test charts.
JackStrawPlot(SAMN15967304_sev_balf, dims = 1:20)
ElbowPlot(SAMN15967304_sev_balf)

#The UMAP to visualize and explore the dataset. As input to UMAP, Seurat suggests using the same PCs as input to cluster analysis
#Run UMAP clustering.
SAMN15967304_sev_balf <- RunUMAP(SAMN15967304_sev_balf, reduction = "pca", dims = 1:16)

#Calculate k.param nearest neighbors for a dataset.
# a função FindNeighbors() leva como entrada a dimensionalidade do conjunto de dados definida anteriormente para a construção do grafo de KNN

SAMN15967304_sev_balf <- FindNeighbors(SAMN15967304_sev_balf, reduction = "pca", dims = 1:16)

#Identify clusters of cells by SNN modularity.
#The parameter resolution is based in paper Liao et al 2020.
#To group the cells, Seurat applies modularity optimization techniques with several algorithms implemented (the default is Louvain and the other options are: Louvain algorithm with multilevel refinement; SLM algorithm and Leiden).
#The resolution parameter defines the "granularity" of clustering with larger values leading to a larger number of clusters. According to Seurat's vignette, setting this parameter between 0.4 - 1.2 usually returns good results for datasets of around 3K cells. The optimal resolution usually increases for larger data sets. The paper Silvin et al 2020 (case study) uses the resolution of 0.3 (I think it's important to test with 0.5).
SAMN15967304_sev_balf <- FindClusters(SAMN15967304_sev_balf, resolution = 1.2)

print("imprimindo objecto Seurat")
print(SAMN15967304_sev_balf)

#convert seurat object to SingleCellExperiment object.
test.sce <- as.SingleCellExperiment(SAMN15967304_sev_balf)


#Find markers for every cluster compared to all remaining cells
SAMN15967304_sev_balf <- FindAllMarkers(SAMN15967304_sev_balf, assay = 'RNA',logfc.threshold = 0.25, only.pos = TRUE, test.use = 'MAST')

#Based on cellassign vignette.
#https://irrationone.github.io/cellassign/articles/introduction-to-cellassign.html#basic-usage
print("imprimindo_test.sce")
print(test.sce)
print("final_impressao_test.sce")

#Matrix construction.
#Markers list.
#Based on cellassign vignette.
#https://irrationone.github.io/cellassign/articles/introduction-to-cellassign.html#constructing-a-marker-gene-matrix
marker_gene_list <- list(
  TCD4naive = c("human----CD3E", "human----CD4", "human----CD45RA", "human----CCR7"),
  TCD4Th1 = c("human----CD3E", "human----CD4","human----TBX21", "human----IFNG", "human----IL2", "human----IL12", "human----TNF"),
  TCD4Th2 = c("human----CD3E", "human----CD4", "human----GATA3", "human----IL4", "human----IL5", "human----IL13"),
  TCD4Th17 = c("human----CD3E", "human----CD4", "human----RORGT", "human----IL17A", "human----IL17F", "human----IL21"),
  TCD4TReg = c("human----CD3E", "human----CD4", "human----FOXP3", "human----TGPB", "human----IL10"),
  TCD8naive = c("human----CD3E", "human----CD4", "human----CD45RA", "human----CCR7"),
  TCD8Citotoxic = c("human----GNLY", "human----PRF1", "human----CD3E", "human----CD8")
)

print(marker_gene_list)


marcadores <- marker_list_to_mat(marker_gene_list)
print(marcadores)
print("Fim dos marcadores")

#identifies markers within the dataset.
#https://nbisweden.github.io/single-cell_sib_scilifelab/session-differential-expression/celltype_assignment.html
shared <- intersect(rownames(marcadores), rownames(test.sce))

#shared <- intersect(rownames(marker_gene_mat), rownames(sce))
print("Shared")
print(shared)
print("######Shared")

#Based on cellassign vignette.
#https://irrationone.github.io/cellassign/articles/introduction-to-cellassign.html#basic-usage
print("Print Size Factors")
test.sce <- scran::computeSumFactors(test.sce)
s1 <- sizeFactors(test.sce)

print(s1)
print("######Print Size Factors")

#cellassign fit
#Based on cellassign vignette.
#https://irrationone.github.io/cellassign/articles/introduction-to-cellassign.html#constructing-a-marker-gene-matrix
fit <- cellassign(exprs_obj = test.sce[shared,], 
                  marker_gene_info = marcadores[shared,], 
                  s = s1, 
                  learning_rate = 1e-2, 
                  shrinkage = TRUE,
                  verbose = TRUE)

print("FIT######")
print(fit)
print("FIT######")

print("CELLLTYPES#####")
celltypes(fit)
print("CELLPROBS######")
cellprobs(fit)

#Plot heatmap
#Baseado no link da vignette do cellassign.
#https://irrationone.github.io/cellassign/articles/introduction-to-cellassign.html#basic-usage
pheatmap::pheatmap(cellprobs(fit))

#https://bioconductor.org/packages/devel/bioc/vignettes/scater/inst/doc/overview.html
#https://bioconductor.riken.jp/packages/3.4/bioc/vignettes/scater/inst/doc/vignette.html#plots-of-expression-values
#https://www.singlecellcourse.org/single-cell-rna-seq-analysis-using-seurat.html#cell-type-annotation-using-singler
#adds a column with the cell types identified by Cellassign in the test.sce object.
colData(test.sce)$cellassign_type <- fit$cell_type

#print the number of identified cell types.
table(test.sce$cellassign_type)

#Plot 2D pelo UMAP.
dpi = 300
png(file="plot_umap.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
plotUMAP(test.sce, colour_by = 'cellassign_type', point_size = 2)

#Plot the expression of a given gene with cellassign_type coloring.
png(file="violin_tnf.png", width = dpi*30, height = dpi*24, units = "px",res = dpi,type='cairo')
plotExpression(test.sce, features = "human----TNF", x = "cellassign_type",colour_by = 'cellassign_type')

