#packages and create directories
library("metacell")
library("dplyr")
if(!dir.exists("severe")) dir.create("severe/")
scdb_init("severe/", force_reinit=T)
if(!dir.exists("results_severe")) dir.create("results_severe/")
scfigs_init("results_severe/")
id = "lung_new"
ord_id = "lung_new_sorted"
library(knitr)

#import object of Seurat
#severe
mcell_import_scmat_tsv(mat_nm = "all", fn = "/scratch/inova-covd19/vanessa.silva/features_not_infected_severe.tsv", dset_nm = "/scratch/inova-covd19/vanessa.silva/metadata_severe.csv", force = F)
mat = scdb_mat("all")
print(dim(mat@mat))

#Filter mitochondrial genes
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
ig_genes = c(grep("ERCC", nms, v=T), 
                grep("^IGH",nms,v=T),
                grep("^IGK", nms, v=T), 
                grep("^IGL", nms, v=T),
                grep("^IGJ", nms, v=T)) ###############
bad_genes = unique(c(ig_genes)) ##############
print(bad_genes)
mcell_mat_ignore_genes(new_mat_id="filtered_matrix", mat_id="all", bad_genes, reverse=F) ############# mudei nome de new_mat_id
############mcell_mat_ignore_small_cells(id, id, 300)

#Calculates gene dataset statistics
#set.seed(27)
#mcell_add_gene_stat(gstat_id=id, mat_id=id, force=T)
#gstat = scdb_gstat(id)
#write.csv(gstat, file="severe/gstat.csv")
#t_vm = 0.4
#mcell_gset_filter_multi(gstat_id=id, gset_id=id, T_tot=50, T_top3=4, T_vm = t_vm, force_new = T)
#save(file="severe/bad_marks.Rda",bad_marks)


#gset_add_genes(gset, genes, subset_id)

###leitura de arquivo


mcell_gset_add_gene(gset_id="test_feats_filtered", genes="CD3,CD4,CD45RA,CCR7high", subset_id = 1) ####function to add genes to a gene set


###################. ALTERNATIVA AO COMANDO ACIMA #################################################
###GENE_SET_FILE = read.delim("NOME DO ARQUIVO COM GENES MARCADORES", sep="\t", stringsAsFactors=F)
####mcell_gset_add_gene(gset_id="test_feats_filtered", genes=GENE_SET_FILE, subset_id = 1) ####function to add genes to a gene set
######################### FIM DA ALTERNATIVA #######################################################


#marker_gset = scdb_gset("test_feats")
#imarker_gset = gset_new_restrict_gset(gset = marker_gset, filt_gset = lateral_gset, inverse = T, desc = "cgraph gene markers w/o lateral genes")
#scdb_add_gset("test_feats_filtered", marker_gset)


#Calculates Knn matrix, resamples and cluster
set.seed(27)
mcell_add_cgraph_from_mat_bknn(mat_id="filtered_matrix",gset_id = "test_feats_filtered" ,graph_id="new_graph",K=100,dsamp=T) ###mudei mat_id, graph id

cgraph = scdb_cgraph(id)
kable(head(cgraph@edges))
#set.seed(27) não devemos resetar o gerador de números pseudoaleatórios a toda hora
mcell_coclust_from_graph_resamp(coc_id="new_coc500",graph_id="new_graph",min_mc_size=20,p_resamp=0.70, n_resamp=500)#### mudei coc_id, graph_id

coclust = scdb_coclust(id)
kable(head(coclust@coclust))
#set.seed(27)
mcell_mc_from_coclust_balanced(coc_id="new_coc500",mat_id="filtered_matrix",mc_id="new_mc",K=100, min_mc_size=30, alpha=2) ####mudei coc_id, mat_id, mc_id

write.csv(coclust@coclust, file="severe/coclust.csv")
mc<- scdb_mc(id)
scdb_add_mc(id,mc)

#Plot outliers
mcell_mc_split_filt(new_mc_id="new_mc_f",  
            mc_id="new_mc",  
            mat_id="filtered_matrix",
            T_lfc=3, plot_mats=F) ####mudei new_mc_id, mc_id, mat_id
			
#Selecting markers genes
marks_colors = read.delim("/scratch/inova-covd19/vanessa.silva/mc_colorize.txt", sep="\t", stringsAsFactors=F)
kable(head(marks_colors))
mc_colorize(new_mc_id = "new_mc_f", marker_colors=marks_colors,override=T)#####mudei new_mc_id, setei mc_id to default (new_mc_id)

#Heatmap plot with metacells
mcell_gset_from_mc_markers(gset_id=paste0(id,"_markers"), mc_id=id)
mcell_mc_plot_marks(mc_id=id, gset_id=paste0(id,"_markers"), mat_id=id,plot_cells = T)

#Plot 2D metacells
mc2d_knn = mcell_mc2d_force_knn(mc2d_id="id_2dproj",mc_id=id, graph_id=id)
tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
mcell_mc2d_plot(mc2d_id="id_2dproj")

#Calculate and plot confusion Matrix
set.seed(27)
mc_hc = mcell_mc_hclust_confu(mc_id=id,graph_id=id)
set.seed(27)
mc_sup = mcell_mc_hierarchy(mc_id=id,mc_hc=mc_hc, T_gap=0.04)
save(file="severe/mc_hc_sup.Rda",mc_hc,mc_sup)
mcell_mc_plot_hierarchy(mc_id=id,graph_id=id,mc_order=mc_hc$order,sup_mc = mc_sup,width=3500, height=3500, min_nmc=2,show_mc_ids = T)
mcell_mc_plot_confusion( mc_id=id,  graph_id=id)

lfp <- log2 (mc@mc_fp)
head (lfp, n=25L)
write.csv(lfp, file="severe/lfp.csv", row.names=TRUE)

#Identification of cell subpopulations
mcell_mc_export_tab(mc_id = id, gstat_id = id, mat_id = "all", T_fold=2, metadata_fields=NULL)
lfp <- log2(mc@mc_fp)

#Gene enrichment plot in metacells
png("results_severe/barplot1.png",h=1000,w=1000);barplot(lfp["IFNG",],col=marks_colors,las=2,main="IFNG",cex.main=3,cex.axis=1,ylab="log2FC",xlab="metacells");dev.off()

