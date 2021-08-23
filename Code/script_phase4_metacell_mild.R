#packages and create directories
library("metacell")
library("dplyr")
if(!dir.exists("mild")) dir.create("mild/")
scdb_init("mild/", force_reinit=T)
if(!dir.exists("results_mild")) dir.create("results_mild/")
scfigs_init("results_mild/")
id = "lung_new"
ord_id = "lung_new_sorted"
library(knitr)

#import objeto do Seurat
#mild
mcell_import_scmat_tsv(mat_nm = "all", fn = "/scratch/inova-covd19/vanessa.silva/features_not_infected_mild.tsv", dset_nm = "/scratch/inova-covd19/vanessa.silva/metadata_mild.csv", force = F)
mat = scdb_mat("all")
print(dim(mat@mat))

#Filter mitochondrial genes
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
ig_genes = c(grep("ERCC", nms, v=T), 
                grep("^IGH",nms,v=T),
                grep("^IGK", nms, v=T), 
                grep("^IGL", nms, v=T),
                grep("^IGL", nms, v=T))
bad_genes = unique(c(grep("^MT-", nms, v=T), grep("^MTMR", nms, v=T), grep("^MTND", nms, v=T), ig_genes))
print(bad_genes)
mcell_mat_ignore_genes(new_mat_id=id, mat_id="all", bad_genes, reverse=F) 
mcell_mat_ignore_small_cells(id, id, 300)

#Calculates gene dataset statistics
set.seed(27)
mcell_add_gene_stat(gstat_id=id, mat_id=id, force=T)
gstat = scdb_gstat(id)
write.csv(gstat, file="mild/gstat.csv")
t_vm = 0.4
mcell_gset_filter_multi(gstat_id=id, gset_id=id, T_tot=50, T_top3=4, T_vm = t_vm, force_new = T)
save(file="mild/bad_marks.Rda",bad_marks)

#Calculates Knn matrix, resamples and cluster
set.seed(27)
mcell_add_cgraph_from_mat_bknn(mat_id=id,gset_id = id,graph_id=id,K=100,dsamp=T)
cgraph = scdb_cgraph(id)
kable(head(cgraph@edges))
set.seed(27)
mcell_coclust_from_graph_resamp(coc_id=id,graph_id=id,min_mc_size=20,p_resamp=0.70, n_resamp=500)
coclust = scdb_coclust(id)
kable(head(coclust@coclust))
set.seed(27)
mcell_mc_from_coclust_balanced(coc_id=id,mat_id= id,mc_id= id,K=100, min_mc_size=30, alpha=2)
write.csv(coclust@coclust, file="mild/coclust.csv")
mc<- scdb_mc(id)
scdb_add_mc(id,mc)

#Plot outliers
mcell_mc_split_filt(new_mc_id=id,
            mc_id=id,
            mat_id=id,
            T_lfc=3, plot_mats=F)
			
#Selecting markers genes
marks_colors = read.delim("/scratch/inova-covd19/vanessa.silva/mc_colorize.txt", sep="\t", stringsAsFactors=F)
kable(head(marks_colors))
mc_colorize(new_mc_id = id, mc_id = id, marker_colors=marks_colors,override=T)

#Heatmap plot with metacells
mcell_gset_from_mc_markers(gset_id=paste0(id,"_markers"), mc_id=id)
mcell_mc_plot_marks(mc_id=id, gset_id=paste0(id,"_markers"), mat_id=id,plot_cells = T)

#Plot 2D metacells
mc2d_knn = mcell_mc2d_force_knn(mc2d_id="id_2dproj",mc_id=id, graph_id=id)
tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
mcell_mc2d_plot(mc2d_id="id_2dproj")

#Calculate and Plot the confusion matrix
set.seed(27)
mc_hc = mcell_mc_hclust_confu(mc_id=id,graph_id=id)
set.seed(27)
mc_sup = mcell_mc_hierarchy(mc_id=id,mc_hc=mc_hc, T_gap=0.04)
save(file="mild/mc_hc_sup.Rda",mc_hc,mc_sup)
mcell_mc_plot_hierarchy(mc_id=id,graph_id=id,mc_order=mc_hc$order,sup_mc = mc_sup,width=3500, height=3500, min_nmc=2,show_mc_ids = T)
mcell_mc_plot_confusion( mc_id=id,  graph_id=id)

lfp <- log2 (mc@mc_fp)
head (lfp, n=25L)
write.csv(lfp, file="mild/lfp.csv", row.names=TRUE)


#Identification of cell subpopulations
mcell_mc_export_tab(mc_id = id, gstat_id = id, mat_id = "all", T_fold=2, metadata_fields=NULL)
lfp <- log2(mc@mc_fp)

#Gene enrichment plot in metacells
png("results_mild/barplot1.png",h=1000,w=1000);barplot(lfp["IFNG",],col=marks_colors,las=2,main="IFNG",cex.main=3,cex.axis=1,ylab="log2FC",xlab="metacells");dev.off()
png("results_mild/barplot2.png",h=1000,w=1000);barplot(lfp["NFKB1",],col=marks_colors,las=2,main="NFKB1",cex.main=3,cex.axis=1,ylab="log2FC",xlab="metacells");dev.off()
