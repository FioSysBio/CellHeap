#packages and create directories
library("metacell")
library("dplyr")
library("knitr")
if(!dir.exists("mild")) dir.create("mild/")
scdb_init("mild/", force_reinit=T)
if(!dir.exists("results_mild")) dir.create("results_mild/")
scfigs_init("results_mild/")
id = "lung_new"
ord_id = "lung_new_sorted"


#import object of Seurat
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
                grep("^IGJ", nms, v=T))#############
bad_genes = unique(c(ig_genes)) ##############
print(bad_genes)
mcell_mat_ignore_genes(new_mat_id="filtered_matrix", mat_id="all", bad_genes, reverse=F) ############# mudei nome de new_mat_id

#Calculates gene dataset statistics
#set.seed(27)
#mcell_add_gene_stat(gstat_id=id, mat_id=id, force=T)
#gstat = scdb_gstat(id)
#write.csv(gstat, file="mild/gstat.csv")
#t_vm = 0.4
#mcell_gset_filter_multi(gstat_id=id, gset_id=id, T_tot=50, T_top3=4, T_vm = t_vm, force_new = T)
#save(file="mild/bad_marks.Rda",bad_marks)

gset = gset_new_gset(vm_set, sprintf("VM %f",T_vm))



scdb_add_gset("test_feats_filtered", gset)
mcell_gset_add_gene(gset_id="test_feats_filtered", genes="CD3,CD4,CD45RA,CCR7high,CD3,CD4,TBX21,IFNG,IL2,IL12,TNF,CD3,CD4,GATA3,IL4,IL5,IL13,CD3,CD4,RORGT,IL17A,IL17F,IL21,CD3,CD4,FOXP3,TGPB,IL10,CD3,CD4,CD45RA,CCR7,GNLY,PRF1,CD3,CD8,CD8,FCGR3B,PI3,G0S2,CPA3,MS4A2,TPSAB1,TPSAB2,STAT1,TNF,IL6,IL1B,CXCL10,CXCL9,IDO1,IRF5,MARCO,TGFBR2,NKG2D,TCRG,CD79A,CD79B,MS4A1,SCGB1A1,SCGB3A1,MSMB,KRT5,AQP3,TP63,CAPS,TPPP3,RSPH1,KRT13,KRT4,SPRR3,KRT8,KRT18,MMP7,SFTPC,SFTPA1,SFTPB", subset_id = 1) ####function to add genes to a gene set

mcell_add_gene_stat(gstat_id="genes", mat_id="filtered_matrix", force=T)

#Calculates Knn matrix, resamples and cluster
set.seed(27)
mcell_add_cgraph_from_mat_bknn(mat_id="filtered_matrix",gset_id = "test_feats_filtered" ,graph_id="new_graph",K=100,dsamp=T)###mudei mat_id, graph id

cgraph = scdb_cgraph(id)
kable(head(cgraph@edges))

mcell_coclust_from_graph_resamp(coc_id="new_coc500",graph_id="new_graph",min_mc_size=20,p_resamp=0.70, n_resamp=500)#### mudei coc_id, graph_id

coclust = scdb_coclust(id)
kable(head(coclust@coclust))

mcell_mc_from_coclust_balanced(coc_id="new_coc500",mat_id="filtered_matrix",mc_id="new_mc",K=100, min_mc_size=30, alpha=2) ####mudei coc_id, mat_id, mc_id

write.csv(coclust@coclust, file="mild/coclust.csv")
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
mcell_mc_plot_marks(mc_id="new_mc_f", gset_id= "test_feats_filtered", mat_id="filtered_matrix",plot_cells = T) ####mudei mc_id, gset_id, mat_id

#Plot 2D metacells
mc2d_knn = mcell_mc2d_force_knn(mc2d_id="id_2dproj",mc_id="new_mc_f", graph_id="new_graph") ####### mudei mc_id, graph_id
tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
mcell_mc2d_plot(mc2d_id="id_2dproj")

#Calculate and Plot the confusion matrix
#set.seed(27)
mc_hcxxx = mcell_mc_hclust_confu(mc_id="new_mc_f",graph_id="new_graph") #####mudei mc_id, graph_id
#set.seed(27)
mc_sup = mcell_mc_hierarchy(mc_id="new_mc_f",mc_hc=mc_hcxxx, T_gap=0.04)####### mudei mc_id
save(file="severe/mc_hc_sup.Rda",mc_hc,mc_sup)
mcell_mc_plot_hierarchy(mc_id="new_mc_f",graph_id="new_graph",mc_order=mc_hcxxx$order,sup_mc = mc_sup,width=3500, height=3500, min_nmc=2,show_mc_ids = T) ###
########  acima mudei mc_id, graph_id

mcell_mc_plot_confusion( mc_id="new_mc_f",  graph_id="new_graph") ######mudei mc_id, graph_id


#Identification of cell subpopulations###############
mcell_mc_export_tab(mc_id ="new_mc_f", gstat_id = "genes", mat_id = "filtered matrix", T_fold=2, metadata_fields=NULL)##################alterei o parâmetro mc_id, gstat_id e mat_id, baseado nas linhas 31 e 56.

#Gene enrichment plot in metacells
lfp <- log2 (mc@mc_fp) ##
head (lfp, n=25L)##
write.csv(lfp, file="severe/lfp.csv", row.names=TRUE)###
png("results_mild/barplot1.png",h=1000,w=1000);barplot(lfp["IFNG",],col=marks_colors,las=2,main="IFNG",cex.main=3,cex.axis=1,ylab="log2FC",xlab="metacells");dev.off()

