#packages and create directories
library("metacell")
library("dplyr")
library("knitr")
if(!dir.exists("severe")) dir.create("severe/")
scdb_init("severe/", force_reinit=T)
if(!dir.exists("results_severe")) dir.create("results_severe/")
scfigs_init("results_severe/")
set.seed(27)

#import object of Seurat
#severe
mcell_import_scmat_tsv(mat_nm = "all", fn = "/scratch/inova-covd19/vanessa.silva/paper/features_not_infected_severe.tsv", dset_nm = "/scratch/inova-covd19/vanessa.silva/paper/metadata_severe.csv", force = F)
mat = scdb_mat("all")
print(dim(mat@mat))
head(mat@mat)


#Filters
nms = c(rownames(mat@mat), rownames(mat@ignore_gmat))
ig_genes = c(grep("ERCC", nms, v=T), 
                grep("^IGH",nms,v=T),
                grep("^IGK", nms, v=T), 
                grep("^IGL", nms, v=T),
                grep("^IGJ", nms, v=T))

bad_genes = unique(c(ig_genes))
print(bad_genes)
mcell_mat_ignore_genes(new_mat_id="filtered_matrix", mat_id="all", bad_genes, reverse=F) 

#Calculates gene set statistics 
mcell_add_gene_stat(gstat_id="test", mat_id="filtered_matrix", force=T)
mcell_gset_filter_varmean(gset_id="test_feats", gstat_id="test", T_vm=0.4, force_new=T)
mcell_gset_filter_cov(gset_id = "test_feats", gstat_id="test", T_tot=30, T_top3=4)
mcell_plot_gstats(gstat_id="test", gset_id="test_feats")

#Calculates Knn matrix, resamples and clustering
mcell_add_cgraph_from_mat_bknn(mat_id="filtered_matrix",gset_id = "test_feats" ,graph_id="new_graph",K=100,dsamp=T) 
mcell_coclust_from_graph_resamp(coc_id="new_coc500",graph_id="new_graph",min_mc_size=20,p_resamp=0.70, n_resamp=500)
mcell_mc_from_coclust_balanced(coc_id="new_coc500",mat_id="filtered_matrix",mc_id="new_mc",K=30, min_mc_size=30, alpha=2) 

#Plot outliers
mcell_mc_split_filt(new_mc_id="new_mc_f",  
            mc_id="new_mc",  
            mat_id="filtered_matrix",
            T_lfc=4, plot_mats=F)

#Hierarchy of the confusion matrix
mc_hc = mcell_mc_hclust_confu(mc_id="new_mc_f",graph_id="new_graph")
mc_sup = mcell_mc_hierarchy(mc_id="new_mc_f",mc_hc=mc_hc, T_gap=0.04)
print(mc_sup[[17]])
print(mc_sup[[3]])

mc_colorize_sup_hierarchy(mc_id="new_mc_f", supmc=mc_sup, supmc_key="/scratch/inova-covd19/vanessa.silva/paper/supmc_new_severe.csv", gene_key=NULL)
mcell_mc_plot_hierarchy(mc_id="new_mc_f",graph_id="new_graph",mc_order=mc_hc$order, sup_mc = mc_sup, width=3500, height=3500, min_nmc=2,show_mc_ids = T) 



#Plot heatmaps with metacells
mcell_mc_plot_marks(mc_id="new_mc_f", gset_id= "test_feats", mat_id="filtered_matrix",plot_cells = T) 

#2D plot of metacells 
mc2d_knn = mcell_mc2d_force_knn(mc2d_id="id_2dproj",mc_id="new_mc_f", graph_id="new_graph") 
tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
mcell_mc2d_plot(mc2d_id="id_2dproj")

