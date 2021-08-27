##### verificar https://tanaylab.github.io/metacell/articles/a-basic_pbmc8k.html


#packages and create directories
library("metacell")
library("dplyr")
library("knitr")
if(!dir.exists("severe")) dir.create("severe/")
scdb_init("severe/", force_reinit=T)
if(!dir.exists("results_severe")) dir.create("results_severe/")
scfigs_init("results_severe/")
id = "lung_new"
ord_id = "lung_new_sorted"


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
#gstat = scdb_gstat(id)
#write.csv(gstat, file="severe/gstat.csv")
#t_vm = 0.4
#mcell_gset_filter_multi(gstat_id=id, gset_id=id, T_tot=50, T_top3=4, T_vm = t_vm, force_new = T)
#save(file="severe/bad_marks.Rda",bad_marks)


#gset_add_genes(gset, genes, subset_id)

###leitura de arquivo

#Em que consiste a variável vm_set?
#o valor T_vm = 0.4, segundo o paper de Viral-track.

############### RESPOSTA #############################################################
### Vamos lá
####
#####O código da chamada está a seguir
#
#' Generating a new gset
#
#' @param sets a vector where names are genes and values are set ids
#' @param desc texual description of the gset
#' @export
#gset_new_gset = function(sets, desc)
#{
#	return(tgGeneSets(sets, desc))
#}
#
#. Para criar um named vetor de set ids que eu vou chamar de vm_set
# vm_set = c( 1 , 1 ).  ######set id - como tenho um único set o id é o mesmo
# names( vm_set ) = c( "NOME_DO_GENE1" , "NOME_DO_GENE2" )
#. 
#
#. O segundo parâmetro é uma descrição de texto qualquer - ex: "Genes da Andrea"
#
# A chamada poderia ficar assim
#
# 
# 
gset = gset_new_gset(vm_set, sprintf("VM %f",T_vm))
vec= c("human----CD3","human----CD4","human----CD45RA","human----CCR7high","human----CD3","human----CD4","human----TBX21","human----IFNG","human----IL2","human----IL12","human----TNF","human----CD3","human----CD4","human----GATA3","human----IL4","human----IL5","human----IL13","human----CD3","human----CD4","human----RORGT","human----IL17A","human----IL17F","human----IL21","human----CD3","human----CD4","human----FOXP3","human----TGPB","human----IL10","human----CD3","human----CD4","human----CD45RA","human----CCR7","human----GNLY","human----PRF1","human----CD3","human----CD8","human----CD8","human----FCGR3B","human----PI3","human----G0S2","human----CPA3","human----MS4A2","human----TPSAB1","human----TPSAB2","human----STAT1","human----TNF","human----IL6","human----IL1B","human----CXCL10","human----CXCL9","human----IDO1","human----IRF5","human----MARCO","human----TGFBR2","human----NKG2D","human----TCRG","human----CD79A","human----CD79B","human----MS4A1","human----SCGB1A1","human----SCGB3A1","human----MSMB","human----KRT5","human----AQP3","human----TP63","human----CAPS","human----TPPP3","human----RSPH1","human----KRT13","human----KRT4","human----SPRR3","human----KRT8","human----KRT18","human----MMP7","human----SFTPC","human----SFTPA1","human----SFTPB") 
gset = gset_new_gset(sets = vec, desc = "genesAndrea")

vec = c("CD3","CD4","CD45RA","CCR7high","CD3","CD4","TBX21","IFNG","IL2","IL12","TNF","CD3","CD4","GATA3","IL4","IL5","IL13","CD3","CD4","RORGT","IL17A","IL17F","IL21","CD3","CD4","FOXP3","TGPB","IL10","CD3","CD4","CD45RA","CCR7","GNLY","PRF1","CD3","CD8","CD8","FCGR3B","PI3","G0S2","CPA3","MS4A2","TPSAB1","TPSAB2","STAT1","TNF","IL6","IL1B","CXCL10","CXCL9","IDO1","IRF5","MARCO","TGFBR2","NKG2D","TCRG","CD79A","CD79B","MS4A1","SCGB1A1","SCGB3A1","MSMB","KRT5","AQP3","TP63","CAPS","TPPP3","RSPH1","KRT13","KRT4","SPRR3","KRT8","KRT18","MMP7","SFTPC","SFTPA1","SFTPB")

gset = gset_new_gset(sets = vec, desc = "genesAndrea")


scdb_add_gset("test_feats_filtered", gset)
mcell_gset_add_gene(gset_id="test_feats_filtered", genes= vec, subset_id = 1) ####function to add genes to a gene set   ###alterei genes = vec


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

cgraph = scdb_cgraph(id="new_graph") ###alterei o id="new_graph"
kable(head(cgraph@edges))
#set.seed(27) não devemos resetar o gerador de números pseudoaleatórios a toda hora
mcell_coclust_from_graph_resamp(coc_id="new_coc500",graph_id="new_graph",min_mc_size=20,p_resamp=0.70, n_resamp=500)#### mudei coc_id, graph_id

coclust = scdb_coclust(id="new_coc500") ###alterei o id="new_coc500"
kable(head(coclust@coclust))
#set.seed(27)
mcell_mc_from_coclust_balanced(coc_id="new_coc500",mat_id="filtered_matrix",mc_id="new_mc",K=100, min_mc_size=30, alpha=2) ####mudei coc_id, mat_id, mc_id

write.csv(coclust@coclust, file="severe/coclust.csv")
mc <- scdb_mc(id="new_mc") ### alterei para id="new_mc"
scdb_add_mc(id="new_mc",mc) ### alterei para id="new_mc"

#Plot outliers
mcell_mc_split_filt(new_mc_id="new_mc_f",  
            mc_id="new_mc",  
            mat_id="filtered_matrix",
            T_lfc=3, plot_mats=F) ####mudei new_mc_id, mc_id, mat_id   ##verificar parâmetro T_lfc.
			
#colorize metacell from selected marker genes
marks_colors = read.delim("/scratch/inova-covd19/vanessa.silva/paper/features_genes.txt", sep="\t", stringsAsFactors=F) ### alterei o caminho
kable(head(marks_colors))
mc_colorize(new_mc_id = "new_mc_f", mc_id = "new_mc", marker_colors=marks_colors,override=T)#####mudei new_mc_id, setei mc_id to default (new_mc_id)
#mc_colorize(new_mc_id, mc_id = new_mc_id, marker_colors = NULL, override = T)
# mc_colorize("new_mc_f", marker_colors=marks_colors)
# mc_colorize_default(mc_id="new_mc", spectrum = NULL)


#Heatmap plot with metacells
#mcell_gset_from_mc_markers(gset_id="test_feats_filtered", mc_id="new_mc") #### alterei o "new_mc"  ###
mcell_mc_plot_marks(mc_id="new_mc", gset_id= "test_feats_filtered", mat_id="filtered_matrix",plot_cells = T) ####mudei mc_id, gset_id, mat_id

#Plot 2D metacells
mc2d_knn = mcell_mc2d_force_knn(mc2d_id="id_2dproj",mc_id="new_mc", graph_id="new_graph") ####### mudei mc_id, graph_id
tgconfig::set_param("mcell_mc2d_height",1000, "metacell")
tgconfig::set_param("mcell_mc2d_width",1000, "metacell")
mcell_mc2d_plot(mc2d_id="id_2dproj")

#Calculate and plot confusion Matrix
#set.seed(27)
mc_hcxxx = mcell_mc_hclust_confu(mc_id="new_mc_f",graph_id="new_graph") #####mudei mc_id, graph_id
#set.seed(27)
mc_sup = mcell_mc_hierarchy(mc_id="new_mc_f",mc_hc=mc_hcxxx, T_gap=0.04)####### mudei mc_id
mcell_mc_plot_hierarchy(mc_id="new_mc_f",graph_id="new_graph",mc_order=mc_hcxxx$order,sup_mc = mc_sup,width=3500, height=3500, min_nmc=2,show_mc_ids = T) ###
########  acima mudei mc_id, graph_id

mcell_mc_plot_confusion( mc_id="new_mc_f",  graph_id="new_graph") ######mudei mc_id, graph_id

#comando necessário para exportar tabela da linha 118.
#A seleção irá considerar somente os genes definidos pela Andrea ou todos?

####### RESPOSTA: se você der como entrada filtered_matrix, as estatísticas serão calculadas em cima de todos os genes
#### para calcular estatísticas limitadas ao genes da Andrea, terá de ser gerada nova matriz contendo apenas genes da Andrea

#A seleção irá considerar somente os genes definidos pela Andrea ou todos? 
mcell_add_gene_stat(gstat_id="genes", mat_id="filtered_matrix", force=T)

#Identification of cell subpopulations###############
mcell_mc_export_tab(mc_id ="new_mc_f", gstat_id = "genes", mat_id = "filtered matrix", T_fold=2, metadata_fields=NULL)##################alterei o parâmetro mc_id, gstat_id e mat_id, baseado nas linhas 31 e 56.

#Gene enrichment plot in metacells
lfp <- log2 (mc@mc_fp) ##
head (lfp, n=25L)##
write.csv(lfp, file="severe/lfp.csv", row.names=TRUE)###
save(file="severe/mc_hc_sup.Rda",mc_hc,mc_sup,lfp)####
png("results_severe/barplot1.png",h=1000,w=1000);barplot(lfp["IFNG",],col=marks_colors,las=2,main="IFNG",cex.main=3,cex.axis=1,ylab="log2FC",xlab="metacells");dev.off()

