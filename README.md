# CellHeap


This repository contains supplementary materials for the paper "CellHeap: a Workflow for Optimizing COVID-19 Single-Cell RNA-seq Data Processing in the Santos Dumont Supercomputer".

CellHeap is a configurable, portable, and robust workflow for scRNA-seq customizable analyses, with quality control throughout the execution steps, ensuring reliable results and deployable on supercomputers.

CellHeap comprises six phases: (1) sample curation, (2) gene count matrix generation, (3) quality control, (4) dimensionality reduction and clustering analysis, (5) cell-level analysis, and (6) gene-level analysis.

Selection criteria for sample curation are the following:

(i) to link to supplementary files on the measurement of genes of experiments;
(ii) to access the sample's raw data through SRA selector links;
(iii) to verify all series samples belong to single species;
(iv) to verify samples come from non-bacterial species;
(v) to verify the source of the description of experiments, protocols, cell line information, cell type, and disease as listed by Experimental Factor Ontology (EFO) or publication data;
(vi) to check if metadata matches names of samples;
(vii) to verify scRNA-Seq experiments use protocols of Smart-seq2, Smart-like, Drop-seq, Seq-well, 10xV2 (3 and 5 prime), or 10xV3 (3 prime).

This repository contains a deployment of CellHeap for a dataset made available by Liao et al. 2020:

Liao, M., Liu, Y., Yuan, J., Wen, Y., Xu, G., Zhao, J., Cheng, L., Li, J., Wang, X., Wang, F., et al.: Single-cell landscape of bronchoalveolar immune cells in patients with covid-19. Nature medicine 26(6), 842–844 (2020)

In that paper, the authors considered 13 patients, of which four patients were controls, 3 presented mild symptoms, and 6 six developed severe symptoms of COVID-19. However, when we applied the dataset selection criteria related to phase 1, one dataset was discarded relative to one of the controls.

The source code of CellHeap for Liao et al. 2020's dataset is in the code directory. <i> Cellranger count </i> output reports are in the Data/Cellranger directory. <i> Metacell </i> 2D projection results are in the Data/Metacell directory.
