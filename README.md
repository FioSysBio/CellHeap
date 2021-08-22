# CellHeap

CellHeap is a configurable, portable, and robust workflow for scRNA-seq customizable analyses, with quality control throughout the execution steps, ensuring reliable results and deployable on supercomputers. 

CellHeap is composed of six phases: sample curation, gene count matrix generation, quality control, dimensionality reduction and clustering analysis, cell level analysis, and gene level analysis.

Selection criteria for sample curation are the following:

(i) to link to supplementary files on the measurement of genes of experiments; 

(ii) to access to the sample's raw data through SRA selector links; 

(iii) to verify all series samples belong to single species; 

(iv) to verify samples come from non-bacterial species; 

(v) to verify the source of description of experiments, protocols, cell line information, cell type, and disease as listed by Experimental Factor Ontology (EFO) or publication data; 

(vi) to check if metadata matches names of samples; 

(vii) to verify scRNA-Seq experiments use protocols of Smart-seq2, Smart-like, Drop-seq, Seq-well, 10xV2 (3 and 5 prime), or 10xV3 (3 prime). 

This repository contains a deployment of CellHeap for a dataset made available by Liao et al. 2020:

Liao, M., Liu, Y., Yuan, J., Wen, Y., Xu, G., Zhao, J., Cheng, L., Li, J., Wang, X.,Wang, F., et al.: Single-cell landscape of bronchoalveolar immune cells in patientswith covid-19. Nature medicine26(6), 842–844 (2020)

In that paper, the authors considered 13 patients, from which four patients were controls, 3 presented mild symptoms, and 6 six patients developed severe symptoms of COVID-19. However, when we applied the dataset selection criteria related to phase 1, one dataset was discarded relative to one of the controls.

Source code of CellHeap for Liao et al. 2020's dataset is in the code directory. Cellranger count outpput reports are in the Data/Cellranger directory.
