# Endo_Immune
## Codes for "Exploring the pathogenesis of atherosclerosis from a perspective of endothelial-immune interaction mechanisms"

To investigate the crosstalk between endothelial and immune cells and its role in atherosclerosis, we obtained scRNA-seq data of vascular cells from the adventitial and medial/endothelial layers of wild-type and Apoe-/- adult mice [(GSE140811)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140811). Raw sequences was downloaded using sratools and aligned to the reference genome GRCm38 using Cellranger v.7.1.0. The analysis procedures were provided as [scRNA-seq analysis.ipynb](./scRNAseq_analysis.ipynb):  
(1) Preprocessing and quality control  
(2) Imputation(MAGIC)  
(3) Normalization and scaling  
(4) PCA and UMAP calculation  
(5) Clustering and annotation  
(6) Pathway analysis (GOBP)  
(7) [Cell-cell communication analysis](./Cellchat_final.ipynb)  
The codes used to generate the figures could be found in [scRNAseq_Plottings.ipynb](./scRNAseq_Plottings.ipynb). The codes for Seurat object conversion (from python to R) were in [Seurat object conversion.ipynb](./Seurat_Conversion.ipynb).

To further confirm the results, we performed RNA-seq analysis on the CD31+ arterial endothelial cells obtained from mice fed with 9 weeks of high-fat diet and 3 extra weeks of Rux/vehicle. Data analysis procedures included:  
(1) [Normalization and DE analysis](./RNAseq_edgeR.R)  
(2) [Visualization](./RNAseq_Plottings.R)  

We also performed serum cytokine assays on mice fed with 9 weeks of high-fat diet and 3 extra weeks of Rux/vehicle. The data was analysed using FCAP 3.0 software:  
(1) Design: All 10 standards included  
(2) Instrument Settings: SSC-H as scatter param, R660-A(FSC-A) as clustering param, Y585-A as reporter param.  
(3) Debris Filtering: Manual gating on SSC-A / FSC-H  
(4) Manual clustering: Not applied (All Automatic)  
(5) Standard settings: 0; 20-5000 (Automatic calculation)  
The codes for visualization could be found in [Cytokine_Boxplot.R](./Cytokine_Boxplot.R).
