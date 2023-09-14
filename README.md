# Overview
This is a benchmarking project for evaluation of deconvolution methods in decoding neuronal composition in GBM context, understanding interaction between neural-tumor cells. 
File naming rules as follow:
* data: references data
* scripts: function or modual code
* jupyter: demo or tutorial

# Datasets
## For signatures and simulation generation

|Short names|Long names|Description|References and links|
|:-:|:-:|---|---|
|YU|GSE117891_Yu_2020|collected 73 surgery points within more than 6000 cells in 10X without prefiltering, but focus on immunological cells. Therefore, we just using count data but do not consider annotation in this study|Yu K, Hu Y, Wu F, Guo Q, Qian Z, Hu W, Chen J, Wang K, Fan X, Wu X, Rasko JE, Fan X, Iavarone A, Jiang T, Tang F, Su XD. Surveying brain tumor heterogeneity by single-cell RNA-sequencing of multi-sector biopsies. Natl Sci Rev. 2020 Aug;7(8):1306-1318. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117891)|
|BH|PRJNA579593_Bhaduri_2020|10X and C1 single cell data, have ensured contain intact cell type(cannot filter before sequecing).This dataset have used to study GBM in neural development perspect as detailed information of subtype cells, thus it will be used in main reference and annotation standard|Bhaduri A, Di Lullo E, Jung D, Müller S, Crouch EE, Espinosa CS, Ozawa T, Alvarado B, Spatazza J, Cadwell CR, Wilkins G, Velmeshev D, Liu SJ, Malatesta M, Andrews MG, Mostajo-Radji MA, Huang EJ, Nowakowski TJ, Lim DA, Diaz A, Raleigh DR, Kriegstein AR. Outer Radial Glia-like Cancer Stem Cells Contribute to Heterogeneity of Glioblastoma. Cell Stem Cell. 2020 Jan 2;26(1):48-63.e6. [SRA](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA579593)|
|LW|GSE138794_Lin Wang_2019|10X GBM sc data, focus on maglinant cells, so do not consider annotation in this study|Wang L, Babikir H, Müller S, Yagnik G, Shamardani K, Catalan F, Kohanbash G, Alvarado B, Di Lullo E, Kriegstein A, Shah S, Wadhwa H, Chang SM, Phillips JJ, Aghi MK, Diaz AA. The Phenotypes of Proliferating Glioblastoma Cells Reside on a Single Axis of Variation. Cancer Discov. 2019 Dec;9(12):1708-1719.[GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi)|

# Deconvolution algorithms
|Name|Classification|References and installation|
|:-:|---|---|
|BayesPrism|multinomial Probabilistic model and updata references based on background signatures to mitigate inflence of maglinant cells|Chu, T., Wang, Z., Pe’er, D. et al. Cell type and gene expression deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and single-cell RNA sequencing in oncology. Nat Cancer 3, 505–517 (2022).[GitHub](https://github.com/Danko-Lab/BayesPrism/tree/main)|
|Bulk2Space|NNLS+$\beta$-VAE deconvolution in a single cell resolution|Liao, J., Qian, J., Fang, Y. et al. De novo analysis of bulk RNA-seq data at spatially resolved single-cell resolution. Nat Commun 13, 6498 (2022).[GitHub](https://github.com/ZJUFanLab/bulk2space)|
|Kassandra|Utility of decision tree counstructed a heirachical framework|Zaitsev A, Chelushkin M, Dyikanov D, Cheremushkin I, Shpak B, Nomie K, Zyrin V, Nuzhdina E, Lozinsky Y, Zotova A, Degryse S, Kotlov N, Baisangurov A, Shatsky V, Afenteva D, Kuznetsov A, Paul SR, Davies DL, Reeves PM, Lanuti M, Goldberg MF, Tazearslan C, Chasse M, Wang I, Abdou M, Aslanian SM, Andrewes S, Hsieh JJ, Ramachandran A, Lyu Y, Galkin I, Svekolkin V, Cerchietti L, Poznansky MC, Ataullakhanov R, Fowler N, Bagaev A. Precise reconstruction of the TME using bulk RNA-seq and a machine learning algorithm trained on artificial transcriptomes. Cancer Cell. 2022 Aug 8;40(8):879-894.e16.[Github](https://github.com/BostonGene/Kassandra/tree/main)|
|BLADE|log-normal Probabilistic model adopting large datasets|Andrade Barbosa, B., van Asten, S.D., Oh, J.W. et al. Bayesian log-normal deconvolution for enhanced in silico microdissection of bulk gene expression data. Nat Commun 12, 6106 (2021). [GitHub](https://github.com/tgac-vumc/BLADE/tree/master)|
|MuSiC/MuSiC2|Weighted non-negative least squares cross-validat references from different samples to reduce noise result in heterogenesis|Wang, X., Park, J., Susztak, K. et al. Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. Nat Commun 10, 380 (2019).[Github](https://github.com/xuranw/MuSiC)/Fan J, Lyu Y, Zhang Q, Wang X, Li M, Xiao R. MuSiC2: cell-type deconvolution for multi-condition bulk RNA-seq data. Brief Bioinform. 2022 Nov 19;23(6):bbac430.[Tutorial](https://xuranw.github.io/MuSiC/articles/pages/MuSiC2.html)|
|CIBERSORTx|nu-SVR (linear) enhance robustness|Newman, A.M., Steen, C.B., Liu, C.L. et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nat Biotechnol 37, 773–782 (2019).[CIBERSORT website](https://cibersortx.stanford.edu/)|
|DECODER(R version)|Unsupervised method using NMF + Regularized linear regression predict fraction of subtype|Peng XL, Moffitt RA, Torphy RJ, Volmar KE, Yeh JJ. De novo compartment deconvolution and weight estimation of tumor samples using DECODER. Nat Commun. 2019 Oct 18;10(1):4729.[Github](https://github.com/laurapeng/decoderr)|
|CDSeqR|a unsupervised model based on multinomial Probabilistic model|Kang K, Huang C, Li Y, Umbach DM, Li L. CDSeqR: fast complete deconvolution for gene expression data from bulk tissues. BMC Bioinformatics. 2021 May 24;22(1):262. doi: 10.1186/s12859-021-04186-5.[GitHub](https://github.com/kkang7/CDSeq_R_Package)|
|ISOpureR|multinomial Probabilistic modelusing k-means clustering, the genes are associated to expression profile across subtypes followed by minimization of the residual squared error|Anghel, C.V., Quon, G., Haider, S. et al. ISOpureR: an R implementation of a computational purification algorithm of mixed tumour profiles. BMC Bioinformatics 16, 156 (2015).[CRAN](https://cran.r-project.org/web/packages/ISOpureR/index.html)|

# Other packages

# Directory structure

# Other references
