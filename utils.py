#Date: 20231002
#Version: V1.0.0
#Author: ZHOU YING



def preprocessing(adata,Sample_key,save_path,batch_key=None,batch_name=None):
    '''running all datasets by this process to filter low-quality cells and doublets
    refer to https://www.sc-best-practices.org/preprocessing_visualiza
    quality control performed on following three covariates:
    1. The number of counts per barcode (count depth)
    2. The number of genes per barconde
    3. The franction of counts from mitochondrial genes per barcode
    Using MAD (median absolute deviations) as threshold to filter all metrics
    MAD=median(|Xi-median(X)|), Xi is the respective QC metric of an observation
    cells are maked as outliers out of 5 MADs [Germain et al., 2020](https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html)'''
    import numpy as np
    import scanpy as sc
    import seaborn as sns
    from scipy.stats import median_abs_deviation
    import anndata as ad
    import pandas as pd
    import os
    from scipy.sparse import csr_matrix
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
    sc.settings.verbosity = 0
    sc.settings.set_figure_params(
        dpi=80,
        facecolor="white",
        frameon=False,
    )
    
    stats=pd.DataFrame(data=list(adata.obs[Sample_key].value_counts()),
             columns=["N.cells in RawData"],
             index=list(adata.obs[Sample_key].value_counts().index))
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes.
    adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))
    sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)

    if batch_name is not None:
        QC_plot_savepath=save_path+batch_name+"_qc_plot.pdf"
    elif batch_key is not None:
        QC_plot_savepath=save_path+adata.obs[batch_key].unique()[0]+"_qc_plot.pdf"
    else:
        print("Please give \"batch_name\" or \"batch_key\" ")
        
    with PdfPages(QC_plot_savepath) as pdf:
        p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
        pdf.savefig()
        # sc.pl.violin(adata, 'total_counts')
        p2 = sc.pl.violin(adata, "pct_counts_mt",show=False)
        pdf.savefig()
        p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",show=False)
        pdf.savefig()
        plt.close()
    
    def is_outlier(adata, metric: str, nmads: int):
        M = adata.obs[metric]
        outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
            np.median(M) + nmads * median_abs_deviation(M) < M
        )
        return outlier

    # log1p_total_counts, log1p_n_genes_by_counts and pct_counts_in_top_20_genes QC covariates each with a threshold of 5 MADs
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )

    # pct_counts_Mt is filtered with 3 MADs. Additionally, cells with a percentage of mitochondrial counts exceeding 8 % are filtered out.
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > 8
    )
    
    print(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
    stats["N.cells passing QC"]=adata.obs[Sample_key].value_counts()
    
    # Doublet Detection using scDblFinder R package
    import anndata2ri
    import logging
    
    import rpy2.rinterface_lib.callbacks as rcb
    import rpy2.robjects as ro
    
    rcb.logger.setLevel(logging.ERROR)
    ro.pandas2ri.activate()
    anndata2ri.activate()
    
    Rcode="""
    doublet <- function(mat){
    library(Seurat)
    library(scater)
    library(scDblFinder)
    library(BiocParallel)
    set.seed(123)
    sce = scDblFinder(
    SingleCellExperiment(
        list(counts=mat),
        ) 
    )
    doublet_score = sce$scDblFinder.score
    doublet_class = sce$scDblFinder.class
    ans=list("score"=doublet_score,"class"=doublet_class)
    return(ans)
    }
    """
    ro.r(Rcode)

    data_mat = adata.X.T
    ans=ro.r['doublet'](data_mat)
    
    adata.obs["scDblFinder_score"] = ans['score']
    adata.obs["scDblFinder_class"] = ans['class']
    adata.obs.scDblFinder_class.value_counts()

    adata=adata[adata.obs.scDblFinder_class=="singlet"]
    stats["N.singlet cells"]=adata.obs[Sample_key].value_counts()

    if batch_name is not None:
        adata.write(save_path+batch_name+"_qc.h5ad")
        stats.to_csv(save_path+batch_name+"_stats.csv")
    elif batch_key is not None:
        adata.write(save_path+adata.obs[batch_key].unique()[0]+"_qc.h5ad")
        stats.to_csv(save_path+adata.obs[batch_key].unique()[0]+"_stats.csv")
    else:
        print("Please give \"batch_name\" or \"batch_key\" ")