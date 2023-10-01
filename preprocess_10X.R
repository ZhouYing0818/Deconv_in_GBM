library(Seurat)
library(celda) # normalization and preprocessing
# dataset integration
library(harmony) 
library(SingleR)
# plots
library(ggplot2)
library(reshape2)
library(cowplot)

file_path1="/data/newdata/deconv/anno_test/PRJNA579593_Bhaduri_2020/count_mat/10X/SF11159_1/outs/filtered_feature_bc_matrix.h5"
file_path2="/data/newdata/deconv/anno_test/PRJNA579593_Bhaduri_2020/count_mat/10X/SF11159_2/outs/filtered_feature_bc_matrix.h5"
sce1 <- Read10X_h5(filename = file_path1)
sce1 <- CreateSeuratObject(counts = sce1)
sce2 <- Read10X_h5(filename = file_path2)
sce2 <- CreateSeuratObject(counts = sce2)
sce<-merge(sce1,sce2,add.cell.ids = c("SF11159_1", "SF11159_2"), project = "SF11159")

## Generic QC function
qc.features <- function(
    obj, # a seurat object, typically the output of as.Seurat
    cc = FALSE
) {
  
  # calculate mitochondrial percetage
  obj[["percent.mito"]] <- PercentageFeatureSet(object = obj, pattern = "^MT-") 
  
  # calculate ribosomal percentage
  obj[["percent.rps"]] <- PercentageFeatureSet(object = obj, pattern = "^RPS") 
  obj[["percent.rpl"]] <- PercentageFeatureSet(object = obj, pattern = "^RPL") 
  
  # calculate MALAT1 (nuclear marker) percentage
  obj[["percent.malat1"]] <- PercentageFeatureSet(object = obj, pattern = "MALAT1") 
  
  # a cell-cycle score
  if(cc) obj <- CellCycleScoring(obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
  
  # return
  return(obj)
}

VlnPlot(sce, features = c("nFeature_RNA", 
                          "nCount_RNA", 
                          "percent.mito",
                          "percent.rps",
                          "percent.rpl"), 
        ncol = 3)

plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


## Functions: normalisation ----

sce <- subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mito < 20)
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- SCTransform(sce, assay = "RNA", return.only.var.genes = FALSE, verbose = FALSE)
sce <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(sce)

all.genes <- rownames(sce)
sce <- ScaleData(sce, features = all.genes)

sce <- RunPCA(sce, features = VariableFeatures(object = sce))
VizDimLoadings(sce, dims = 1:2, reduction = "pca")
DimPlot(sce, reduction = "pca")
DimHeatmap(sce, dims = 1:20, cells = 500, balanced = TRUE)
ElbowPlot(sce)

sce <- FindNeighbors(sce, dims = 1:10)
sce <- FindClusters(sce, resolution = 0.5)

sce <- RunUMAP(sce, dims = 1:10)
DimPlot(sce, reduction = "umap")

norm.bySeurat <- function(
    obj, # a Seurat object
    min.nUMI = 500, # remove cells/nuclei with nUMI less than this
    max.percent.mito = 20, # remove cells/nuclei with percent.mito greater than this
    max.nUMI.percentile = 0.998 # remove cells/nuclei with nUMI greater than this percentile of all cells/nuclei
) {
  
  # subset cells
  topCount <- quantile(obj@meta.data$nCount_RNA, probs = max.nUMI.percentile)
  obj <- subset(x = obj, subset = (nCount_RNA > min.nUMI & nCount_RNA < topCount & percent.mito < max.percent.mito))
  
  # obj <- subset(x = obj, subset = nCount_RNA > 500 & percent.mito < 5)
  
  # normalise expression levels
  #obj <- NormalizeData(object = obj, normalization.method = "LogNormalize", scale.factor = 10000) # standard parameters for Seurat
  
  # find variable genes (i.e. features)
  obj <- FindVariableFeatures(object = obj, selection.method = "vst", nfeatures = 2000)
  
  # scale data, and regress out covariates of nUMI and percent.mito. takes ~2 min. 
  # note that FindVariableFeatures returns identical results whether or not the data are scaled first
  obj <- SCTransform(ref_mat, assay = "RNA", return.only.var.genes = FALSE, verbose = FALSE)
  
  # output
  return(obj)
}

