# Utility functions for single-cell analysis
library(Seurat)
library(tidyverse)
library(SingleR)
library(celldex)

load_sc_data <- function(data_path) {
  # Load 10X data
  data <- Read10X(data_path)
  seurat_obj <- CreateSeuratObject(counts = data, project = "scRNA")
  return(seurat_obj)
}

quality_control <- function(seurat_obj, mt_pattern = "^MT-", rb_pattern = "^RP[SL]") {
  # Calculate QC metrics
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = rb_pattern)
  
  return(seurat_obj)
}

plot_qc_metrics <- function(seurat_obj) {
  # Plot QC metrics
  p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  return(list(violin = p1, mt_scatter = p2, feature_scatter = p3))
}
