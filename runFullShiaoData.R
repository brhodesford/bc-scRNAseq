# check if packages installed and install if not
# Function to check if a package is installed, and install if missing
if (!requireNamespace("Seurat", quietly = TRUE)) {
  remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
}
if (!requireNamespace("SeuratData", quietly = TRUE)) {
  devtools::install_github('satijalab/seurat-data')
}

if (!requireNamespace("SeuratWrappers", quietly = TRUE)) {
  remotes::install_github('satijalab/seurat-wrappers')
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}


if (!requireNamespace("patchwork", quietly = TRUE)) {
  devtools::install_github("thomasp85/patchwork", repos = "https://cloud.r-project.org")
}


if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate", repos = "https://cloud.r-project.org")
}


if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}

if (!requireNamespace("vidger", quietly = TRUE)) {
  BiocManager::install("vidger")
}

if (!requireNamespace("factoextra", quietly = TRUE)) {
  install.packages("factoextra", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("plyr", quietly = TRUE)) {
  install.packages("plyr", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}


# Read in libraries

library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(reticulate)
library(pheatmap)
library(DESeq2)
library(vidger)
library(factoextra)
library(plyr)
library(org.Hs.eg.db)
library(biomaRt)

options(future.globals.maxSize = 2e10)

source("/c4/home/brford/general_scripts/working_with_10x_data_forRMD.R")

setwd('/c4/home/brford/breastCancerProject/Shiao')


# Read in raw data
shiao.imm.seurat <- readRDS("shiao_imm_Seurat.rds")
shiao.nonimm.seurat <- readRDS("shiao_nonimm_Seurat.rds")

shiao.list <- list(shiao.imm.seurat, shiao.nonimm.seurat)

# Make qc table
qc.table <- RunAndVizQC(list.of.seurat.objs = shiao.list, 
                        outputdir = "./"
)

# Add Pre QC number of Cells Value to table. Change defaults only if data quality suggests it is worth it.
qc.table$"QC_Filter_nFeature_RNA" <- 200 # default is 200
qc.table$"QC_Filter_nCount_RNA" <- 20000 # default is 20000
qc.table$"QC_Filter_percent.mt" <- 10 # default is 10

# rerun RunAndVizQC with info table provided to apply filters
shiao.list <- RunAndVizQC(list.of.seurat.objs = shiao.list, 
                          outputdir = "./",
                          info.file = qc.table
)


shiao.list <- lapply(shiao.list, function(x){
  x <- SCTransform(x)
  x <- RunPCA(x, assay = "SCT")
  return(x)
})

shiao.list[[1]] <- FindNeighbors(shiao.list[[1]], dims = 1:30, reduction = "pca", assay = "SCT")
shiao.list[[1]] <- FindClusters(shiao.list[[1]], resolution = 1, cluster.name = "umap", assay = "SCT")
shiao.list[[1]] <- RunUMAP(shiao.list[[1]], dims = 1:30, reduction = "pca", reduction.name = "umap", assay = "SCT")

shiao.list[[2]] <- FindNeighbors(shiao.list[[2]], dims = 1:15, reduction = "pca", assay = "SCT")
shiao.list[[2]] <- FindClusters(shiao.list[[2]], resolution =0.5, cluster.name = "umap", assay = "SCT")
shiao.list[[2]] <- RunUMAP(shiao.list[[2]], dims = 1:15, reduction = "pca", reduction.name = "umap", assay = "SCT")

# Make UMAPs with clusters and with comparisons 
pdf("UMAP.pdf")
lapply(shiao.list, function(x){
  DimPlot(x, reduction = "unintegrated_clusters")
})
dev.off()

pdf("UMAPs_integrationCompare.pdf", width = 10, height = 5 )
DimPlot(bc.down, reduction = "umap", group.by = c("orig.ident", "umap"))
dev.off()

pdf("UMAPs_integrationCompare_patients.pdf", width = 10, height = 5 )
lapply(shiao.list, function(x){
  DimPlot(x, reduction = "umap", group.by = "patient_id")
})
dev.off()

pdf("UMAPs_integrationCompare_celltype.pdf", width = 6, height = 5)
lapply(shiao.list, function(x){
  DimPlot(bc.down, reduction = "umap", group.by = "cellType")
})
dev.off()

imm.genes <- c("CD3D", "CD3E", "CD2", "LYZ", "CD68", "TYROBP", "CD79A",  "MZB1", "MS4A1", "JCHAIN", "TPSAB1", "TPSB2", "CPA3", "NCAM1", "KLRB1")
nonimm.genes <- c("COL1A1", "DCN", "C1R", "CD24", "KRT19", "SCGB2A2", "EPCAM", "CLDN5", "FLT1", "RAMP2", "PECAM1", "TPSAB1", "TPSB2", "CPA3")

imm.exp <- AverageExpression(shiao.down.list[[1]], features = imm.genes)
imm.exp <- imm.exp[[2]]

nonimm.exp <- AverageExpression(shiao.down.list[[2]], features = nonimm.genes)
nonimm.exp <- nonimm.exp[[2]]

pdf("AnnotationHeatmap_Imm.pdf")
print(pheatmap(imm.exp, scale = "row"))
dev.off()

pdf("AnnotationHeatmap_Nonimm.pdf")
print(pheatmap(nonimm.exp, scale = "row"))
dev.off()

saveRDS(shiao.list, "FullShiaoData_clustered.rds")