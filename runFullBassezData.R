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

setwd('/c4/home/brford/breastCancerProject/Bassez')


# Read in raw data
cohort1 <- readRDS("TNBC1.rds")
cohort1 <- CreateSeuratObject(counts = cohort1)
metadata.cohort1 <- read.csv("cohort1_meta.csv", row.names = 1)
cohort1@meta.data <- metadata.cohort1
cohort1@meta.data$orig.ident <- "cohort1"

cohort2 <-readRDS("TNBC.rds")
metadata.cohort2 <- read.csv("cohort2_meta.csv", row.names = 1)
cohort2@meta.data <- metadata.cohort2
cohort2@meta.data$orig.ident <- "cohort2"

bc.list <- list(cohort1, cohort2)

# Make qc table
qc.table <- RunAndVizQC(list.of.seurat.objs = bc.list, 
                        outputdir = "./"
)

# Add Pre QC number of Cells Value to table. Change defaults only if data quality suggests it is worth it.
qc.table$"QC_Filter_nFeature_RNA" <- 200 # default is 200
qc.table$"QC_Filter_nCount_RNA" <- 20000 # default is 20000
qc.table$"QC_Filter_percent.mt" <- 10 # default is 10

# rerun RunAndVizQC with info table provided to apply filters
bc.list <- RunAndVizQC(list.of.seurat.objs = bc.list, 
                       outputdir = "./",
                       info.file = qc.table
)


cohort1 <- bc.list[[1]]
cohort2 <- bc.list[[2]]

bc.down <- merge(cohort1, cohort2)

bc.down <- SCTransform(bc.down)
bc.down <- RunPCA(bc.down, assay = "SCT")

bc.down <- FindNeighbors(bc.down, dims = 1:40, reduction = "pca", assay = "SCT")
bc.down <- FindClusters(bc.down, resolution =0.5, cluster.name = "unintegrated_clusters", assay = "SCT")
bc.down <- RunUMAP(bc.down, dims = 1:40, reduction = "pca", reduction.name = "unintegrated_clusters", assay = "SCT")

# Make UMAPs with clusters and with comparisons 
pdf("UMAP.pdf")
DimPlot(bc.down, reduction = "unintegrated_clusters", group.by = c("orig.ident", "unintegrated_clusters"))
dev.off()

pdf("UMAPs_integrationCompare.pdf", width = 10, height = 5 )
DimPlot(bc.down, reduction = "unintegrated_clusters", group.by = c("orig.ident", "unintegrated_clusters"))
dev.off()

pdf("UMAPs_integrationCompare_patients.pdf", width = 10, height = 5 )
DimPlot(bc.down, reduction = "unintegrated_clusters", group.by = "patient_id")
dev.off()

pdf("UMAPs_integrationCompare_celltype.pdf", width = 6, height = 5 )
DimPlot(bc.down, reduction = "unintegrated_clusters", group.by = "cellType")
dev.off()

ref.genes <- c("CD3D", "CD3E", "CD2", "COL1A1", "DCN", "C1R", "LYZ", "CD68", "TYROBP", "CD24", "KRT19", "SCGB2A2", "EPCAM", "CD79A",  "MZB1", "MS4A1", "JCHAIN", "CLDN5", "FLT1", "RAMP2", "PECAM1", "CPA3", "TPSAB1", "TPSB2", "LILRA4", "CXCR3", "IRF7", "HBA1", "HBA2", "HBB")

avg.exp <- AverageExpression(bc.down, features = ref.genes)
avg.exp <- avg.exp[[2]]
pdf("AnnotationHeatmap.pdf")
pheatmap(avg.exp, scale = "row")
dev.off()

saveRDS(bc.down, "Full_Bassez_SeuratObj.RDS")