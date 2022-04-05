#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.1.1 linux
library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(magrittr)

source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")


#=========== read Van Hove 2019 brainimmuneatlas counts and annotations ================================
# wget https://www.brainimmuneatlas.org/data_files/toDownload/filtered_feature_bc_matrix_MouseTransplantedGBM.zip
# wget https://www.brainimmuneatlas.org/data_files/toDownload/annot_Mouse_GBM_Full.csv

counts <- Read10X(data.dir = "data/annotation_references/GSE128854/filtered_feature_bc_matrix_MouseTransplantedGBM/filtered_gene_bc_matrices/mm10")
meta.data = read.csv("data/annotation_references/GSE128854/annot_Mouse_GBM_Full.csv")
table(colnames(counts) %in% meta.data$cell)
rownames(meta.data) = meta.data$cell
counts = counts[,meta.data$cell]
VHove <- CreateSeuratObject(counts = counts,min.cells = 0,min.features = 0,meta.data = meta.data)
VHove %<>% NormalizeData()
VHove %<>% FindVariableFeatures()
VHove %<>% ScaleData()
VHove %<>% RunPCA()
VHove %<>% RunUMAP(dims = 1:15)
UMAPPlot(VHove,group.by = "cluster")
VHove[["umap"]]@cell.embeddings[,"UMAP_1"] = VHove$x
VHove[["umap"]]@cell.embeddings[,"UMAP_2"] = VHove$y

# wget https://www.brainimmuneatlas.org/data_files/toDownload/filtered_gene_bc_matrices_WT_wholeBrain.zip
# wget https://www.brainimmuneatlas.org/data_files/toDownload/annot_K10.csv
# unzip filtered_gene_bc_matrices_WT_wholeBrain.zip
# mv filtered_gene_bc_matrices filtered_gene_bc_matrices_WT_wholeBrain


