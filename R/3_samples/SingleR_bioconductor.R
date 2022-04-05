#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0.3 linux
library(SingleR)
library(Seurat)
library(celldex)
library(SingleCellExperiment)
library(magrittr)
library(data.table)

source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

# ====== load single cell =============
object <- readRDS("data/GBM_PBMC_3_20210908.rds")

sce <- SingleCellExperiment(list(logcounts=object[["RNA"]]@data),
                                colData=DataFrame(object@meta.data))
rm(object);GC()

references = c("blue_encode","VHove2019","PAntunes2021","VHove2019_v2")[4]
# ====== load reference =============
if(references == "blue_encode"){
    blue_encode <- BlueprintEncodeData()
    rownames(blue_encode) = Hmisc::capitalize(tolower(rownames(blue_encode)))
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(blue_encode)
    ))
    length(common)
    table(blue_encode$label.fine)
    system.time(trained <- trainSingleR(ref = blue_encode[common,],
                                        labels=blue_encode$label.fine))
    system.time(pred <- classifySingleR(sce[common,], trained))
    # elapsed 4872.846 sec
    saveRDS(object = pred, file = "output/GBM_PBMC_20210907_singleR_pred.rds")
}

if(references == "VHove2019"){
    counts <- Read10X(data.dir = "data/annotation_references/GSE128854/filtered_feature_bc_matrix_MouseTransplantedGBM/filtered_gene_bc_matrices/mm10")
    meta.data = read.csv("data/annotation_references/GSE128854/annot_Mouse_GBM_Full.csv")
    table(colnames(counts) %in% meta.data$cell)
    rownames(meta.data) = meta.data$cell
    counts = counts[,meta.data$cell]
    libsizes <- colSums(counts)
    size.factors <- libsizes/mean(libsizes)
    VHove2019 <- SingleCellExperiment(list(logcounts=log1p(t(t(counts)/size.factors))),
                                 colData=DataFrame(meta.data))
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(VHove2019)
    ))
    length(common)
    table(VHove2019$cluster)
    system.time(trained <- trainSingleR(ref = VHove2019[common,],
                                        labels=VHove2019$cluster))
    system.time(pred <- classifySingleR(sce[common,], trained))
    # elapsed 4872.846 sec
    saveRDS(object = pred, file = "output/GBM_PBMC_20220329_VHove2019_singleR_pred.rds")
}

if(references == "PAntunes2021"){
    meta.data <- fread("data/annotation_references/GSE163120/GSE163120_annot.Mouse.GBM.KO1_2_3.WT1_2_3.csv.gz")
    counts <- fread("data/annotation_references/GSE163120/GSE163120_Mouse.GBM.KO1_2_3.WT1_2_3.filtered.gene.bc.matrix.csv.gz")
    rownames(counts) = counts$V1
    counts[,V1:=NULL]
    table(colnames(counts) %in% meta.data$cell)
    rownames(meta.data) = meta.data$cell
    counts = setcolorder(counts,meta.data$cell)
    libsizes <- colSums(counts)
    size.factors <- libsizes/mean(libsizes)
    PAntunes2021 <- SingleCellExperiment(list(logcounts=log1p(t(t(counts)/size.factors))),
                                      colData=DataFrame(meta.data))
    rownames(PAntunes2021) = rownames(counts)
    
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(PAntunes2021)
    ))
    length(common)
    table(PAntunes2021$cluster)
    system.time(trained <- trainSingleR(ref = PAntunes2021[common,],
                                        labels=PAntunes2021$cluster))
    system.time(pred <- classifySingleR(sce[common,], trained))
    # elapsed 4872.846 sec
    saveRDS(object = pred, file = "output/GBM_PBMC_20220401_PAntunes2021_singleR_pred.rds")
}


if(references == "VHove2019_v2"){
    counts <- Read10X(data.dir = "data/annotation_references/GSE128854/filtered_gene_bc_matrices_WT_wholeBrain/mm10")
    meta.data = read.csv("data/annotation_references/GSE128854/annot_K10.csv")
    colnames(counts) %<>% gsub("-.*","",.)
    counts = counts[,meta.data$cell]
    table(colnames(counts) %in% meta.data$cell)
    meta.data$cell[!(meta.data$cell %in% colnames(counts) )]
    
    rownames(meta.data) = meta.data$cell
    counts = counts[,meta.data$cell]
    libsizes <- colSums(counts)
    size.factors <- libsizes/mean(libsizes)
    VHove2019 <- SingleCellExperiment(list(logcounts=log1p(t(t(counts)/size.factors))),
                                      colData=DataFrame(meta.data))
    common <- Reduce(intersect, list(rownames(sce),
                                     rownames(VHove2019)
    ))
    length(common)
    table(VHove2019$cluster)
    system.time(trained <- trainSingleR(ref = VHove2019[common,],
                                        labels=VHove2019$cluster))
    system.time(pred <- classifySingleR(sce[common,], trained))
    # elapsed 4872.846 sec
    saveRDS(object = pred, file = "output/GBM_PBMC_20220404_VHove2019_v2_singleR_pred.rds")
}
