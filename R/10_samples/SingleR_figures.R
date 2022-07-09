# conda activate r4.1.1
library(Seurat)
library(magrittr)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================

##############################
# create singleR data frame
###############################
object = readRDS("data/GBM_SCT_10_20220707.rds")
pred = readRDS(file = "output/GBM_10_20220707_blue_encode_singleR_pred.rds")
pred = readRDS(file = "output/GBM_10_20220707_VHove2019_singleR_pred.rds")
pred = readRDS(file = "output/GBM_10_20220707_PAntunes2021_singleR_pred.rds")
pred = readRDS(file = "output/GBM_10_20220707_VHove2019_v2_singleR_pred.rds")


##############################
# adjust final_celltype cell label
##############################
# combine cell types
singlerDF = data.frame("label.fine" = pred$pruned.labels,
                       row.names = rownames(pred))
table(rownames(pred) == rownames(object@meta.data))
table(is.na(singlerDF$label.fine))
singlerDF$label.fine[is.na(singlerDF$label.fine)]= "unknown"

object@meta.data %<>% cbind(label.fine = singlerDF[,c("label.fine")])


singlerDF$label.main = singlerDF$label.fine


singlerDF[grep("CD4+",singlerDF$label.main),"label.main"] ="T_cells:CD4+"
singlerDF[grep("CD8+",singlerDF$label.main),"label.main"] ="T_cells:CD8+"
singlerDF[grep("B-cells",singlerDF$label.main),"label.main"] ="B_cells"

singlerDF$label.main %<>% gsub("Tregs","T_cells:regs",.)
singlerDF$label.main %<>% gsub("Adipocytes|Adipocytes|Astrocytes|Fibroblasts|mv Endothelial cells|Endothelial cells","Nonhematopoietic cells",.)
singlerDF$label.main %<>% gsub("Mesangial cells|Neurons|Skeletal muscle","Nonhematopoietic cells",.)
saveRDS(object@meta.data, file = "output/GBM_10_20220707_metadata.rds")


##############################
# adjust VHove2019 cell label
##############################
object$VHove2019_label = singlerDF$label.fine
UMAPPlot.1(object, group.by = "VHove2019_label",label = T, label.repel = T,do.print = T)

object$VHove2019_v2 = singlerDF$label.fine

##############################
# adjust PAntunes2021 cell label
##############################
object$PAntunes2021 = singlerDF$label.fine
UMAPPlot.1(object, group.by = "PAntunes2021",label = T, label.repel = T,do.print = T)


##############################
# process color scheme
##############################
table(colnames(object) == rownames(singlerDF))
object@meta.data %<>% cbind(singlerDF[,c("label.fine")])
object$cell.label = singlerDF$cell.label
UMAPPlot.1(object = object, label = T, label.repel = T,group.by = "cell.label",
           no.legend = T,cols = Singler.colors,
           pt.size = 0.1,label.size = 5,alpha = 0.85,
           do.print = T,do.return = F,
           title ="labeling by blue_encode")

saveRDS(object, file = "data/GBM_PBMC_3_20210908.rds")
