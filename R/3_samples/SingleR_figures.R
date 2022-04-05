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
object = readRDS("data/GBM_PBMC_3_20210908.rds")
pred = readRDS(file = "output/GBM_PBMC_20210907_singleR_pred.rds")
pred = readRDS(file = "output/GBM_PBMC_20220329_VHove2019_singleR_pred.rds")
pred = readRDS(file = "output/GBM_PBMC_20220401_PAntunes2021_singleR_pred.rds")
pred = readRDS(file = "output/GBM_PBMC_20220404_VHove2019_v2_singleR_pred.rds")
singlerDF = data.frame("label.fine" = pred$pruned.labels,
                       row.names = rownames(pred))
table(rownames(pred) == rownames(object@meta.data))
table(is.na(singlerDF$label.fine))
singlerDF$label.fine[is.na(singlerDF$label.fine)]= "unknown"

##############################
# adjust final_celltype cell label
##############################
# combine cell types
singlerDF$final_celltype = object$final_celltype

singlerDF[grep("^Microglia$",singlerDF$final_celltype),"label.fine"] ="Microglia"
singlerDF[grep("^Activated Microglia$",singlerDF$final_celltype),"label.fine"] ="Activated Microglia"
singlerDF$cell.label = singlerDF$label.fine

singlerDF[grep("CD4+",singlerDF$cell.label),"cell.label"] ="T_cells:CD4+"
singlerDF[grep("CD8+",singlerDF$cell.label),"cell.label"] ="T_cells:CD8+"
singlerDF[grep("B-cells",singlerDF$cell.label),"cell.label"] ="B_cells"

singlerDF$cell.label %<>% gsub("Tregs","T_cells:regs",.)
singlerDF$cell.label %<>% gsub("Adipocytes|Adipocytes|Astrocytes|Fibroblasts|mv Endothelial cells|Endothelial cells","Nonhematopoietic cells",.)
singlerDF$cell.label %<>% gsub("Mesangial cells|Neurons|Skeletal muscle","Nonhematopoietic cells",.)

singlerDF %<>% cbind(object[["umap"]]@cell.embeddings)
singlerDF[singlerDF$label.fine %in% "Macrophages" & singlerDF$UMAP_1 < -6,"cell.label"] ="Macrophages cluster 1"


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
object@meta.data %<>% cbind(singlerDF[,c("label.fine","cell.label")])
object$cell.label = singlerDF$cell.label
UMAPPlot.1(object = object, label = T, label.repel = T,group.by = "cell.label",
    no.legend = T,cols = Singler.colors,
    pt.size = 0.1,label.size = 5,alpha = 0.85,
    do.print = T,do.return = F,
    title ="labeling by blue_encode")

saveRDS(object, file = "data/GBM_PBMC_3_20210908.rds")


# by barplot
cell_Freq <- table(object$label.fine) %>% as.data.frame
cell_Freq$Percent <- prop.table(cell_Freq$Freq) %>% round(digits = 2) %>% scales::percent()
cols = ExtractMetaColor(object)
cell_Freq$cols = cols[cell_Freq$Var1]
cell_Freq = cell_Freq[order(cell_Freq$Var1),]

cell_Freq = cell_Freq[order(cell_Freq$Freq,decreasing = T),]
cell_Freq$Var1 %<>% factor(levels = as.character(cell_Freq$Var1))
colnames(cell_Freq)[1:2] = c("Cell_Type", "Cell_Number")

jpeg(paste0(path,"cell_type_numbers.jpeg"), units="in", width=6, height=6,res=600)
ggbarplot(cell_Freq, "Cell_Type", "Cell_Number",
          fill = "Cell_Type", color = "black",xlab = "",
          palette = cell_Freq$col,x.text.angle = 90,
          ylab = "Cell Number",
          label = cell_Freq$Percent,
          lab.size = 3,
          sort.val = "desc",
          width = 1, size = 0.5,
          title = "Numbers of cell types in total 6 samples")+NoLegend()+
    theme(plot.title = element_text(hjust = 0.5,size=15),
          axis.text.x = element_text(vjust = 0.5))
dev.off()
