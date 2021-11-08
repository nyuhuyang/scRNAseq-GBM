# conda activate r4.0.3
library(Seurat)
library(magrittr)
library(pheatmap)
library(kableExtra)
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
object = readRDS("data/GBM_SCT_6_20210922.rds")
pred = readRDS(file = "output/GBM_6_20210922_singleR_pred.rds")

singlerDF = data.frame("label.fine" = pred$pruned.labels,
                       row.names = rownames(pred))
table(rownames(pred) == rownames(object@meta.data))
table(is.na(singlerDF$label.fine))
singlerDF$label.fine[is.na(singlerDF$label.fine)]= "unknown"

##############################
# adjust cell label
##############################
# combine cell types
singlerDF$label.major = singlerDF$label.fine

singlerDF[grep("CD4+",singlerDF$label.major),"label.major"] ="T_cells:CD4+"
singlerDF[grep("CD8+",singlerDF$label.major),"label.major"] ="T_cells:CD8+"
singlerDF[grep("B-cells",singlerDF$label.major),"label.major"] ="B_cells"
singlerDF[grep("Macrophages",singlerDF$label.major),"label.major"] ="Macrophages"
singlerDF$label.major %<>% gsub("Tregs","T_cells:regs",.)

singlerDF$label.major %<>% gsub("MEP|CLP|HSC|CMP|GMP|MPP","HSC/progenitors",.)
singlerDF$label.major %<>% gsub("DC|Megakaryocytes","other myloid cells",.)
singlerDF$label.major %<>% gsub("Skeletal muscle|Smooth muscle|Pericytes|Mesangial cells","muscle cells",.)
singlerDF$label.major %<>% gsub("mv Endothelial cells","Endothelial cells",.)


##############################
# process color scheme
##############################
table(colnames(object) == rownames(singlerDF))
object@meta.data %<>% cbind(singlerDF[,c("label.fine","label.major")])
object$label.major = singlerDF$label.major
UMAPPlot.1(object = object, label = T, label.repel = T,group.by = "label.major",
    no.legend = T,cols = Singler.colors,
    pt.size = 0.1,label.size = 5,alpha = 0.85,
    do.print = T,do.return = F,
    title ="labeling by blue_encode")

saveRDS(object, file = "data/GBM_SCT_6_20210922.rds")


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
