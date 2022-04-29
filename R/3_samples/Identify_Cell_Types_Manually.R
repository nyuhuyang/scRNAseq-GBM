library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(magrittr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
#===========

object <- readRDS("data/GBM_PBMC_3_20210908.rds")
object@meta.data$Samples %<>% factor(levels = c("No tumor","LGG","HGG"))

s.genes <- Hmisc::capitalize(tolower(cc.genes$s.genes))
g2m.genes <- Hmisc::capitalize(tolower(cc.genes$g2m.genes))
s.genes %<>% gsub("Mlf1ip","Cenpu",.)
g2m.genes %<>% plyr::mapvalues(from = c("Fam64a", "Hn1"),
                               to = c("Pimreg","Jpt1"))

object %<>% CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
resolutions = c(seq(0.8,1.5, by = 0.1),seq(2,5, by = 0.5))
for(i in 1:length(resolutions)){
    object %<>% FindClusters(resolution = resolutions[i])
    object@meta.data[,paste0("RNA_snn_res.",resolutions[i])] %<>% as.factor()
}

df_samples <- readxl::read_excel("doc/20220421_RNA_ssn_res.1.5.xlsx")
df_samples = as.data.frame(df_samples)

for (label in c("Cell_type","Cell type layer 2","Cell type layer 1")){
    object@meta.data[,label] = as.character(plyr::mapvalues(object@meta.data[,"RNA_snn_res.1.5"],
                                          from = df_samples$RNA_snn_res.1.5,
                                          to = df_samples[,label]))
}

change = object$final_celltype %in% c("Activated Microglia","CD8 T cells", "CD4 T cells","NK cells")
object@meta.data[change,"Cell_type"] = as.character(object@meta.data[change,"final_celltype"])
change = object$final_celltype %in% c("NK cells")
object@meta.data[change,"Cell type layer 2"] = as.character(object@meta.data[change,"final_celltype"])


length(unique(object$Cell_type))
source("https://raw.githubusercontent.com/nyuhuyang/multiomicShinyCell/master/shinyApp_ad/util_palette.R")
source("https://raw.githubusercontent.com/nyuhuyang/multiomicShinyCell/master/shinyApp_ad/util.R")


object %<>% AddMetaColor(label = "Cell_type", colors = color_generator("colorBlind",length(unique(object$Cell_type))))

saveRDS(object@meta.data, file = "data/GBM_PBMC_3_20210908_metadata.rds")
