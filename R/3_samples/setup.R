library(Seurat) # Seurat 4
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

object <- readRDS("data/GBM_PBMC_3_20210908.rds")

object$orig.ident = gsub("_.*","",rownames(object@meta.data))

s.genes <- Hmisc::capitalize(tolower(cc.genes$s.genes)) %>% FilterGenes(object,.)
g2m.genes <- Hmisc::capitalize(tolower(cc.genes$g2m.genes)) %>% FilterGenes(object,.)
object %<>% CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


Idents(object) = "Samples"
g2 <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(features){
    VlnPlot(object = object, features = features, ncol = 1, pt.size = 0.01)+
        theme(axis.text.x = element_text(size=15,angle = 0,hjust = 0.5),legend.position="none")
})

jpeg(paste0(path,"S1_nGene.jpeg"), units="in", width=10, height=7,res=600)
print(g2[[1]]+ggtitle("nFeature_RNA after filteration")+
                    scale_y_log10(limits = c(100,10000))+
                    theme(plot.title = element_text(hjust = 0.5)))
dev.off()
jpeg(paste0(path,"S1_nUMI.jpeg"), units="in", width=10, height=7,res=600)
print(g2[[2]]+ggtitle("nCount_RNA after filteration")+
                    scale_y_log10(limits = c(500,100000))+
                    theme(plot.title = element_text(hjust = 0.5)))
dev.off()
jpeg(paste0(path,"S1_mito.jpeg"), units="in", width=10, height=7,res=600)
print(g2[[3]]+ggtitle("mito % after filteration")+
                    ylim(c(0,15))+
                    theme(plot.title = element_text(hjust = 0.5)))
dev.off()