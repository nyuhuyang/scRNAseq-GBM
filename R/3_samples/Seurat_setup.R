########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.0.3
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","sctransform","harmony","magrittr"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)


########################################################################
#
#  1 Seurat Alignment
#
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# read sample summary list
# read sample summary list
df_samples <- readxl::read_excel("doc/20210822_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()

#======1.2 load  Seurat =========================
object = readRDS(file = "data/GBM_PBMC_3_20210908.rds")

table(df_samples$sample %in% object$orig.ident)
meta.data = object@meta.data
for(i in 1:length(df_samples$sample)){
    cells <- meta.data$orig.ident %in% df_samples$sample[i]
    print(df_samples$sample[i])
    print(table(cells))
    meta.data[cells,"Samples"] = df_samples$progress[i]
    meta.data[cells,"Type"] = df_samples$type[i]
    }
meta.data$orig.ident %<>% factor(levels = df_samples$sample)
table(rownames(object@meta.data) == rownames(meta.data))
table(colnames(object) == rownames(meta.data))
object@meta.data = meta.data
Idents(object) = "orig.ident"
#======1.6  Normalization =========================
set.seed(100)
object %<>% 

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(object) <- "integrated"

# Run the standard cca workflow for umap & tsne visualization
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(npcs = 100, verbose = FALSE)
jpeg(paste0(path,"ElbowPlot.jpeg"), units="in", width=10, height=7,res=600)
print(ElbowPlot(object,ndims = 100))
dev.off()
npcs =100
object %<>% JackStraw(num.replicate = 100,dims = npcs)
object %<>% ScoreJackStraw(dims = 1:npcs)

a = 0:9*10+1
for(i in 1:length(a)){
    b = a[i]+9
    jpeg(paste0(path,"JackStrawPlot_",a[i],"_",b,".jpeg"), units="in", width=10, height=7,res=600)
    print(JackStrawPlot(object, dims = a[i]:b))
    dev.off()
    Progress(i,length(a))
}
npcs = 38
jpeg(paste0(path,"S1_RunHarmony.jpeg"), units="in", width=10, height=7,res=600)
system.time(object %<>% RunHarmony(group.by = "orig.ident", dims.use = 1:npcs))
dev.off()

object %<>% RunUMAP(reduction = "harmony", dims = 1:npcs)
system.time(object %<>% RunTSNE(reduction = "harmony", dims = 1:npcs))

object %<>% FindNeighbors(reduction = "harmony",dims = 1:2)
object %<>% FindClusters(resolution = 0.8)

saveRDS(object, file = "data/GBM_PBMC_3_20210908.rds")