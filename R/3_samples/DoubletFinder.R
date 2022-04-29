########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.0.3
library(Seurat)
library(dplyr)
library(cowplot)
library(magrittr)
library(DoubletFinder)
require(fields)
require(parallel)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("R/util.R")

# run below line in cluster

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)
########################################################################
#
#  2. DoubletFinder
#
# ######################################################################

# samples

object = readRDS(file = "data/MCL_61_20220331.rds")

(samples = unique(object$orig.ident))
object_list <- SplitObject(object,split.by = "orig.ident")
rm(object);GC()
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
Sys.setenv("OMP_NUM_THREADS" = 16)

npcs = 100
sweep.res_list <- list()
for (i in 1:length(object_list)) {
    sweep.res_list[[i]] <- paramSweep_v4(object_list[[i]], PCs = 1:npcs, sct = T)
    Progress(i,length(object_list))
}
saveRDS(sweep.res_list,file = "output/MCL_61_20220331_sweep.res_list.rds")
sweep.res_list = readRDS(file = "output/MCL_61_20220331_sweep.res_list.rds")
sweep_list <- lapply(sweep.res_list, function(x) summarizeSweep(x, GT = FALSE))
bcmvn_list <- lapply(sweep_list,find.pK)

(maximal_pk <- sapply(bcmvn_list,function(x) {
    as.numeric(as.character(x[find.localMaxima(x$BCmetric),"pK"]))
    }))
maximal_pk

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
for(i in 1:length(object_list)){
    print(paste("processing",unique(object_list[[i]]$orig.ident)))
    homotypic.prop <- modelHomotypic(object_list[[i]]@meta.data$cell.types)
    nExp_poi <- round(Multiplet_Rate(object_list[[i]])*length(colnames(object_list[[i]])))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    object_list[[i]] <- doubletFinder_v3(object, PCs = 1:59,
                                         pN = 0.25, pK = maximal_pk,
                                         nExp = nExp_poi, reuse.pANN = FALSE, sct = F,
                                         annotations = object$label.fine)
    object_list[[i]] <- doubletFinder_v3(object_list[[i]], PCs = 1:50,
                                         pN = 0.25, pK = maximal_pk[i],
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = grep("pANN",colnames(object_list[[i]]@meta.data),value = T),
                                         sct = FALSE)
    colName = colnames(object_list[[i]]@meta.data)
    colName[grep("DF.classifications",colName)] = c("Low_confident_doublets",
                                                    "High_confident_doublets")
    colnames(object_list[[i]]@meta.data) = colName
    Progress(i,length(object_list))
}

for(i in 1:length(object_list)){
    object_list[[i]]@meta.data$row.names = rownames(object_list[[i]]@meta.data)
}
meta.data_list <- lapply(object_list, function(x) {
    temp <- x@meta.data
    temp$row.names = rownames(temp)
    return(temp)
    })
meta.data = bind_rows(meta.data_list)
rownames(meta.data) = meta.data$row.names

object = readRDS(file = "data/MCL_61_20220331.rds")
meta.data = meta.data[rownames(object@meta.data),]
meta.data$doublets = gsub("Doublet","Doublet-Low Confidence",meta.data$Low_confident_doublets)
meta.data[meta.data$High_confident_doublets %in% "Doublet","doublets"] = "Doublet-High Confidence"
meta.data = cbind(object@meta.data,meta.data$doublets)
colnames(meta.data)[ncol(meta.data)] = "Doublets"
table(meta.data$Doublets)
object@meta.data = meta.data
saveRDS(object,file=paste0("data/MCL_61_20220331.rds"))

UMAPPlot.1(object, group.by = "Doublets",cols = c("red","orange","black"),raster=FALSE,
           title = "Singlets and possible Doublets", do.print = T,pt.size = 0.3)
