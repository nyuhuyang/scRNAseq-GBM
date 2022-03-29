#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0.3 linux
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(magrittr)

source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

# ====== load single cell =============
object <- readRDS("data/GBM_SCT_6_20210922.rds")

sce <- SingleCellExperiment(list(logcounts=object[["SCT"]]@data),
                                colData=DataFrame(object@meta.data))
rm(object);GC()

# ====== load reference =============
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
saveRDS(object = pred, file = "output/GBM_6_20210922_singleR_pred.rds")
