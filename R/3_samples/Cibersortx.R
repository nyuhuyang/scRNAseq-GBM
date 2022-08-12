########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################
#conda activate r4.0.3 
# devtools::install_github('satijalab/seurat-data') #3.1.5.9900
# remotes::install_github("mojaveazure/seurat-disk")
invisible(lapply(c("Seurat","dplyr","cowplot","magrittr","tidyr","tibble","data.table"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
save_path <- "Yang/Deconvolution/"
if(!dir.exists(save_path))dir.create(save_path, recursive = T)

object <- readRDS("data/GBM_PBMC_3_20210908.rds")

for(s in c("No tumor","LGG","HGG")){
    print(s)
    sub_object = subset(object,Samples == s)
    data = sub_object[["RNA"]]@data
    cell_ID = sub_object$Cell_type
    colnames(data) = cell_ID
    RowSums <- rowSums(data)
    data = data[RowSums > 0 ,]
    print(dim(data))
    #ColSums <- colSums(data)
    #scale_data = sweep(data, 2, ColSums,"/")*100000
    #ColSums <- colSums(scale_data)
    #range(ColSums)
    
    GeneSymbol = data.frame("GeneSymbol" =  rownames(data))
    mat = as.data.table(cbind(GeneSymbol,data),check.names = FALSE)

    fwrite(mat,file = paste0(save_path,"scRNAseq_GBM_",s,"_reference.txt"),quote = FALSE,
           row.names = FALSE,sep = "\t")
}

#========================= aggregate data with k cells =====================
# generate expression txt file for Cibersortx 
#' @param object Seurat object
#' @param k an integer for the number of folds. createFolds argment
#' @param do.return TRUE/FALSE
#' @param group.by aggregate by. column name in meta.data
#' @example PrepareCiberSortx(object, k = 10, group.by = "Cell_type")
PrepareCiberSortx <- function(object, k = 10, do.return = FALSE, group.by = NULL,
                              file.name = NULL,...){
    try(if(!(group.by %in% colnames(object@meta.data))) stop("group.by is not within meta.data column"))
    
    set.seed(201)
    if(k > 1){
        #Split object meta.data into 2 group =====
        split_meta.data <- split(object@meta.data, f = object@meta.data[,group.by])
        #Split each meta.data into k group =====
        for(i in 1:length(split_meta.data)){
            meta_index <- caret::createDataPartition(split_meta.data[[i]][,group.by],
                                                     groups = min(k, length(split_meta.data[[i]][,group.by])))
            for(n in 1:length(meta_index)){
                split_meta.data[[i]][meta_index[[n]],"GSEA"] = 
                    paste(split_meta.data[[i]][meta_index[[n]],ident],
                          n, sep = "_")
            }
        }
    }
    meta.data = bind_rows(split_meta.data)
    object@meta.data = meta.data[rownames(object@meta.data),]
    
    Idents(object) = "GSEA"
    if(!is.null(continuous.label)) {
        Idents(object) %<>% factor(levels = paste0(rep(continuous.label,each = k),
                                                   "_",rep(1:k)))
    }
    
    
    if(k == 1 & !is.null(continuous.label)) {
        Idents(object) %<>% factor(levels = continuous.label)
    }
    
    
    print("#====Calculate Average Expression======")
    GSEA_expr <- AverageExpression(object, assays = DefaultAssay(object))
    GSEA_expr = GSEA_expr[[1]]
    GSEA_name <- data.frame("NAME" = rownames(GSEA_expr),
                            "DESCRIPTION" = rep(NA,nrow(GSEA_expr)),
                            stringsAsFactors = F)
    GSEA_expr <- cbind.data.frame(GSEA_name,GSEA_expr)
    
    samples <- gsub("_([0-9]+).*$", "", colnames(GSEA_expr)[-c(1,2)])
    
    if(is.null(continuous.label)) {
        cls_list <- list(c(length(samples),length(unique(samples)), 1),
                         paste(c("#",unique(samples)), collapse = " "),
                         paste(samples, collapse = " "))
    } else 
        if(all(continuous.label %in% unique(object@ident))){
            numeric <- match(samples,continuous.label)
            cls_list <- list("#numeric",
                             paste(c("#",unique(samples)), collapse = "."),
                             paste(numeric, collapse = " "))
        } else stop("Incorrect continuous.label!")
    
    if(rownames(GSEA_expr)[1] == 
       Hmisc::capitalize(tolower(rownames(GSEA_expr)[1]))){
        if(Mouse2Human == "biomaRt"){
            print("#====Replace gene names using BiomaRt ======")
            rownames.GSEA_expr = rownames(GSEA_expr)
            human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
            mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
            
            genesV2 = biomaRt::getLDS(attributes = c("mgi_symbol"), #filters = "mgi_symbol",
                                      values = rownames(GSEA_expr) , mart = mouse,
                                      attributesL = c("hgnc_symbol"), 
                                      martL = human, uniqueRows=T)
            rm = duplicated(genesV2[,1])
            genesV2 = genesV2[!rm,]
            colnames(genesV2) = c("gene","NAME")
            colnames(GSEA_expr)[1] = "gene"
            GSEA_expr <- merge(genesV2,GSEA_expr,by = "gene")
            GSEA_expr = GSEA_expr[,-1]
        } else if(Mouse2Human == "toupper") GSEA_expr$NAME %<>% toupper()
        
    }
    if(is.null(file.name)) file.name = paste(unique(Idents(object)), collapse = "_")
    if(!is.null(continuous.label)) file.name <- paste(continuous.label,collapse = "_")
    if(is.null(save.path)) save.path <- paste0("output/",gsub("-","",Sys.Date()))
    if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
    write.table(GSEA_expr, file = paste0(save.path, file.name,
                                         "_",k,".txt"),
                sep = "\t", quote = FALSE,row.names = FALSE)
    
    fn = paste0(save.path, file.name,"_",k,".cls")
    if (file.exists(fn)) file.remove(fn)
    lapply(cls_list, cat, "\n", file=fn,append=TRUE)
    
    if(do.return) return(GSEA_expr)
    
}