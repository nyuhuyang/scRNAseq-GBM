add_to_seurat.1 <- function (seurat_obj = NULL, infercnv_output_path, top_n = 10,
                             bp_tolerance = 2e+06)
{
    assay = DefaultAssay(seurat_obj)
    lfiles <- list.files(infercnv_output_path, full.names = FALSE)
    if (!file.exists(paste(infercnv_output_path, "run.final.infercnv_obj",
                           sep = .Platform$file.sep))) {
        flog.warn(sprintf("::Could not find \"run.final.infercnv_obj\" file at: %s",
                          paste(infercnv_output_path, "run.final.infercnv_obj",
                                sep = .Platform$file.sep)))
        stop()
    }
    infercnv_obj = readRDS(paste(infercnv_output_path, "run.final.infercnv_obj",
                                 sep = .Platform$file.sep))
    if (is.null(seurat_obj)) {
        flog.info("No Seurat object provided, will only write metadata matrix.")
    } else if (!(setequal(row.names(seurat_obj@meta.data), colnames(infercnv_obj@expr.data)) ||
               setequal(colnames(seurat_obj[[assay]]), colnames(infercnv_obj@expr.data)))) {
        flog.warn("::Cell names in Seurat object and infercnv results do not match")
        stop()
    }
    analysis_mode_pattern = "rand_trees.hmm_mode-subclusters"
    if (!is.null(infercnv_obj@options$analysis_mode)) {
        analysis_mode_pattern = gsub("[\"]", "", infercnv_obj@options$analysis_mode)
    }
    if (any(grep(lfiles, pattern = paste0("HMM_CNV_predictions.HMM.*",
                                          analysis_mode_pattern, ".Pnorm_0.[0-9]+")))) {
        scaling_factor = 1
        if (any(grep(lfiles, pattern = paste0("HMM_CNV_predictions.HMMi6.*",
                                              analysis_mode_pattern, ".Pnorm_0.[0-9]+")))) {
            center_state = 1
        }
        else if (any(grep(lfiles, pattern = paste0("HMM_CNV_predictions.HMMi3.*",
                                                   analysis_mode_pattern, ".Pnorm_0.[0-9]+")))) {
            center_state = 1
        }
        else {
            flog.warn("::Found filtered HMM predictions output, but they do not match any known model type.")
            stop()
        }
        regions = read.table(paste(infercnv_output_path, sort(lfiles[grep(lfiles,
                                                                          pattern = paste0("HMM_CNV_predictions.HMMi[36].*",
                                                                                           analysis_mode_pattern, ".Pnorm_0.[0-9]+.pred_cnv_regions.dat"))])[1],
                                   sep = .Platform$file.sep), sep = "\t", header = TRUE,
                             check.names = FALSE)
        hmm_genes = read.table(paste(infercnv_output_path, sort(lfiles[grep(lfiles,
                                                                            pattern = paste0("HMM_CNV_predictions.HMMi[36].*",
                                                                                             analysis_mode_pattern, ".Pnorm_0.[0-9]+.pred_cnv_genes.dat"))])[1],
                                     sep = .Platform$file.sep), sep = "\t", header = TRUE,
                               check.names = FALSE)
    } else if (any(grep(lfiles, pattern = paste0("17_HMM_predHMM.*",
                                               analysis_mode_pattern)))) {
        scaling_factor = 2
        if (any(grep(lfiles, pattern = paste0("17_HMM_predHMMi6.*",
                                              analysis_mode_pattern)))) {
            center_state = 3
        }
        else if (any(grep(lfiles, pattern = paste0("17_HMM_predHMMi3.*",
                                                   analysis_mode_pattern)))) {
            center_state = 2
        }
        else {
            flog.warn("::Found HMM predictions output, but they do not match any known model type")
            stop()
        }
        regions = read.table(paste(infercnv_output_path, lfiles[grep(lfiles,
                                                                     pattern = paste0("17_HMM_predHMMi[36].*", analysis_mode_pattern,
                                                                                      ".pred_cnv_regions.dat"))][1], sep = .Platform$file.sep),
                             sep = "\t", header = TRUE, check.names = FALSE)
        hmm_genes = read.table(paste(infercnv_output_path, lfiles[grep(lfiles,
                                                                       pattern = paste0("17_HMM_predHMMi[36].*", analysis_mode_pattern,
                                                                                        ".pred_cnv_genes.dat"))][1], sep = .Platform$file.sep),
                               sep = "\t", header = TRUE, check.names = FALSE)
    } else {
        flog.warn(sprintf("::Could not find any HMM predictions outputs at: %s",
                          infercnv_output_path))
        stop()
    }
    features_to_add <- infercnv:::.get_features(infercnv_obj = infercnv_obj,
                                     infercnv_output_path = infercnv_output_path, regions = regions,
                                     hmm_genes = hmm_genes, center_state = center_state,
                                     scaling_factor = scaling_factor, top_n = top_n, bp_tolerance = bp_tolerance)
    if (!is.null(seurat_obj)) {
        for (lv in levels(infercnv_obj@gene_order$chr)) {
            seurat_obj@meta.data[[paste0("has_cnv_", lv)]] = features_to_add$feature_vector_chrs_has_cnv[[lv]]
            seurat_obj@meta.data[[paste0("has_loss_", lv)]] = features_to_add$feature_vector_chrs_has_loss[[lv]]
            seurat_obj@meta.data[[paste0("has_dupli_", lv)]] = features_to_add$feature_vector_chrs_has_dupli[[lv]]
            seurat_obj@meta.data[[paste0("proportion_cnv_",
                                         lv)]] = features_to_add$feature_vector_chrs_gene_cnv_proportion[[lv]]
            seurat_obj@meta.data[[paste0("proportion_loss_",
                                         lv)]] = features_to_add$feature_vector_chrs_gene_loss_proportion[[lv]]
            seurat_obj@meta.data[[paste0("proportion_dupli_",
                                         lv)]] = features_to_add$feature_vector_chrs_gene_dupli_proportion[[lv]]
            seurat_obj@meta.data[[paste0("proportion_scaled_cnv_",
                                         lv)]] = features_to_add$feature_vector_chrs_gene_cnv_proportion_scaled[[lv]]
            seurat_obj@meta.data[[paste0("proportion_scaled_loss_",
                                         lv)]] = features_to_add$feature_vector_chrs_gene_loss_proportion_scaled[[lv]]
            seurat_obj@meta.data[[paste0("proportion_scaled_dupli_",
                                         lv)]] = features_to_add$feature_vector_chrs_gene_dupli_proportion_scaled[[lv]]
        }
        for (n in names(features_to_add)[grep(names(features_to_add),
                                              pattern = "top_")]) {
            seurat_obj@meta.data[[n]] = features_to_add[[n]]
        }
    }
    out_mat = matrix(NA, ncol = ((9 * length(levels(infercnv_obj@gene_order$chr))) +
                                     length(features_to_add) - 9), nrow = ncol(infercnv_obj@expr.data))
    out_mat_feature_names = vector("character", ((9 * length(levels(infercnv_obj@gene_order$chr))) +
                                                     length(features_to_add) - 9))
    i = 1
    for (lv in levels(infercnv_obj@gene_order$chr)) {
        out_mat[, i] = features_to_add$feature_vector_chrs_has_cnv[[lv]]
        out_mat[, i + 1] = features_to_add$feature_vector_chrs_has_loss[[lv]]
        out_mat[, i + 2] = features_to_add$feature_vector_chrs_has_dupli[[lv]]
        out_mat[, i + 3] = features_to_add$feature_vector_chrs_gene_cnv_proportion[[lv]]
        out_mat[, i + 4] = features_to_add$feature_vector_chrs_gene_loss_proportion[[lv]]
        out_mat[, i + 5] = features_to_add$feature_vector_chrs_gene_dupli_proportion[[lv]]
        out_mat[, i + 6] = features_to_add$feature_vector_chrs_gene_cnv_proportion_scaled[[lv]]
        out_mat[, i + 7] = features_to_add$feature_vector_chrs_gene_loss_proportion_scaled[[lv]]
        out_mat[, i + 8] = features_to_add$feature_vector_chrs_gene_dupli_proportion_scaled[[lv]]
        out_mat_feature_names[i:(i + 8)] = c(paste0("has_cnv_",
                                                    lv), paste0("has_loss_", lv), paste0("has_dupli_",
                                                                                         lv), paste0("proportion_cnv_", lv), paste0("proportion_loss_",
                                                                                                                                    lv), paste0("proportion_dupli_", lv), paste0("proportion_scaled_cnv_",
                                                                                                                                                                                 lv), paste0("proportion_scaled_loss_", lv), paste0("proportion_scaled_dupli_",
                                                                                                                                                                                                                                    lv))
        i = i + 9
    }
    for (n in names(features_to_add)[grep(names(features_to_add),
                                          pattern = "top_")]) {
        out_mat[, i] = features_to_add[[n]]
        out_mat_feature_names[i] = n
        i = i + 1
    }
    colnames(out_mat) = out_mat_feature_names
    row.names(out_mat) = colnames(infercnv_obj@expr.data)
    write.table(out_mat, paste(infercnv_output_path, "map_metadata_from_infercnv.txt",
                               sep = .Platform$file.sep), quote = FALSE, sep = "\t")
    return(seurat_obj)
}




#remove NA columns and NA rows, remove duplicate Gene_name
CleanUp <- function(df){

        rm_NA_col <- df[which(df[,1] == "Gene_name"),] %>% is.na %>% as.vector
        df = df[,!rm_NA_col]
        rm_NA_row <- apply(df,1, function(x) all(is.na(x)))
        df = df[!rm_NA_row,]
        colnames(df) = df[which(df[,1] == "Gene_name"),] %>% as.character
        df = df[-which(df[,1] == "Gene_name"),]

        rm_col <- colnames(df) %in% c("Gene_id","biotype")
        df = df[,!rm_col]
        df = RemoveDup(df)

        return(df)
}

# modify doubletFinder_v3 to commendate SCT




doubletFinder_v3 <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE, annotations = NULL) {
    require(Seurat); require(fields); require(KernSmooth)
    ## Generate new list of doublet classificatons from existing pANN vector to save time
    if (reuse.pANN != FALSE ) {
        pANN.old <- seu@meta.data[ , reuse.pANN]
        classifications <- rep("Singlet", length(pANN.old))
        classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
        seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
        return(seu)
    }

    assay = DefaultAssay(seu)
    if (reuse.pANN == FALSE) {
        ## Make merged real-artifical data
        real.cells <- rownames(seu@meta.data)
        data <- seu[[assay]]@counts[, real.cells]
        n_real.cells <- length(real.cells)
        n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
        print(paste("Creating", n_doublets, "artificial doublets...", sep = " "))
        real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
        real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
        doublets <- (data[, real.cells1] + data[, real.cells2])/2
        colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
        data_wdoublets <- cbind(data, doublets)
        # Keep track of the types of the simulated doublets
        if(!is.null(annotations)){
            stopifnot(typeof(annotations)=="character")
            stopifnot(length(annotations)==length(Cells(seu)))
            stopifnot(!any(is.na(annotations)))
            annotations <- factor(annotations)
            names(annotations) <- Cells(seu)
            doublet_types1 <- annotations[real.cells1]
            doublet_types2 <- annotations[real.cells2]
        }
        ## Store important pre-processing information
        orig.commands <- seu@commands
        ## Pre-process Seurat object

        if (sct == FALSE) {
            print("Creating Seurat object...")
            seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

            print("Normalizing Seurat object...")
            seu_wdoublets <- NormalizeData(seu_wdoublets,
                                           normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                           scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                           margin = orig.commands$NormalizeData.RNA@params$margin)

            print("Finding variable genes...")
            seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
                                                  selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                                                  loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                                                  clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                                                  mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                                                  dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                                                  num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                                                  binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                                                  nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                                                  mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                                                  dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)

            print("Scaling data...")
            seu_wdoublets <- ScaleData(seu_wdoublets,
                                       features = orig.commands$ScaleData.RNA$features,
                                       model.use = orig.commands$ScaleData.RNA$model.use,
                                       do.scale = orig.commands$ScaleData.RNA$do.scale,
                                       do.center = orig.commands$ScaleData.RNA$do.center,
                                       scale.max = orig.commands$ScaleData.RNA$scale.max,
                                       block.size = orig.commands$ScaleData.RNA$block.size,
                                       min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)

            print("Running PCA...")
            seu_wdoublets <- RunPCA(seu_wdoublets,
                                    features = orig.commands$ScaleData.RNA$features,
                                    npcs = length(PCs),
                                    rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                                    weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                                    verbose=FALSE)
            pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
            cell.names <- rownames(seu_wdoublets@meta.data)
            nCells <- length(cell.names)
            rm(seu_wdoublets); gc() # Free up memory
        }
        if (sct == TRUE) {
            if (sct == TRUE) {
                require(sctransform)
                print("Creating Seurat object...")
                seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

                print("Running SCTransform...")
                seu_wdoublets <- SCTransform(seu_wdoublets)

                print("Running PCA...")
                seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
                pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
                cell.names <- rownames(seu_wdoublets@meta.data)
                nCells <- length(cell.names)
                rm(seu_wdoublets); gc()
            }
        }
        ## Compute PC distance matrix
        print("Calculating PC distance matrix...")
        dist.mat <- fields::rdist(pca.coord)

        ## Compute pANN
        print("Computing pANN...")
        pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
        if(!is.null(annotations)){
            neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = length(levels(doublet_types1))))
        }
        rownames(pANN) <- real.cells
        colnames(pANN) <- "pANN"
        k <- round(nCells * pK)
        for (i in 1:n_real.cells) {
            neighbors <- order(dist.mat[, i])
            neighbors <- neighbors[2:(k + 1)]
            neighbor.names <- rownames(dist.mat)[neighbors]
            pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
        }
        print("Classifying doublets..")
        classifications <- rep("Singlet", n_real.cells)
        classifications[order(pANN$pANN[1:n_real.cells], decreasing = TRUE)[1:nExp]] <- "Doublet"
        seu@meta.data[, paste("pANN", pN, pK, nExp, sep = "_")] <- pANN[rownames(seu@meta.data),
                                                                        1]
        seu@meta.data[, paste("DF.classifications", pN, pK,
                              nExp, sep = "_")] <- classifications
        return(seu)
    }
}

# DoubletFinder
# find histgram local maximam
find.localMaxima <- function(x) {
    # Use -Inf instead if x is numeric (non-integer)
    y <- diff(c(-.Machine$integer.max, x)) > 0L
    rle(y)$lengths
    y <- cumsum(rle(y)$lengths)
    y <- y[seq.int(1L, length(y), 2L)]
    if (x[[1]] == x[[2]]) {
        y <- y[-1]
    }
    which(x == max(x[y]))
}

# re-calculate the cdr3 frequency in data.frame
Frequency <- function(df, key = "cdr3",top=NULL,remove.na =T,
                      remove.dup =T, verbose = FALSE){

        #colnames(df) = tolower(colnames(df))
        if(any(colnames(df) %in% "frequency")){
                df = df[,-which(colnames(df) %in% "frequency")]
        }
        cdr3_table <- table(df[,key]) %>% as.data.frame
        if(verbose) table(cdr3_table$Freq) %>% print
        Sum <- sum(cdr3_table$Freq)
        if(colnames(df)[1] == "ngene"|colnames(df)[1] == "nGene"){
        cdr3_table$frequency = cdr3_table$Freq/Sum
        }
        cdr3_table = cdr3_table[order(cdr3_table$Freq,decreasing = T),]
        if(remove.na) cdr3_table = cdr3_table[cdr3_table$Var1!="na",]
        colnames(cdr3_table)[1] = "cdr3"
        df_new <- left_join(df, cdr3_table,by = "cdr3")
        if(remove.dup) df_new = df_new[!duplicated(df_new$cdr3),]
        df_new$Freq[is.na(df_new$Freq)]=0
        df_new$frequency[is.na(df_new$frequency)]=0
        df_new = df_new[order(df_new$freq,decreasing = T),]
        if(is.null(top)) top = nrow(df_new)
        top = min(top,nrow(df_new))
        return(df_new[1:top,])
}

# DouletFinder
Multiplet_Rate <- function(object, numBatches = 1, num10xRuns = 1){

    numCellsRecovered = 1.0 * ncol(object)
    m = 4.597701e-06
    r = 0.5714286

    numCellsLoaded = numCellsRecovered / r
    multipletRate = m * numCellsLoaded / num10xRuns

    singletRate = 1.0 - multipletRate;
    numSinglet = singletRate * numCellsRecovered
    numMultiplet = numCellsRecovered - numSinglet
    numIdentMultiplet = numMultiplet * (numBatches - 1) / numBatches
    numNonIdentMultiplet = numMultiplet - numIdentMultiplet
    numCells = numSinglet + numNonIdentMultiplet

    return(numNonIdentMultiplet/numCells)
}

# doubletFinder
# http://rstudio-pubs-static.s3.amazonaws.com/329613_f53e84d1a18840d5a1df55efb90739d9.html
qplot_2axis <- function(data,x = "pK", y1 = "MeanBC", y2 = "BCmetric"){
    if(class(data[,x]) == "factor") data[,x] <- as.numeric(as.character(data[,x]))
    data_y1 <- data[,y1]
    data_y2 <- data[,y2]
    a <- range(data_y1)
    b <- range(data_y2)
    scale_factor <- diff(a)/diff(b)
    data_y2 <- ((data_y2 - b[1]) * scale_factor) + a[1]
    trans <- ~ ((. - a[1]) / scale_factor) + b[1]

    g <- ggplot(data = data, aes_string(x = x, y = y1))+
        geom_line()+geom_point()+
        geom_point(aes(y = data_y2),colour = "blue")+
        geom_line(aes(y = data_y2),colour = "blue")+
        scale_y_continuous(name = y1,
                           sec.axis = sec_axis(trans=trans, name=y2))+
        theme(axis.text.y.right = element_text(color = "blue"))

    g

}

#remove duplicate rownames with lower rowsumns
#' @param mat input as data.frame with gene name
#' @export mat matrix with gene as rownames, no duplicated genes
RemoveDup <- function(mat){
        gene_id <- as.matrix(mat[,1])
        mat <- mat[,-1]
        if(!is.matrix(mat)) mat <- sapply(mat,as.numeric)
        rownames(mat) <- 1:nrow(mat)
        mat[is.na(mat)] = 0
        mat <- cbind(mat, "rowSums" = rowSums(mat))
        mat <- mat[order(mat[,"rowSums"],decreasing = T),]
        gene_id <- gene_id[as.numeric(rownames(mat))]
        remove_index <- duplicated(gene_id)
        mat <- mat[!remove_index,]
        rownames(mat) <- gene_id[!remove_index]
        return(mat[,-ncol(mat)])
}


#' select 1/4 of cell from control
ScaleDown <- function(object, control=c("BH","DJ","MD","NZ")){

        normal_cells = lapply(control, function(x){
        rownames(object@meta.data)[(object@meta.data$orig.ident %in% x)]
        })
        set.seed(101)
        remove_normal_cells = lapply(normal_cells, function(x) sample(x, size = length(x)*3/4)) %>% unlist
        table(object@cell.names %in% remove_normal_cells)
        cell.use <- object@cell.names[!(object@cell.names %in% remove_normal_cells)]
        object <- SubsetData(object, cells.use = cell.use)
        object@meta.data$orig.ident = gsub(paste(control,collapse = "|"),"Normal",object@meta.data$orig.ident)

        return(object)
}

#' select 10% of cell from control
#' Seurat 3
ScaleDown.1 <- function(object, control=c("BH","DJ","MD","NZ"), scale.down=0.1){

    normal_cells = lapply(control, function(x){
        rownames(object@meta.data)[(object@meta.data$orig.ident %in% x)]
    })
    set.seed(101)
    remove_normal_cells = lapply(normal_cells, function(x) {
        sample(x, size = length(x)*(1-scale.down))
        }) %>% unlist
    table(colnames(object) %in% remove_normal_cells)
    cell.use <- colnames(object)[!(colnames(object) %in% remove_normal_cells)]
    object <- object[, cell.use]
    object@meta.data$orig.ident = gsub(paste(control,collapse = "|"),"Normal",object@meta.data$orig.ident)

    return(object)
}


#' Generate Trends of top 20 TCR clones in longitudinal samples
#' @param TCR_data TCR data.frame report. must contain columns like c("orig.ident","cdr3","frequency")
#' @param order.by character vector indicates sample longitudinal order
#' @param group character, the summary of all samples. For output figure titile only
#' @param size ggplot point size
#' @param top integer, include top N TCR clones from each time points
#' @param x.margin integer, blank margin on x axis
#' @example order.by = c("Pt17_C03_TCRB","Pt17_C07_TCRB","Pt17_C31_TCRB")
#' TCR_Longitud_Plot(Pt_17_TCR,order.by, "Pt-17 bulk")
TCR_Longitud_Plot <- function(TCR_data, order.by, group,color=NULL,
size=2,colors = singler.colors,x.margin=0.125,
do.return = TRUE, do.print = FALSE,do.log =TRUE, top=20){

    colnames(TCR_data) = tolower(colnames(TCR_data))
    TCR_data$orig.ident = as.character(TCR_data$orig.ident)
    TCR_data <- split(TCR_data, TCR_data$orig.ident) %>%
    lapply(function(x) Frequency(x,top=top)) %>% bind_rows
    keep_col <- c("orig.ident","cdr3","frequency")
    TCR_data = TCR_data[order(TCR_data$frequency,decreasing = T),]
    TCR_data <- TCR_data[,keep_col]
    # fill up 0
    TCR_data %<>% spread(key="orig.ident", value="frequency",fill = 0)
    TCR_data %<>% gather(key="orig.ident",value="frequency",-cdr3)

    TCR_data$orig.ident = factor(TCR_data$orig.ident,levels = order.by)
    TCR_data$samples <- as.integer(TCR_data$orig.ident)
    g1 <- ggplot(data=TCR_data, aes(x=samples, y=frequency, group=cdr3,colour=cdr3)) +
    geom_line(size = size) +
    geom_point(size = size)+
    #scale_fill_manual(values = colors)+
    theme_minimal()+
    xlab("longitudinal samples")+
    scale_x_continuous(limits=c(min(TCR_data$samples)-x.margin,
    max(TCR_data$samples)+x.margin),
    breaks=(min(TCR_data$samples)):(max(TCR_data$samples)),
    labels=order.by)+
    scale_y_continuous(labels = scales::percent)+
    ggtitle(paste("Trends of top",top,"TCR clones in longitudinal",group,
    "samples"))+
    theme(text = element_text(size=15),
    legend.position="none",
    plot.title = element_text(hjust = 0.5),
    axis.text.x  = element_text(angle=70, vjust=0.5))
    if(do.log) g1 = g1 + scale_y_log10()
    if(do.print){
        path <- paste0("output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,group,"_TCR.jpeg"), units="in", width=10, height=7,res=600)
        print(g1)
        dev.off()
    }
    if(do.return) return(g1)
}


#' Generate Trends of top 20 TCR clones in longitudinal samples
#' @param TCR_data TCR data.frame report. must contain columns like c("orig.ident","cdr3","frequency")
#' @param order.by length =2character vector indicates pariwise sample at x-axis and y-axis
#' @param top integer, include top N TCR clones from each time points
#' @param size ggplot point size
#' @param margin integer, blank margin on both x and y axis
#' @example TCRPairPlot(Pt_17_TCR,c("Pt17_C03_TCRB","Pt17_C07_TCRB"))
TCRPairPlot <- function(TCR_data, order.by, size=2,margin=0.1,
do.return = TRUE, do.print = FALSE,top=20){
    colnames(TCR_data) = tolower(colnames(TCR_data))
    TCR_data = TCR_data[(TCR_data$orig.ident %in% order.by),]
    TCR_data <- split(TCR_data, TCR_data$orig.ident) %>%
    lapply(function(x) Frequency(x,top = top)) %>% bind_rows
    keep_col <- c("orig.ident","cdr3","frequency")
    TCR_data = TCR_data[order(TCR_data$frequency,decreasing = T),]
    TCR_data <- TCR_data[,keep_col]
    # fill up 0
    TCR_data %<>% spread(key="orig.ident", value="frequency",fill = 0)
    TCR_data = TCR_data[,1:3]
    colnames(TCR_data)[2:3] =c("sample1","sample2")
    TCR_data$counts = as.numeric(TCR_data$sample1>0)+
    as.numeric(TCR_data$sample2>0)

    g1 <- ggplot(data=TCR_data, aes(x=sample1, y=sample2,
    group=counts, color=counts)) +
    geom_point(size = 2)+
    scale_fill_manual(values = singler.colors)+
    theme_minimal()+
    xlab(paste("frequency in",order.by[1]))+
    ylab(paste("frequency in",order.by[2]))+
    scale_x_sqrt(labels = scales::percent,expand=c(0,0),
    limits=c(0,max(TCR_data$sample1)*(1+margin)))+
    scale_y_sqrt(labels = scales::percent,expand=c(0,0),
    limits=c(0,max(TCR_data$sample2)*(1+margin)))+
    ggtitle(paste("Top",top,"TCR clones between",
    paste(order.by,collapse = " and ")))+
    theme(text = element_text(size=15),
    legend.position="none",
    plot.title = element_text(hjust = 0.5),
    axis.text.x  = element_text(angle=70, vjust=0.5))

    if(do.print){
        path <- paste0("output/",gsub("-","",Sys.Date()),"/")
        if(!dir.exists(path)) dir.create(path, recursive = T)
        jpeg(paste0(path,paste0(order.by,collapse = "_"),"_pairTCR.jpeg"), units="in", width=10, height=7,res=600)
        print(g1)
        dev.off()
    }
    if(do.return) return(g1)
}


show.velocity.on.embedding.cor.1 <- function (emb, vel, n = 100, cell.colors = NULL, arrow.colors = par("fg"),
                                              corr.sigma = 0.05,
          show.grid.flow = FALSE, grid.n = 20, grid.sd = NULL, min.grid.cell.mass = 1,
          min.arrow.size = NULL, arrow.scale = 1, max.grid.arrow.length = NULL,
          fixed.arrow.length = FALSE, plot.grid.points = FALSE, scale = "log",
          nPcs = NA, arrow.lwd = 1, xlab = "", ylab = "", n.cores = velocyto.R:::defaultNCores(),
          do.par = T, show.cell = NULL, cell.border.alpha = 0.3, cc = NULL,
          return.details = FALSE, expression.scaling = FALSE, ...)
{
    randomize <- FALSE
    if (do.par)
        par(mfrow = c(1, 1), mar = c(3.5, 3.5, 2.5, 1.5), mgp = c(2,
                                                                  0.65, 0), cex = 0.85)
    celcol <- "white"
    if (is.null(show.cell)) {
        celcol <- cell.colors[rownames(emb)]
    }
    plot(emb, bg = celcol, pch = 21, col = ac(1, alpha = cell.border.alpha),
         xlab = xlab, ylab = ylab,fg = gray(0.7), ...)
    em <- as.matrix(vel$current)
    ccells <- intersect(rownames(emb), colnames(em))
    em <- em[, ccells]
    emb <- emb[ccells, ]
    nd <- as.matrix(vel$deltaE[, ccells])
    cgenes <- intersect(rownames(em), rownames(nd))
    nd <- nd[cgenes, ]
    em <- em[cgenes, ]
    if (randomize) {
        nd <- t(apply(nd, 1, function(x) (rbinom(length(x),
                                                 1, 0.5) * 2 - 1) * abs(sample(x))))
    }
    if (is.null(cc)) {
        cat("delta projections ... ")
        if (scale == "log") {
            cat("log ")
            cc <- velocyto.R:::colDeltaCorLog10(em, (log10(abs(nd) + 1) *
                                            sign(nd)), nthreads = n.cores)
        }
        else if (scale == "sqrt") {
            cat("sqrt ")
            cc <- velocyto.R:::colDeltaCorSqrt(em, (sqrt(abs(nd)) * sign(nd)),
                                  nthreads = n.cores)
        }
        else if (scale == "rank") {
            cat("rank ")
            cc <- velocyto.R:::colDeltaCor((apply(em, 2, rank)), (apply(nd,
                                                           2, rank)), nthreads = n.cores)
        }
        else {
            cat("linear ")
            cc <- velocyto.R:::colDeltaCor(em, nd, nthreads = n.cores)
        }
        colnames(cc) <- rownames(cc) <- colnames(em)
        diag(cc) <- 0
    }
    cat("knn ... ")
    if (n > nrow(cc)) {
        n <- nrow(cc)
    }
    emb.knn <- velocyto.R:::balancedKNN(t(emb), k = n, maxl = nrow(emb),
                           dist = "euclidean", n.threads = n.cores)
    diag(emb.knn) <- 1
    cat("transition probs ... ")
    tp <- exp(cc/corr.sigma) * emb.knn
    tp <- t(t(tp)/Matrix::colSums(tp))
    tp <- as(tp, "dgCMatrix")
    cat("done\n")
    if (!is.null(show.cell)) {
        i <- match(show.cell, rownames(emb))
        if (is.na(i))
            stop(paste("specified cell", i, "is not in the embedding"))
        points(emb, pch = 19, col = ac(val2col(tp[rownames(emb),
                                                  show.cell], gradient.range.quantile = 1), alpha = 0.5))
        points(emb[show.cell, 1], emb[show.cell, 2], pch = 3,
               cex = 1, col = 1)
        di <- t(t(emb) - emb[i, ])
        di <- di/sqrt(Matrix::rowSums(di^2)) * arrow.scale
        di[i, ] <- 0
        dir <- Matrix::colSums(di * tp[, i])
        dic <- Matrix::colSums(di * (tp[, i] > 0)/sum(tp[, i] >
                                                          0))
        dia <- dir - dic
        suppressWarnings(arrows(emb[colnames(em)[i], 1], emb[colnames(em)[i],
                                                             2], emb[colnames(em)[i], 1] + dic[1], emb[colnames(em)[i],
                                                                                                       2] + dic[2], length = 0.05, lwd = 1, col = "blue"))
        suppressWarnings(arrows(emb[colnames(em)[i], 1], emb[colnames(em)[i],
                                                             2], emb[colnames(em)[i], 1] + dir[1], emb[colnames(em)[i],
                                                                                                       2] + dir[2], length = 0.05, lwd = 1, col = "red"))
        suppressWarnings(arrows(emb[colnames(em)[i], 1] + dic[1],
                                emb[colnames(em)[i], 2] + dic[2], emb[colnames(em)[i],
                                                                      1] + dir[1], emb[colnames(em)[i], 2] + dir[2],
                                length = 0.05, lwd = 1, lty = 1, col = "grey50"))
        suppressWarnings(arrows(emb[colnames(em)[i], 1], emb[colnames(em)[i],
                                                             2], emb[colnames(em)[i], 1] + dia[1], emb[colnames(em)[i],
                                                                                                       2] + dia[2], length = 0.05, lwd = 1, col = "black"))
    }
    else {
        cat("calculating arrows ... ")
        arsd <- data.frame(t(velocyto.R:::embArrows(emb, tp, arrow.scale,
                                       n.cores)))
        rownames(arsd) <- rownames(emb)
        if (expression.scaling) {
            tpb <- tp > 0
            tpb <- t(t(tpb)/colSums(tpb))
            es <- as.matrix(em %*% tp) - as.matrix(em %*% as.matrix(tpb))
            pl <- pmin(1, pmax(0, apply(as.matrix(vel$deltaE[,
                                                             colnames(es)]) * es, 2, sum)/sqrt(colSums(es *
                                                                                                           es))))
            arsd <- arsd * pl
        }
        ars <- data.frame(cbind(emb, emb + arsd))
        colnames(ars) <- c("x0", "y0", "x1", "y1")
        colnames(arsd) <- c("xd", "yd")
        rownames(ars) <- rownames(emb)
        cat("done\n")
        if (show.grid.flow) {
            cat("grid estimates ... ")
            rx <- range(c(range(ars$x0), range(ars$x1)))
            ry <- range(c(range(ars$y0), range(ars$y1)))
            gx <- seq(rx[1], rx[2], length.out = grid.n)
            gy <- seq(ry[1], ry[2], length.out = grid.n)
            if (is.null(grid.sd)) {
                grid.sd <- sqrt((gx[2] - gx[1])^2 + (gy[2] -
                                                         gy[1])^2)/2
                cat("grid.sd=", grid.sd, " ")
            }
            if (is.null(min.arrow.size)) {
                min.arrow.size <- sqrt((gx[2] - gx[1])^2 + (gy[2] -
                                                                gy[1])^2) * 0.01
                cat("min.arrow.size=", min.arrow.size, " ")
            }
            if (is.null(max.grid.arrow.length)) {
                max.grid.arrow.length <- sqrt(sum((par("pin")/c(length(gx),
                                                                length(gy)))^2)) * 0.25
                cat("max.grid.arrow.length=", max.grid.arrow.length,
                    " ")
            }
            garrows <- do.call(rbind, lapply(gx, function(x) {
                cd <- sqrt(outer(emb[, 2], -gy, "+")^2 + (x -
                                                              emb[, 1])^2)
                cw <- dnorm(cd, sd = grid.sd)
                gw <- Matrix::colSums(cw)
                cws <- pmax(1, Matrix::colSums(cw))
                gxd <- Matrix::colSums(cw * arsd$xd)/cws
                gyd <- Matrix::colSums(cw * arsd$yd)/cws
                al <- sqrt(gxd^2 + gyd^2)
                vg <- gw >= min.grid.cell.mass & al >= min.arrow.size
                cbind(rep(x, sum(vg)), gy[vg], x + gxd[vg],
                      gy[vg] + gyd[vg])
            }))
            colnames(garrows) <- c("x0", "y0", "x1", "y1")
            if (fixed.arrow.length) {
                suppressWarnings(arrows(garrows[, 1], garrows[,
                                                              2], garrows[, 3], garrows[, 4], length = 0.05,
                                        lwd = arrow.lwd))
            }
            else {
                alen <- pmin(max.grid.arrow.length, sqrt(((garrows[,
                                                                   3] - garrows[, 1]) * par("pin")[1]/diff(par("usr")[c(1,
                                                                                                                        2)]))^2 + ((garrows[, 4] - garrows[, 2]) *
                                                                                                                                       par("pin")[2]/diff(par("usr")[c(3, 4)]))^2))
                suppressWarnings(lapply(1:nrow(garrows), function(i) arrows(garrows[i,1], garrows[i, 2], garrows[i, 3], garrows[i,4],
                                                                            length = alen[i], lwd = arrow.lwd, col = arrow.colors)))
            }
            if (plot.grid.points)
                points(rep(gx, each = length(gy)), rep(gy, length(gx)),
                       pch = ".", cex = 0.1, col = ac(1, alpha = 0.4))
            cat("done\n")
            if (return.details) {
                cat("expression shifts .")
                scale.int <- switch(scale, log = 2, sqrt = 3,
                                    1)
                if (!expression.scaling) {
                    tpb <- tp > 0
                    tpb <- t(t(tpb)/colSums(tpb))
                    es <- as.matrix(em %*% tp) - as.matrix(em %*%
                                                               as.matrix(tpb))
                }
                cat(".")
                gs <- do.call(cbind, parallel::mclapply(gx,
                                                        function(x) {
                                                            cd <- sqrt(outer(emb[, 2], -gy, "+")^2 +
                                                                           (x - emb[, 1])^2)
                                                            cw <- dnorm(cd, sd = grid.sd)
                                                            gw <- Matrix::colSums(cw)
                                                            cws <- pmax(1, Matrix::colSums(cw))
                                                            cw <- t(t(cw)/cws)
                                                            gxd <- Matrix::colSums(cw * arsd$xd)
                                                            gyd <- Matrix::colSums(cw * arsd$yd)
                                                            al <- sqrt(gxd^2 + gyd^2)
                                                            vg <- gw >= min.grid.cell.mass & al >= min.arrow.size
                                                            if (any(vg)) {
                                                                z <- es %*% cw[, vg]
                                                            }
                                                            else {
                                                                NULL
                                                            }
                                                        }, mc.cores = n.cores, mc.preschedule = T))
                if (scale == "log") {
                    nd <- (log10(abs(nd) + 1) * sign(nd))
                }
                else if (scale == "sqrt") {
                    nd <- (sqrt(abs(nd)) * sign(nd))
                }
                cat(".")
                gv <- do.call(cbind, parallel::mclapply(gx,
                                                        function(x) {
                                                            cd <- sqrt(outer(emb[, 2], -gy, "+")^2 +
                                                                           (x - emb[, 1])^2)
                                                            cw <- dnorm(cd, sd = grid.sd)
                                                            gw <- Matrix::colSums(cw)
                                                            cws <- pmax(1, Matrix::colSums(cw))
                                                            cw <- t(t(cw)/cws)
                                                            gxd <- Matrix::colSums(cw * arsd$xd)
                                                            gyd <- Matrix::colSums(cw * arsd$yd)
                                                            al <- sqrt(gxd^2 + gyd^2)
                                                            vg <- gw >= min.grid.cell.mass & al >= min.arrow.size
                                                            if (any(vg)) {
                                                                z <- nd %*% cw[, vg]
                                                            }
                                                            else {
                                                                NULL
                                                            }
                                                        }, mc.cores = n.cores, mc.preschedule = T))
                cat(". done\n")
                return(invisible(list(tp = tp, cc = cc, garrows = garrows,
                                      arrows = as.matrix(ars), vel = nd, eshifts = es,
                                      gvel = gv, geshifts = gs, scale = scale)))
            }
        }
        else {
            apply(ars, 1, function(x) {
                if (fixed.arrow.length) {
                    suppressWarnings(arrows(x[1], x[2], x[3],
                                            x[4], length = 0.05, lwd = arrow.lwd))
                }
                else {
                    ali <- sqrt(((x[3] - x[1]) * par("pin")[1]/diff(par("usr")[c(1,
                                                                                 2)]))^2 + ((x[4] - x[2]) * par("pin")[2]/diff(par("usr")[c(3,
                                                                                                                                            4)]))^2)
                    suppressWarnings(arrows(x[1], x[2], x[3], col = arrow.colors,
                                            x[4], length = min(0.05, ali), lwd = arrow.lwd))
                }
            })
        }
    }
    return(invisible(list(tp = tp, cc = cc)))
}
