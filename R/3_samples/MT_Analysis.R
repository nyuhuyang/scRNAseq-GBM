# Libraries and Parameters ------------------------------------------------

library(Seurat)
library(ggplot2)
library(plyr)
library(dplyr)
library(cowplot)
library(ggplotify)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(SingleCellExperiment)
library(ggsci)
library(ggthemes)

# Need to increase the memory allocated to run some of the code
memory.limit(size=56000)


# Load Seurat Object -------------------------------------------

all.samples.integrated <- readRDS("GBM_PBMC_20210908.rds")

# Not doing this name change
#Make final_celltype character and then rename Unknown
all.samples.integrated$final_celltype <-  as.character(all.samples.integrated$final_celltype)
 
nameEryth <- (all.samples.integrated$final_celltype %in% "Unknown") == TRUE
all.samples.integrated$final_celltype[nameEryth] <- "Immature Erythrocytes"
# 
# nameMacrophages <- (all.samples.integrated$final_celltype %in% "Macrophage Cluster 1") == TRUE
# all.samples.integrated$final_celltype[nameMacrophages] <- "Macrophage Cluster 1"
# 
# nameAAM <- (all.samples.integrated$final_celltype %in% "Macrophage Cluster 2") == TRUE
# all.samples.integrated$final_celltype[nameAAM] <- "Macrophage Cluster 2"

#add levels to the final cell type of the Seurat object
all.samples.integrated$final_celltype <- factor(x = all.samples.integrated$final_celltype, 
                                                levels = c("Microglia",
                                                           "Activated Microglia",
                                                           "Macrophage Cluster 1",
                                                           "Macrophage Cluster 2",
                                                           "Monocytes",
                                                           "T cells",
                                                           "CD8 T cells",
                                                           "CD4 T cells",
                                                           "NK cells",
                                                           "B cells",
                                                           "Immature Erythrocytes"))


# Subsets Used in Analysis --------------------------------------------

# Tumor
subset_tumor<-subset(all.samples.integrated, subset = Samples == c('HGG','LGG'))

# Macrophage Cluster 1 and Microglia
subset <- subset(all.samples.integrated, subset = final_celltype %in% c("Macrophage Cluster 1",
                                                                        "Macrophage Cluster 2",
                                                                        "Microglia",
                                                                        "Activated Microglia",
                                                                        "Monocytes"))
# CD4 T Cells
subset_Cd4 <- subset(all.samples.integrated, subset = final_celltype == "CD4 T cells")

# CD8 T Cells
subset_Cd8<-subset(all.samples.integrated, subset = final_celltype == "CD8 T cells")


# T Cells
subset_T <- subset(all.samples.integrated, subset = final_celltype == "T cells")

#Total T Cells
subset_TotalTCells <- subset(all.samples.integrated, subset = final_celltype %in% 
                               c("T cells","CD4 T cells","CD8 T cells"))
                
# Macrophage Cluster 1
subset_Macrophages <- subset(all.samples.integrated, subset = final_celltype == "Macrophage Cluster 1")

# Macrophage Cluster 2
subset_AAM <-subset(all.samples.integrated, subset = final_celltype == "Macrophage Cluster 2")

# Microglia 
subset_Microglia <- subset(all.samples.integrated, subset = final_celltype == "Microglia")

# Activated Microglia
subset_am <- subset(all.samples.integrated, subset = final_celltype == "Activated Microglia")

# Tumor Macrophage Cluster 1
subset_Macrophages_tumor <- subset(subset_tumor, subset = final_celltype == "Macrophage Cluster 1")

# Tumor Macrophage Cluster 2
subset_AAM_tumor <- subset(subset_tumor, subset = final_celltype == "Macrophage Cluster 2")

# Tumor Microglia
subset_microglia_tumor<-subset(subset_tumor, subset = final_celltype == "Microglia")

# Tumor Activated Microglia
subset_active_microglia_tumor<-subset(subset_tumor, subset = final_celltype == "Activated Microglia")

# Tumor All T Cells
subset_totalT_tumor<-subset(subset_tumor, subset = final_celltype %in% c("T cells","CD4 T cells","CD8 T cells"))

# Tumor CD4 T Cells
subset_cd4T_tumor<-subset(subset_tumor, subset = final_celltype == "CD4 T cells")

# Tumor CD8 T Cells
subset_cd8T_tumor<-subset(subset_tumor, subset = final_celltype == "CD8 T cells")

# HGG and No Tumor Groups
subset_Hgg_no_tumor <- subset(all.samples.integrated, subset = Samples == c('HGG','No tumor'))

# LGG and No Tumor Groups
subset_Lgg_no_tumor <- subset(all.samples.integrated, subset = Samples == c('LGG','No tumor'))

# Version of Macrophage Cluster 1 & 2 that are characterized by sample identity
all.samples.integrated_sample_ident <-SetIdent(all.samples.integrated, value = "Samples")
subset_Macrophages_all<-subset(all.samples.integrated_sample_ident, subset = final_celltype == "Macrophage Cluster 1")
subset_AAM_all<-subset(all.samples.integrated_sample_ident, subset = final_celltype == "Macrophage Cluster 2")

# Macrophage subset of HGG/No Tumor Subset
subset_Macrophages_high_no <-subset(subset_Hgg_no_tumor, subset = final_celltype == "Macrophage Cluster 1")
subset_AAM_high_no <-subset(subset_Hgg_no_tumor, subset = final_celltype == "Macrophage Cluster 2")

# Macrophage subset of LGG/No Tumor Subset
subset_Macrophages_low_no <-subset(subset_Lgg_no_tumor, subset = final_celltype == "Macrophage Cluster 1")
subset_AAM_low_no <-subset(subset_Lgg_no_tumor, subset = final_celltype == "Macrophage Cluster 2")

# Macrophage Clusters with final_celltype identity
subset_all_macrophages<-subset(all.samples.integrated, 
                           subset = final_celltype == c('Macrophage Cluster 1','Macrophage Cluster 2'))
subset_all_macrophages<-SetIdent(subset_all_macrophages, value = "final_celltype")

# Macrophage Cluster 2 & Microglia
subset_macromicro<-subset(all.samples.integrated, 
                          subset = final_celltype %in% c('Microglia','Activated Microglia','Macrophage Cluster 2'))
subset_macromicro<-SetIdent(subset_macromicro, value = "final_celltype")

#Macrophage Cluster 2 Defined by Samples
subset_mac2<-subset(all.samples.integrated, subset = final_celltype == c('Macrophage Cluster 2'))
subset_mac2<-SetIdent(subset_mac2, value = "Samples")

#Macrophage Cluster 2 Defined by Samples
subset_mac1<-subset(all.samples.integrated, subset = final_celltype == c('Macrophage Cluster 1'))
subset_mac1<-SetIdent(subset_mac1, value = "Samples")



# Integrated DimPlot ------------------------------------------------------

col_palette <- c("#F8766D","#DB8E00","#00C1A7","#64B200","#00BD5C","#AEA200",
                 "#00BADE","#00A6FF","#B385FF","#EF67EB","#FF63B6" ) 

png("all_samples_integrated_revised_celltype.png", units="in", res = 1000, w = 6, h = 6)
integrated_plot<-DimPlot(all.samples.integrated, 
                         reduction = "umap", 
                         group.by = "final_celltype", 
                         label = F, repel = T,
                         cols = col_palette)  +
  ggtitle("") + 
  theme(title = element_text(face = "italic", size = 16, hjust = 0.5), 
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        legend.position = "bottom",
        legend.text = element_text(size=8)) + 
  NoAxes()
integrated_plot
dev.off()

# p <- Seurat::DimPlot(all.samples.integrated) # Generate the tSNE plot, but save it as an object
# pbuild <- ggplot2::ggplot_build(p) # Use ggplot_build to deconstruct the ggplot object
# pdata <- pbuild$data[[1]] # Pull the data used for the plot
# pdata <-  pdata[order(pdata$group), ] # Order the plot data by group
# ucols <- unique(pdata$colour) # Get a vector of unique colors
# names(ucols) <- unique(pdata$group)

# DimPlot Split by Samples ----------------------------------------------

png("all_integrated_revised_celltype_splitby_samples.png", units="in", res = 1000, w = 15, h = 6)

DimPlot(all.samples.integrated, reduction = "umap", 
        group.by = "final_celltype", 
        label = F, 
        repel = T, 
        split.by = "Samples",
        cols = col_palette)  + 
  ggtitle("") + 
  theme(title = element_text(face = "italic", size = 16, hjust = 0.5), 
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        legend.position = "bottom") + 
  NoAxes()
dev.off()

# Features for Plots ------------------------------------------------------

features_for_plots<-c("Gpr34", "Tmem119", "Il1b", "Ccl3", "Adgre1",'Cd68', 
                      "Ly6c1", "Cd3d", "Cd3e","Cd4", "Tnfrsf4","Cd8a", "Cd8b1", 
                      "Ncr1","Cd79a", "Cd19")
# Flipped
features_for_plots_f<-c("Cd19","Cd79a", "Ncr1", "Cd8b1","Cd8a", "Tnfrsf4","Cd4", 
                      "Cd3e", "Cd3d", "Ly6c1",'Cd68', "Adgre1", "Ccl3", "Il1b",
                      "Tmem119","Gpr34")


# DotPlot of Selected Features --------------------------------------------

all.samples.integrated<-SetIdent(all.samples.integrated, value = "final_celltype")   #set ident to "final_celltype"
png("marker_gene_expression_all_samples.png", units="in", res = 1000, w = 8.5, h = 5)
dotplot<-DotPlot(all.samples.integrated, 
                 features = features_for_plots, 
                 col.max = 5, 
                 col.min =0 , 
                 cols = c("Blue", "Red"),
                 idents = c("Microglia",
                            "Activated Microglia",
                            "Macrophage Cluster 1",
                            "Macrophage Cluster 2",
                            "Monocytes",
                            "T cells",
                            "CD8 T cells",
                            "CD4 T cells",
                            "NK cells",
                            "B cells"), 
                 dot.scale = 4) + 
  RotatedAxis() + 
  xlab("") + 
  ylab("")
print(dotplot)
dev.off()

png("marker_gene_expression_all_samples_flipped_axis.png", units="in", res = 1000, w = 8.5, h = 5)
dotplot<-DotPlot(all.samples.integrated, 
                 features = features_for_plots_f, 
                 col.max = 5, 
                 col.min =0 , 
                 cols = c("Blue", "Red"),
                 idents = c("Microglia",
                            "Activated Microglia",
                            "Macrophage Cluster 1",
                            "Macrophage Cluster 2",
                            "Monocytes",
                            "T cells",
                            "CD8 T cells",
                            "CD4 T cells",
                            "NK cells",
                            "B cells"), 
                 dot.scale = 4) + 
  RotatedAxis() + 
  xlab("") + 
  ylab("")
print(dotplot)
dev.off()
# ViolinPlot of Selected Features -----------------------------------------

png("marker_genes_violin_plot.png",units="in",res = 1000, w =5, h = 15)
violin_plot<-VlnPlot(all.samples.integrated, group.by = "final_celltype",
                     idents = c("Microglia","Activated Microglia",
                                "Macrophage Cluster 1","Macrophage Cluster 2",
                                "Monocytes","T cells","CD4 T cells",
                                "CD8 T cells","NK cells","B cells"),
                     stack = TRUE, 
                     features = features_for_plots, 
                     flip = TRUE, 
                     pt.size = 1, 
                     fill.by = "ident", 
                     sort = FALSE) + 
  NoLegend() + 
  theme(axis.text.x = element_text(angle = 90), 
        axis.text.y = element_text(angle = 90)) + 
  xlab("")
violin_plot
dev.off()

png("marker_genes_violin_plot_flipped.png",units="in",res = 1000, w =5, h = 15)
violin_plot<-VlnPlot(all.samples.integrated, group.by = "final_celltype",
                     idents = c("Microglia","Activated Microglia",
                                "Macrophage Cluster 1","Macrophage Cluster 2",
                                "Monocytes","T cells","CD4 T cells",
                                "CD8 T cells","NK cells","B cells"),
                     stack = TRUE, 
                     features = features_for_plots_f, 
                     flip = TRUE, 
                     pt.size = 1, 
                     fill.by = "ident", 
                     sort = FALSE) + 
  NoLegend() + 
  theme(axis.text.x = element_text(angle = 90), 
        axis.text.y = element_text(angle = 90)) + 
  xlab("")
violin_plot
dev.off()

# Integrated Feature Plots ------------------------------------------------

png("feature_plot_cellwise.png", units="in", res = 2500, w = 8.5, h =5.2)   ## prepare feature plot for individual gene and then merge it using cowplot
plot1<-FeaturePlot(all.samples.integrated, reduction = "umap", features = "Cd19",coord.fixed = TRUE, cols = c("lightgrey", "Brown")) + NoAxes() + theme(plot.title = element_text(face = "italic", size = 8, color = "chocolate1"), legend.position = "bottom", legend.text = element_text(size = 5))
plot2<-FeaturePlot(all.samples.integrated, reduction = "umap", features = "Cd79a",coord.fixed = TRUE, cols = c("lightgrey", "Brown")) + NoAxes() + theme(plot.title = element_text(face = "italic", size = 8, color = "chocolate1"))
plot3<-FeaturePlot(all.samples.integrated, reduction = "umap", features = "Cd3e",coord.fixed = TRUE, cols = c("lightgrey", "Brown")) + NoAxes() + theme(plot.title = element_text(face = "italic", size = 8, color = "chocolate1"))
plot9<-FeaturePlot(all.samples.integrated, reduction = "umap", features = "Cd3d",coord.fixed = TRUE, cols = c("lightgrey", "Brown")) + NoAxes() + theme(plot.title = element_text(face = "italic", size = 8, color = "chocolate1"))
plot4<-FeaturePlot(all.samples.integrated, reduction = "umap", features = "Cd4",coord.fixed = TRUE, cols = c("lightgrey", "Brown")) + NoAxes() + theme(plot.title = element_text(face = "italic", size = 8, color = "chocolate1"))
plot5<-FeaturePlot(all.samples.integrated, reduction = "umap", features = "Ncr1",coord.fixed = TRUE, cols = c("lightgrey", "Brown")) + NoAxes() + theme(plot.title = element_text(face = "italic", size = 8, color = "chocolate1"))
plot6<-FeaturePlot(all.samples.integrated, reduction = "umap", features = "Tmem119",coord.fixed = TRUE, cols = c("lightgrey", "Brown")) + NoAxes() + theme(plot.title = element_text(face = "italic", size = 8, color = "chocolate1"))
plot7<-FeaturePlot(all.samples.integrated, reduction = "umap", features = "P2ry12",coord.fixed = TRUE, cols = c("lightgrey", "Brown")) + NoAxes() + theme(plot.title = element_text(face = "italic", size = 8, color = "chocolate1"))
plot8<-FeaturePlot(all.samples.integrated, reduction = "umap", features = "Ly6c1",coord.fixed = TRUE, cols = c("lightgrey", "Brown")) + NoAxes() + theme(plot.title = element_text(face = "italic", size = 8, color = "chocolate1"))
plot10<-FeaturePlot(all.samples.integrated, reduction = "umap", features = "Adgre1",coord.fixed = TRUE, cols = c("lightgrey", "Brown")) + NoAxes() + theme(plot.title = element_text(face = "italic", size = 8, color = "chocolate1"))
plot16<-FeaturePlot(all.samples.integrated, reduction = "umap", features = "Cd68",coord.fixed = TRUE, cols = c("lightgrey", "Brown")) + NoAxes() + theme(plot.title = element_text(face = "italic", size = 8, color = "chocolate1"))
plot12<-FeaturePlot(all.samples.integrated, reduction = "umap", features = "Cd8a",coord.fixed = TRUE, cols = c("lightgrey", "Brown")) + NoAxes() + theme(plot.title = element_text(face = "italic", size = 8, color = "chocolate1"))
plot13<-FeaturePlot(all.samples.integrated, reduction = "umap", features = "Ccl3",coord.fixed = TRUE, cols = c("lightgrey", "Brown")) + NoAxes() + theme(plot.title = element_text(face = "italic", size = 8, color = "chocolate1"))
#plot14<-FeaturePlot(all.samples.integrated, reduction = "umap", features = "Tnf",coord.fixed = TRUE, cols = c("lightgrey", "Brown")) + NoAxes() + theme(plot.title = element_text(face = "italic", size = 8, color = "chocolate1"))
plot15<-FeaturePlot(all.samples.integrated, reduction = "umap", features = "Il1b",coord.fixed = TRUE, cols = c("lightgrey", "Brown")) + NoAxes() + theme(plot.title = element_text(face = "italic", size = 8, color = "chocolate1"))
leg<-get_legend(plot1)
gg<-ggplot() + theme_minimal()
feature_plot<-plot_grid(plot6 + theme(legend.position = "none"), 
                        plot7 + theme(legend.position = "none"), 
                        plot13 + theme(legend.position = "none"),
                        plot15 + theme(legend.position = "none"),
                        plot10 + theme(legend.position = "none"), 
                        plot16 + theme(legend.position = "none"), 
                        plot8 + theme(legend.position = "none"), 
                        plot3 + theme(legend.position = "none"),
                        plot9 + theme(legend.position = "none"), 
                        plot4 + theme(legend.position = "none"), 
                        plot12 + theme(legend.position = "none"),
                        plot5 + theme(legend.position = "none"), 
                        plot1 + theme(legend.position = "none"),
                        plot2 + theme(legend.position = "none"), 
                        labels = c("                      Microglia","","       Activated Microglia","","                 Macrophages","", "    Monocytes","                         T Cells","","   CD4 T Cells","   CD8 T Cells","      NK Cells", "                         B Cells", ""), ncol = 6, rel_widths = 1, label_size = 8, label_fontface = "italic",label_colour = "brown4", hjust = -0.6)
lege<-plot_grid(gg, gg, leg, gg, gg, ncol = 5, rel_widths = 1)
final_feature_plot<-plot_grid(feature_plot, lege, ncol = 1,nrow = 2, rel_heights = c(5, 0.9))
final_feature_plot
dev.off()


# DimPlot for Myeloid Cells -----------------------------------------------

smap<-levels(subset$Samples)
colr<-c("green4", "firebrick1", "darkred")   ####### to set a colour for samples
names(colr)<-smap
png("myeloid_cells_subset_dimesnionplot.png", units="in", res = 1000, w = 6, h = 6)    ########### dimension plot for myeloid subset objects
dimp<-DimPlot(subset, reduction = "umap", 
              label = F, 
              repel = F, 
              group.by = "Samples", 
              cols = colr) + 
  ggtitle("Myeloid cell subsetted samples") + 
  theme(title = element_text(face = "italic", 
                             size = 14, 
                             hjust = 0.5),
        legend.position = "bottom", 
        panel.border = element_rect(colour = "black", fill=NA, size=1)) + NoAxes()
dimp
dev.off()

# TotalTCells Top 20 Heatmap ----------------------------------------------

subset_TotalTCells<-SetIdent(subset_TotalTCells, value = "Samples")   ###### set ident to samples to calculate DE genes
sample_markers<-FindAllMarkers(subset_TotalTCells, 
                               only.pos = TRUE, 
                               min.pct = 0.1, 
                               test.use = "wilcox",  
                               logfc.threshold = 0.25) ###### find DE genes between samples
top20TotalTCells<- as.data.frame(sample_markers %>% 
                                   group_by(cluster) %>% 
                                   top_n(n = 20, 
                                         wt = avg_log2FC))    ######## get top 20 genes from each samples based on log2FC 
top20TotalTCells_genes<-top20TotalTCells[,c(6,7)] 
rownames(top20TotalTCells_genes)<-top20TotalTCells_genes$gene  ### to get top 20 from each samples

png("TotalTCell Heatmap.png", units="in", res = 1000, w = 10, h = 8)
DoHeatmap(subset_TotalTCells, 
          features = top20TotalTCells$gene, 
          group.by = "Samples",
          group.colors = c("HGG" = "green4", 
                           "LGG" = "firebrick1", 
                           "No Tumor" = "blue")) + 
  scale_fill_gradientn(colors = c("blue", "white", "red"))+
  ggtitle("Combined T Cell Clusters Genes")
dev.off()
write.csv(top20TotalTCells,"top20TotalTCells.csv", row.names = TRUE)

# FindAllMarkers CSV Files for IPA ----------------------------------------

subset_Macrophages_tumor<-SetIdent(subset_Macrophages_tumor, value = "Samples")   ###### set ident to samples to calculate DE genes
sample_markers<-FindAllMarkers(subset_Macrophages_tumor, 
                               only.pos = TRUE, 
                               min.pct = 0.1, 
                               test.use = "wilcox",  
                               logfc.threshold = 0) ###### find DE genes between samples
subset_AAM_tumor<-SetIdent(subset_AAM_tumor, value = "Samples")   ###### set ident to samples to calculate DE genes
sample_markers2_andneg_01<-FindAllMarkers(subset_AAM_tumor, 
                                          only.pos = FALSE, 
                                          min.pct = 0.1, 
                                          test.use = "wilcox",  
                                          logfc.threshold = 0) ###### find DE genes between samples

sample_markers_andneg_01<-FindAllMarkers(subset_Macrophages_tumor, 
                                         only.pos = FALSE, 
                                         min.pct = 0.1, 
                                         test.use = "wilcox",  
                                         logfc.threshold = 0) ###### find DE genes between samples


top20Macrophagestumor<- as.data.frame(sample_markers %>% 
                                 group_by(cluster) %>% 
                                 top_n(n = 20, 
                                       wt = avg_log2FC))    ######## get top 20 genes from each samples based on log2FC 
top20Macrophagestumor_genes<-top20Macrophagestumor[,c(6,7)] 
rownames(top20Macrophagestumor_genes)<-top20Macrophagestumor_genes$gene  ### to get top 20 from each samples

png("Mac1_Tumor_only_DE_heatmap.png", width = 15, height = 10, units = "in", res = 1000)
DoHeatmap(subset_Macrophages_tumor, 
          features = top20Macrophagestumor$gene, 
          group.by = "Samples",
          group.colors = c("HGG" = "green4", 
                           "LGG" = "firebrick1", 
                           "No Tumor" = "blue")) + 
  scale_fill_gradientn(colors = c("blue", "white", "red"))+
  ggtitle("Macrophage Cluster 1 Genes")
dev.off()

write.csv(sample_markers_andneg_01,"tumor_only_macrophage_cluster_1_sample_markers_and_negative_01.csv", row.names = TRUE)
write.csv(sample_markers2_andneg_01,"tumor_only_macrophage_cluster_2_sample_markers_and_negative_01.csv", row.names = TRUE)

subset_microglia_tumor <-SetIdent(subset_microglia_tumor, value = "Samples")
subset_active_microglia_tumor <-SetIdent(subset_active_microglia_tumor, value = "Samples")
subset_totalT_tumor <-SetIdent(subset_totalT_tumor, value = "Samples")
subset_cd4T_tumor <-SetIdent(subset_cd4T_tumor, value = "Samples")
subset_cd8T_tumor <-SetIdent(subset_cd8T_tumor, value = "Samples")

sample_markers_updown_microglia <-FindAllMarkers(subset_microglia_tumor, 
                                                 only.pos = FALSE, 
                                                 min.pct = 0.1, 
                                                 test.use = "wilcox",  
                                                 logfc.threshold = 0) ###### find DE genes between samples
sample_markers_updown_active_microglia<-FindAllMarkers(subset_active_microglia_tumor, 
                                                       only.pos = FALSE, 
                                                       min.pct = 0.1, 
                                                       test.use = "wilcox",  
                                                       logfc.threshold = 0) ###### find DE genes between samples
sample_markers_updown_totalT<-FindAllMarkers(subset_totalT_tumor, 
                                             only.pos = FALSE, 
                                             min.pct = 0.1, 
                                             test.use = "wilcox",  
                                             logfc.threshold = 0) ###### find DE genes between samples
sample_markers_updown_cd4T<-FindAllMarkers(subset_cd4T_tumor, 
                                           only.pos = FALSE, 
                                           min.pct = 0.1, 
                                           test.use = "wilcox",  
                                           logfc.threshold = 0) ###### find DE genes between samples
sample_markers_updown_CD8T<-FindAllMarkers(subset_cd8T_tumor, 
                                           only.pos = FALSE, 
                                           min.pct = 0.1, 
                                           test.use = "wilcox",  
                                           logfc.threshold = 0) ###### find DE genes between samples



write.csv(sample_markers_updown_microglia,"sample_markers_updown_microglia.csv", row.names = TRUE)
write.csv(sample_markers_updown_active_microglia,"sample_markers_updown_active_microglia.csv", row.names = TRUE)
write.csv(sample_markers_updown_totalT,"sample_markers_updown_totalT.csv", row.names = TRUE)
write.csv(sample_markers_updown_cd4T,"sample_markers_updown_cd4T.csv", row.names = TRUE)
write.csv(sample_markers_updown_CD8T,"sample_markers_updown_CD8T.csv", row.names = TRUE)

all.samples.integrated_sample_ident <-SetIdent(all.samples.integrated, value = "Samples")
subset_tumor <-SetIdent(subset_tumor, value = "Samples")
subset_Hgg_no_tumor <-SetIdent(subset_Hgg_no_tumor, value = "Samples")
subset_Lgg_no_tumor <-SetIdent(subset_Lgg_no_tumor, value = "Samples")


subset_Macrophages_all<-subset(all.samples.integrated_sample_ident, subset = final_celltype == "Macrophage Cluster 1")
subset_AAM_all<-subset(all.samples.integrated_sample_ident, subset = final_celltype == "Macrophage Cluster 2")

subset_Macrophages_tumor<-subset(subset_tumor, subset = final_celltype == "Macrophage Cluster 1")
subset_AAM_tumor<-subset(subset_tumor, subset = final_celltype == "Macrophage Cluster 2")

subset_Macrophages_high_no <-subset(subset_Hgg_no_tumor, subset = final_celltype == "Macrophage Cluster 1")
subset_AAM_high_no <-subset(subset_Hgg_no_tumor, subset = final_celltype == "Macrophage Cluster 2")

subset_Macrophages_low_no <-subset(subset_Lgg_no_tumor, subset = final_celltype == "Macrophage Cluster 1")
subset_AAM_low_no <-subset(subset_Lgg_no_tumor, subset = final_celltype == "Macrophage Cluster 2")

sample_markers_updown_Macrophages_all <-FindAllMarkers(subset_Macrophages_all, 
                                                only.pos = FALSE, 
                                                min.pct = 0.1, 
                                                test.use = "wilcox",  
                                                logfc.threshold = 0) ###### find DE genes between samples
sample_markers_updown_AAM_all <-FindAllMarkers(subset_AAM_all, 
                                                only.pos = FALSE, 
                                                min.pct = 0.1, 
                                                test.use = "wilcox",  
                                                logfc.threshold = 0) ###### find DE genes between samples

sample_markers_updown_Macrophages_tumor <-FindAllMarkers(subset_Macrophages_tumor, 
                                                  only.pos = FALSE, 
                                                  min.pct = 0.1, 
                                                  test.use = "wilcox",  
                                                  logfc.threshold = 0) ###### find DE genes between samples
sample_markers_updown_AAM_tumor <-FindAllMarkers(subset_AAM_tumor, 
                                                  only.pos = FALSE, 
                                                  min.pct = 0.1, 
                                                  test.use = "wilcox",  
                                                  logfc.threshold = 0) ###### find DE genes between samples

sample_markers_updown_Macrophages_high_no <-FindAllMarkers(subset_Macrophages_high_no, 
                                                    only.pos = FALSE, 
                                                    min.pct = 0.1, 
                                                    test.use = "wilcox",  
                                                    logfc.threshold = 0) ###### find DE genes between samples
sample_markers_updown_AAM_high_no <-FindAllMarkers(subset_AAM_high_no, 
                                                    only.pos = FALSE, 
                                                    min.pct = 0.1, 
                                                    test.use = "wilcox",  
                                                    logfc.threshold = 0) ###### find DE genes between samples

sample_markers_updown_Macrophages_low_no <-FindAllMarkers(subset_Macrophages_low_no, 
                                                   only.pos = FALSE, 
                                                   min.pct = 0.1, 
                                                   test.use = "wilcox",  
                                                   logfc.threshold = 0) ###### find DE genes between samples
sample_markers_updown_AAM_low_no <-FindAllMarkers(subset_AAM_low_no, 
                                                   only.pos = FALSE, 
                                                   min.pct = 0.1, 
                                                   test.use = "wilcox",  
                                                   logfc.threshold = 0) ###### find DE genes between samples

write.csv(sample_markers_updown_Macrophages_all,"sample_markers_updown_Macrophages_all.csv", row.names = TRUE)
write.csv(sample_markers_updown_AAM_all,"sample_markers_updown_AAM_all.csv", row.names = TRUE)
write.csv(sample_markers_updown_Macrophages_tumor,"sample_markers_updown_Macrophages_tumor.csv", row.names = TRUE)
write.csv(sample_markers_updown_AAM_tumor,"sample_markers_updown_AAM_tumor.csv", row.names = TRUE)
write.csv(sample_markers_updown_Macrophages_high_no,"sample_markers_updown_Macrophages_high_no.csv", row.names = TRUE)
write.csv(sample_markers_updown_AAM_high_no,"sample_markers_updown_AAM_high_no.csv", row.names = TRUE)
write.csv(sample_markers_updown_Macrophages_low_no,"sample_markers_updown_Macrophages_low_no.csv", row.names = TRUE)
write.csv(sample_markers_updown_AAM_low_no,"sample_markers_updown_AAM_low_no.csv", row.names = TRUE)


# Heatmaps for Figure Panels ----------------------------------------------------------------

# Differences between macrophage clusters
features_macs <- c('Tgfbi','Cd14','Itga4','Cd163','Ifitm3','Adgre1','Trem2','Cd9',
                   'Cd81','Hvcn1','Spp1','Ctsb','Ctsd','Ctsz','Lgals3','Cd63',
                   'Cd68','Syngr1','Timp2','C1qa','C1qb','C1qc') 

subset_macrophages <- subset(all.samples.integrated,
                             subset = final_celltype == c('Macrophage Cluster 1',
                                                          'Macrophage Cluster 2'))

subset_macrophages <- SetIdent(subset_macrophages, 
                               value = "final_celltype")

DoHeatmap(subset_macrophages,features = features_macs, label = FALSE) + 
          scale_fill_gradientn(colors = c("blue", "white", "red")) +
          ggtitle("Macrophage Genes")

# T cells

subset_TotalTCells <- SetIdent(subset_macrophages, 
                               value = "Samples")

features_macs <- c('Ccl8','Cxcr6','Ifng','Ccl5','Fabp5','Apoe','Apoc1','Hmox1') 

DoHeatmap(subset_TotalTCells,features = features_macs, label = FALSE) + 
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  ggtitle("T Cell Genes")

# Heatmaps for Supplemental Figures ---------------------------------------

# Macrophage Cluster 1
subset_Macrophages<-SetIdent(subset_Macrophages, value = "Samples")   ###### set ident to samples to calculate DE genes
sample_markers<-FindAllMarkers(subset_Macrophages, 
                               only.pos = TRUE, 
                               min.pct = 0.1, 
                               test.use = "wilcox",  
                               logfc.threshold = 0.25) ###### find DE genes between samples
top20Macrophages<- as.data.frame(sample_markers %>% 
                            group_by(cluster) %>% 
                            top_n(n = 20, 
                                  wt = avg_log2FC))    ######## get top 20 genes from each samples based on log2FC 
top20Macrophages_genes<-top20Macrophages[,c(6,7)] 
rownames(top20Macrophages_genes)<-top20Macrophages_genes$gene  ### to get top 20 from each samples
top20Macrophages <- top20Macrophages[c(41:60,21:40,1:20),]

DoHeatmap(subset_Macrophages, 
          features = top20Macrophages$gene, 
          group.by = "Samples",
          group.colors = c("HGG" = "green4", 
                           "LGG" = "firebrick1", 
                           "No Tumor" = "blue")) + 
  scale_fill_gradientn(colors = c("blue", "white", "red"))+
  ggtitle("Macrophage Genes")

# Macrophage Cluster 2
subset_AAM<-SetIdent(subset_AAM, value = "Samples")   ###### set ident to samples to calculate DE genes
sample_markers<-FindAllMarkers(subset_AAM, 
                               only.pos = TRUE, 
                               min.pct = 0.1, 
                               test.use = "wilcox",  
                               logfc.threshold = 0.25) ###### find DE genes between samples
top20AAM<- as.data.frame(sample_markers %>% 
                                   group_by(cluster) %>% 
                                   top_n(n = 20, 
                                         wt = avg_log2FC))    ######## get top 20 genes from each samples based on log2FC 
top20AAM_genes<-top20AAM[,c(6,7)] 
rownames(top20AAM_genes)<-top20AAM_genes$gene  ### to get top 20 from each samples
top20AAM <- top20AAM[c(41:60,21:40,1:20),]

DoHeatmap(subset_AAM, 
          features = top20AAM$gene, 
          group.by = "Samples",
          group.colors = c("HGG" = "green4", 
                           "LGG" = "firebrick1", 
                           "No Tumor" = "blue")) + 
  scale_fill_gradientn(colors = c("blue", "white", "red"))+
  ggtitle("Alternatively Activated Macrophage Genes")

# T Cells
subset_TotalTCells<-SetIdent(subset_TotalTCells, value = "Samples")   ###### set ident to samples to calculate DE genes
sample_markers<-FindAllMarkers(subset_TotalTCells, 
                               only.pos = TRUE, 
                               min.pct = 0.1, 
                               test.use = "wilcox",  
                               logfc.threshold = 0.25) ###### find DE genes between samples
top20TotalTCells<- as.data.frame(sample_markers %>% 
                           group_by(cluster) %>% 
                           top_n(n = 20, 
                                 wt = avg_log2FC))    ######## get top 20 genes from each samples based on log2FC 
top20TotalTCells_genes<-top20TotalTCells[,c(6,7)] 
rownames(top20TotalTCells_genes)<-top20TotalTCells_genes$gene  ### to get top 20 from each samples
top20TotalTCells <- top20TotalTCells[c(41:60,21:40,1:20),]

DoHeatmap(subset_TotalTCells, 
          features = top20TotalTCells$gene, 
          group.by = "Samples",
          group.colors = c("HGG" = "green4", 
                           "LGG" = "firebrick1", 
                           "No Tumor" = "blue")) + 
  scale_fill_gradientn(colors = c("blue", "white", "red"))+
  ggtitle("T Cell Genes")



# CSV Files for Cluster Identification ------------------------------------

#Overall Differential Expression
sample_markers<-FindAllMarkers(all.samples.integrated, 
                               only.pos = TRUE, 
                               min.pct = 0.1, 
                               test.use = "wilcox",  
                               logfc.threshold = 0.25) ###### find DE genes between samples
write.csv(sample_markers,"sample_markers_total.csv", row.names = TRUE)

#Macrophage Differential Expression
# Macrophage Clusters with final_celltype identity
subset_all_macrophages<-subset(all.samples.integrated, 
                               subset = final_celltype == c('Macrophage Cluster 1','Macrophage Cluster 2'))
sample_markers<-FindAllMarkers(subset_all_macrophages, 
                               only.pos = TRUE, 
                               min.pct = 0.1, 
                               test.use = "wilcox",  
                               logfc.threshold = 0.25) ###### find DE genes between samples
write.csv(sample_markers,"sample_markers_macrophages.csv", row.names = TRUE)

#Total T Cells
subset_TotalTCells <- subset(all.samples.integrated, subset = final_celltype %in% 
                               c("T cells","CD4 T cells","CD8 T cells"))
subset_TotalTCells<-SetIdent(subset_TotalTCells, value = "Samples")   ###### set ident to samples to calculate DE genes
sample_markers<-FindAllMarkers(subset_TotalTCells, 
                               only.pos = TRUE, 
                               min.pct = 0.1, 
                               test.use = "wilcox",  
                               logfc.threshold = 0.25) ###### find DE genes between samples
write.csv(sample_markers,"sample_markers_TotalTCells.csv", row.names = TRUE)


# Generate SingleCellExperiment Object ------------------------------------

asi.sce <- as.SingleCellExperiment(all.samples.integrated)

all.samples.integrated[["RNA"]]@counts
all.samples.integrated[["RNA"]]@data



# add column for percent of hemoglobin genes
all.samples.integrated[["percent.hb"]] <- PercentageFeatureSet(all.samples.integrated, pattern = "^Hb[ab]")

# check summary stats and visualize if any cells have high hemoglobin (indicative of RBCs)
summary(seurat.object$percent.hb)
# or
VlnPlot(all.samples.integrated, features = "percent.hb")

# subset if there are cells that seem to have high hemoglobin
# for example...
seurat.object.subset <- subset(seurat.object, subset = percent.hb < 50)

