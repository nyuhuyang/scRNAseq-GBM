# **Single-cell RNA Sequencing Reveals Immunosuppressive Myeloid Cell Diversity During Malignant Progression in Glioma**

[![arXiv shield](https://img.shields.io/badge/biorxiv-2021.09.24.461735-red.svg?style=flat)](https://www.biorxiv.org/content/10.1101/2021.09.24.461735v1.full)

This project provides the code developed for the study in [Rajendran S et al. Single-cell RNA Sequencing Reveals Immunosuppressive Myeloid Cell Diversity and Restricted Cytotoxic Effector Cell Trafficking and Activation During Malignant Progression in Glioma. 2021](https://www.biorxiv.org/content/10.1101/2021.09.24.461735v1.full)

# Abstract
Low grade gliomas (LGG) 
account for about two-thirds of all glioma diagnoses in adolescents and young adults (AYA) and malignant progression of these patients leads to dismal outcomes. Recent studies have shown the importance of the dynamic tumor microenvironment in high-grade gliomas (HGG), yet its role is still poorly understood in low-grade glioma malignant progression. Here, we investigated the heterogeneity of the immune microenvironment using a platelet-derived growth factor (PDGF)-driven RCAS (replication-competent ASLV long terminal repeat with a splice acceptor) glioma model that recapitulates the malignant progression of low to high-grade glioma in humans and also provides a model system to characterize immune cell trafficking and evolution. To illuminate changes in the immune cell landscape during tumor progression, we performed single-cell RNA sequencing on immune cells isolated from animals bearing no tumor (NT), LGG and HGG, with a particular focus on the myeloid cell compartment, which is known to mediate glioma immunosuppression. LGGs demonstrated significantly increased infiltrating T cells, CD4 T cells, CD8 T cells, B cells, and natural killer cells in the tumor microenvironment, whereas HGGs significantly abrogated this infiltration. Our study identified two distinct macrophage clusters in the tumor microenvironment; one cluster appeared to be bone marrow-derived while another was defined by overexpression of Trem2, a known anti-tumor immunity marker in myeloid cell subpopulations. Our data demonstrates that these two distinct macrophage clusters show an immune-activated phenotype (Stat1, Tnf, Cxcl9 and Cxcl10) in LGG which evolves to an immunosuppressive state (Lgals3, Apoc1 and Id2) in HGG that restricts T cell recruitment and activation. We identified CD74 and macrophage migration inhibition factor (MIF) as potential targets for these distinct macrophage populations. Interestingly, these results were mirrored by our analysis of the TCGA dataset, which demonstrated a statistically significant association between CD74 overexpression and decreased overall survival in AYA patients with grade II LGGs.   Targeting immunosuppressive myeloid cells and intra-tumoral macrophages within this therapeutic window may ameliorate mechanisms associated with immunosuppression before and during malignant progression.


## **Requirements**

* R (tested in R version 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out")
* R librarys: Seurat (v4.0.0), SingleR (v1.4.0), SingleCellExperiment (v1.12.0), ShinyCell (v2.0.0)

## **shinyApp Data**

[Click here to load shinyApp](https://weillcornellmed.shinyapps.io/3_samples_ShinyCell)

![](https://github.com/nyuhuyang/scRNAseq-GBM/blob/main/Figures/UMAP.jpg)

[shinyApp source codes and data](https://www.dropbox.com/sh/s7ewv1s5clmpjua/AAALKvlMATgbxhcrlDEhqiqqa)

## **Reproduce results**

#### **1~2. Data preprocess**
[1 QC.R](https://github.com/nyuhuyang/scRNAseq-GBM/blob/main/R/3_samples/QC.R)
Read 10X matrix and perform quality control.

[2 Seurat_setup](https://github.com/nyuhuyang/scRNAseq-GBM/blob/main/R/3_samples/Seurat_setup.R)
Interagted dataset using harmony and prepare UMAP and tSNE plot.

#### **3~4. Identify cell by SingleR and prepare figures**
[3 SingleR_bioconductor.R](https://github.com/nyuhuyang/scRNAseq-GBM/blob/main/R/3_samples/SingleR_bioconductor.R)
[4 SingleR_figures.R](https://github.com/nyuhuyang/scRNAseq-GBM/blob/main/R/3_samples/SingleR_figures.R)

This script uses Blueprint and Encode Data to to identify cell types automatically

#### **5. Prepare figures**
[5 MT_Analysis.R](https://github.com/nyuhuyang/scRNAseq-GBM/blob/main/R/3_samples/MT_Analysis.R)

We use this script to prepare part of figures. We use above shinyApp to prepare other figures like volcano plots.
