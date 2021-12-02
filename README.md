# **Single-cell RNA Sequencing Reveals Immunosuppressive Myeloid Cell Diversity During Malignant Progression in Glioma**

This project provides the code developed for the study in [Rajendran S et al. Single-cell RNA Sequencing Reveals Immunosuppressive Myeloid Cell Diversity and Restricted Cytotoxic Effector Cell Trafficking and Activation During Malignant Progression in Glioma. 2021](https://www.biorxiv.org/content/10.1101/2021.09.24.461735v1.full)


## **Requirements**

* R (tested in R version 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out")
* R librarys: Seurat (v4.0.0), SingleR (v1.4.0), SingleCellExperiment (v1.12.0), ShinyCell (v2.0.0)

## **shinyApp Data**

[shinyApp](https://weillcornellmed.shinyapps.io/3_samples_ShinyCell)

![](https://github.com/nyuhuyang/scRNAseq-GBM/blob/main/Figures/UMAP.jpg)
![](https://github.com/nyuhuyang/scRNAseq-MouseSkinEpithelia/blob/master/Figs/F1B_dotplot.jpeg)

[shinyApp source codes and data for above ](https://www.dropbox.com/sh/s7ewv1s5clmpjua/AAALKvlMATgbxhcrlDEhqiqqa)

## **Reproduce results**

#### **1. Data preprocess**
[1 setup.R](https://github.com/nyuhuyang/scRNAseq-GBM/blob/main/R/3_samples/setup.R)

#### **2-3. Identify cell by SingleR and prepare figures**
[2 SingleR_bioconductor.R](https://github.com/nyuhuyang/scRNAseq-GBM/blob/main/R/3_samples/SingleR_bioconductor.R)
[3 SingleR_figures.R](https://github.com/nyuhuyang/scRNAseq-GBM/blob/main/R/3_samples/SingleR_figures.R)

This script uses Blueprint and Encode Data to to identify cell types automatically

