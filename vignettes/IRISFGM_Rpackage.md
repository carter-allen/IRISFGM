---
title: "IRIS-FGM vignette"
subtitle: IRIS-FGM 
abstract: IRIS-FGM, integrative and interpretation system for co-expression module analysisa biclustering-based gene regulation inference and cell type prediction method for single-cell RNA-Seq data. This R package integrates in-house computational tools and provides two types of analysis, including QUBIC2 based co-expression analysis and LTMG (left-truncated mixture Gaussian model) based scRNA-Seq analysis (quick mode). IRIS-FGM contains fourfour  major steps; (i) data preprocessing and regulatory signal modelingd LTMG modeling; (ii) co-regulated expression gene module identification; (iii) cell clustering; (iv) co-expression module and differentially expressed gene analysis. 
author: Yuzhou Chang
date: "12 Aug, 2020"
output:
  BiocStyle::html_document:
    number_sections: no
    toc: yes
    highlight: pygments
vignette: >
  %\VignetteIndexEntry{IRIS-FGM vignette}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---


# Intorduction to IRIS-FGM

## General introduction 
IRIS-FGM integrates in-house and state-of-the-art computational tools and provides two analysis strategies, including bicluster-based co-expression gene analysis [(Xie, et al., 2020)](https://academic.oup.com/bioinformatics/article-abstract/36/4/1143/5567116?redirectedFrom=fulltext) and LTMG (left-truncated mixture Gaussian model)-embedded scRNA-Seq analysis [(Wan, et al., 2019)](https://academic.oup.com/nar/article/47/18/e111/5542876).

## Main function

The main idea of IRIS-FGM consists of two major strategies: 

* (i)  biclustering 
* (ii) quick mode (cluster) 



## Obejct structure 
The computational frame is constructed under the S4 data structure in R. The structure of `BRIC object` is: 

- -**BRIC_Object:** *name of object is called BRIC.*
  - -**raw_count:** *raw data matrix (gene in row, cell in columns, prefer using gene symbol).*
  - -**processed_count:** *normalized and imputation (default: FALSE).*
  - -**Meta_info:** *cell classification information based on LTMG and Bicluster.*
  - -**Discretization:** *discretized matrix based on qubic 1.0, which prepares for microarry and bulk RNA-Seq analysis.*
  - -**LTMG:** *LTMG slot is for storing relative results from first strategy.*
    - -**LTMG_discrete:** *Condition assigned matrix, which is generating from LTMG model.*
    - -**LTMG_BinarySingleSignal:** *binary matrix based on gene “on /off.”*
    - -**LTMG_BinaryMultisignal:** *binary matrix based on multiple condition.*
    - -**DimReduce:** *include three dimension loading score, including PCA, Tsne, and UMAP*
    - -**MarkerGene:** *Marker gene based on cell type identified by Seurat clustering method.*
    - -**Pathway:** *based on marker gene.*
    - -**tmp.Seurat:** *temporary Seurat object. In this Seurat Object, the starting matrix is LTMG signalling matrix.*
  - -**Bicluster:** **
    - -**Coreg_gene:** *co-regulatory genes are stored in this slot as dataframe; the first column is gene name and the second column is module number.*
    - -**CoCond_cell:** *co-condition cell are stored in this slot as dataframe; the first column is cell name and the second column is module number.*
    - -**MarkerGene:** *Marker gene based on cell type identified by Markov chain clustering algorithm.*
    - -**Pathway:** *genes based on co-expression gene of gene module (from Coreg_gene).*



# Requirements
## Environment

We recommend user to install IRIS-FGM on large memory (32GB) based linux operation system if user aims at analyzing bicluster-based co-expression analysis; if user aims at analyzing data by quick mode, we recommend to install IRIS-FGM on small memeory (8GB) based Windows or linux operation system; IRIS-FGM does not support MAC. 
We will assum you have the following installed:

* R (equal or greater than 3.5)

Pre-install packge
`install.packages(c('BiocManager','devtools', 'AdaptGauss', "pheatmap", 'mixtools','MCL', 'anocva', "Polychrome", 'qgraph','Rtools','ggpubr',"ggraph"))`
                   
`BiocManager::install(c('org.Mm.eg.db','multtest', 'org.Hs.eg.db','clusterProfiler','DEsingle', 'DrImpute', 'scater', 'scran'))`
                       
`devtools::install_github(repo = 'satijalab/seurat')`

## Input

1. The input to IRIS-FGM is the single-cell RNA-seq expression matrix:

+ Rows correspond to genes and columns correspond to cells.
+ Expression units: the preferred expression values are RPKM/FPKM/CPM. 
+ The data file should be tab delimited.

2. IRIS-FGM also accepts output files from 10X CellRanger, includinhg a folder which contains three individual files and h5 file. 

## Others

When you perform co-expression analysis, it will output several intermediate files, thus please make sure that you have write permission to the folder where IRIS-FGM is located. 

# Installation

For installation, simply type the following command in your R console, please select option 3 when R asks user to update packages:

`devtools::install_github("BMEngineeR/IRISFGM", force = TRUE)`


# Example dataset

This tutorial run on a real dataset to illustrate the results obtained at each step.

As example, we will use Yan's data, a dataset containing 90 cells and 20,214 genes from human embryo, to conduct cell type prediction.

> Yan, L. et al. Single-cell RNA-Seq profiling of human preimplantation embryos and embryonic stem cells. Nat. Struct. Mol. Biol. 20, 1131-1139 (2013)

The original expression matrix was downloaded from <https://s3.amazonaws.com/scrnaseq-public-datasets/manual-data/yan/nsmb.2660-S2.csv>. The expression is provided as RPKM value. For convenience, we removed the space in the column names and deleted the second column(Transcript_ID). The processed data is available at <https://bmbl.bmi.osumc.edu/downloadFiles/Yan_expression.txt>.


# 1. Input data, create IRISCEM object, add meta information, and preprocessing. 

IRIS-FGM can accepted 10X chromium input files, including a folder (contain gene name, cell name, and sparse matrix) and .h5 file.

## Input data

1. set working directory and import library

```r
# dir.create("your working directory",recursive = TRUE)
# setwd("your working directory")
library(IRISFGM)
```

2. Read from .h5 file. 
`ReadFrom10X_h5("~/5k_pbmc_protein_v3_filtered_feature_bc_matrix.h5")`

3. Read from 10X folder, which should contain three files (barcode, gene name, and sparse matrix)

`ReadFrom10X_folder("~/hg19/")`

4. Read from .csv or .txt file

First, we should download data from the link, then we will use this data set as example to run the pipeline.


```r
download.file("https://bmbl.bmi.osumc.edu/downloadFiles/Yan_cell_label.txt",destfile = "./Yan_cell_label.txt")
InputMatrix <- read.table(url("https://bmbl.bmi.osumc.edu/downloadFiles/Yan_expression.txt"),
                          header = TRUE, 
                          row.names = 1,
                          check.names = FALSE)
```

## Add meta information

1. For the computational efficiency, we will use subsampling data, and create IRIS-FGM object.


```r
set.seed(123)
seed_idx <- sample(1:nrow(InputMatrix),3000)
InputMatrix_sub <- InputMatrix[seed_idx,]
object <- CreateIRISFGMObject(InputMatrix_sub)
```

```
## Creating IRISCEM object. 
## The original input file contains 90 cells and 3000 genes 
## Removed 97 genes that total expression value is equal or less than 0
## Removed 0 cells that number of expressed gene is equal or less than 0
```

2. Addmeta: this step can add customized cell label by user, the format of file passing to `meta.info` is data frame of which row name should be cell ID, and column name should be cell type.    































