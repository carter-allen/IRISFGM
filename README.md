# IRIS-FGM
Status update on Dec/28/2020:

The new package has passed the Bioconductor automatic check on local computers without any errors and warnings 
<img src="https://user-images.githubusercontent.com/26455910/103252357-80aa8b00-494a-11eb-8dac-973658ba78e9.png" alt="IRISFGM_status1">
<img src="https://user-images.githubusercontent.com/26455910/103252384-a2a40d80-494a-11eb-80a3-dc76af43e6a1.png" alt="IRISFGM_status2">


---<span style="color:green;">I</span>ntegrative sc<span style="color: red;">R</span>NA-Seq <span style="color: red;">I</span>nterpretation <span style="color: red;">S</span>ystem for <span style="color: red;">F</span>unctional <span style="color: red;">G</span>ene <span style="color: red;">M</span>odule analysis

## Abstract
Single-cell RNA-Seq data is useful in discovering cell heterogeneity and signature genes in specific cell populations in cancer and other complex diseases. Specifically, the investigation of functional gene modules (FGM) can help to understand gene interactive networks and complex biological processes. QUBIC2 is recognized as one of the most efficient and effective tools for FGM identification from scRNA-Seq data. However, its availability is limited to a C implementation, and its applicative power is affected by only a few downstream analyses functionalities. We developed an R package named IRIS-FGM (integrative scRNA-Seq interpretation system for functional gene module analysis) to support the investigation of FGMs and cell clustering using scRNA-Seq data. Empowered by QUBIC2, IRIS-FGM can identify co-expressed and co-regulated FGMs, predict types/clusters, identify differentially expressed genes, and perform functional enrichment analysis. It is noteworthy that IRIS-FGM also applies Seurat objects that can be easily used in the Seurat vignettes.

## Workflow
<img src="https://user-images.githubusercontent.com/26455910/103300172-4550a080-49cc-11eb-9b32-15747ea45f79.png" alt="IRISFGM_working_flow">

## Environment

1. System environment:

We recommend users to install IRIS-FGM on large memory (32GB) based UNIX/Linux operation system if the user aims at analyzing bicluster-based co-expression analysis. If the user seeks to analyze data by quick mode, we recommend installing IRIS-FGM on a small memory (8GB) based Windows or Linux operation system; IRIS-FGM does not support MAC. 

2. R environment:
R (equal or greater than 3.5)

# Usage

## Pre-installation
1. Install required packages from CRAN: 

```
install.packages(c('BiocManager','devtools', 'AdaptGauss', "pheatmap", 'mixtools','MCL', 'anocva', 'qgraph','Rtools','ggpubr',"ggraph","Seurat","Polychrome"))
```
                   
2. Install required packages from Bioconductor:  

```
BiocManager::install(c('org.Mm.eg.db','multtest','org.Hs.eg.db','clusterProfiler','DEsingle', 'DrImpute','scater', 'scran'))
```

3. Install IRISFGM from github:

```
devtools::install_github("BMEngineeR/IRISFGM")
```

Now the MetaQUBIC is successfully installed and ready for use. 
***
## Data preparation

# 1. Input data, create IRISCEM object, add meta information, and preprocessing. 

IRIS-FGM can accepted 10X chromium input files, including a folder (contain gene name, cell name, and sparse matrix) and .h5 file.

## Input data

1. set working directory and import library
```{r setwd, eval =TRUE, echo = TRUE}
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

```{r txt, eval= TRUE, echo = TRUE}
InputMatrix <- read.table(url("https://bmbl.bmi.osumc.edu/downloadFiles/Yan_expression.txt"),
                          header = TRUE, 
                          row.names = 1,
                          check.names = FALSE)
```

## Add meta information

1. For the computational efficiency, we will use subsampling data, and create IRIS-FGM object.

```{r create_object, eval= TRUE, echo = TRUE,message=TRUE}
set.seed(123)
seed_idx <- sample(1:nrow(InputMatrix),3000)
InputMatrix_sub <- InputMatrix[seed_idx,]
object <- CreateIRISFGMObject(InputMatrix_sub)
```

2. Addmeta: this step can add customized cell label by user, the format of file passing to `meta.info` is data frame of which row name should be cell ID, and column name should be cell type.    
```{r add_metadata, eval= TRUE, echo = TRUE}
my_meta <- read.table(url("https://bmbl.bmi.osumc.edu/downloadFiles/Yan_cell_label.txt"),header = TRUE,row.names = 1)
object <- AddMeta(object, meta.info = my_meta)
```

3. plotmeta: plot meta information based on RNA count and Feature number. This step is for the following subset step in terms of filtering out low quality data.    
```{r plot_metadata,eval= TRUE, echo = TRUE}
PlotMeta(object)
```
<img src="https://user-images.githubusercontent.com/26455910/103784216-ac430a80-5007-11eb-8709-c22fa746ef8d.png" alt="IRISFGM_plotmeta">

4. remove low quality data based on the previous plot.
```{r subset_data,eval= TRUE, echo =  TRUE}
object <- SubsetData(object , nFeature.upper=2000,nFeature.lower=250)
```

## Preprocesing 

User can choose perform normalization or imputation based on their need. The normalization method has two options, one is the simplist CPM normalization (default `normalization = 'cpm'`). The other is from package scran and can be opened by using parameter `normalization = 'scran'`, . The imputation method is from package DrImpute and can be opened by using parameter `IsImputation = TRUE` (default as closed).
```{r ProcessData,echo = TRUE, eval= TRUE}
object <- ProcessData(object, normalization = "cpm", IsImputation = FALSE)
```

# 2. Run LTMG

The argument `Gene_use = 500` is  top 500 highlt variant genes which are selected to run LTMG. For quick mode, we recommend to use top 2000 gene (here we use top 500 gene for saving time). On the contrary, for co-expression gene analysis, we recommend to use all gene by changing `Gene_use = "all"`. 
```{r run_LTMG, echo = TRUE,eval = FALSE}
# do not show progress bar
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
# demo only run top 500 gene for saving time.
object <- quiet(RunLTMG(object, Gene_use = "500"))
# you can get LTMG signal matrix
LTMG_Matrix <- GetLTMGmatrix(object)
LTMG_Matrix[1:5,1:5]
```
| ID       | OoCyte_1  | OoCyte_2 | OoCyte_3 | Zygote_2 | 2_Cell_embryo_1_Cell_1 |
|----------|-----------|----------|----------|----------|------------------------|
| MTRNR2L2 | 1         | 2        | 2        | 2        | 2                      |
| TPT1     | 2         | 2        | 2        | 2        | 3                      |
| UBB      | 2         | 2        | 2        | 2        | 2                      |
| FABP5    | 3         | 2        | 2        | 2        | 2                      |
| PTTG1    | 2         | 2        | 2        | 2        | 2                      |

# 3. Biclustering

IRIS-FGM can provide biclustering function, which is based on our in-house novel algorithm, 
QUBIC2 (<https://github.com/maqin2001/qubic2>). Here we will show the basic biclustering 
usage of IRIS-FGM using a $500\times 87$ expression matrix generated from previous top 500 variant genes. 
However, we recommend user should use `Gene_use = all` to generate LTMG matrix. 

## LTMG-discretized bicluster (recommend for small single cell RNA-seq data)
User can type the following command to run discretization (LTMG) + biclustering directly:
```{r biclustering_basedLTMG,eval= FALSE,echo = TRUE}
object <- CalBinaryMultiSignal(object)
object <- RunBicluster(object, DiscretizationModel = "LTMG",OpenDual = FALSE,
                       NumBlockOutput = 100, BlockOverlap = 0.5, BlockCellMin = 25)

```

## Quantile-discretized bicluster (recommend for bulk RNA-Seq, microarray data, or large single cell RNA-Seq data)

This will output several files, and among them you will find one named  `tmp_expression.txt.chars.blocks`,which contains the predicted biclusters.
Or, user may use first version discretization strategy provided by QUBIC 1.0.
```{r biclustering_basedQuantile,eval=FALSE,echo = TRUE}
object <- RunDiscretization(object)
object <- RunBicluster(object, DiscretizationModel = "Quantile",OpenDual = FALSE,
                       NumBlockOutput = 100, BlockOverlap = 0.5, BlockCellMin = 25)
```

(The default parameters in IRIS-FGM are BlockCellMin=15, BlockOverlap=0.7,
Extension=0.90, NumBlockOutput=100 you may use other parameters as you like, just specify them in the argument)


# 4. Cell clustering

## 4.1 Perform dimension Reduction and implement Seurat clustering method.
User can use `reduction = "umap"` or `reductopm = "tsne"` to perform dimension reduction. 
```{r Run_dimReduce, eval= FALSE, echo = TRUE}
# demo only run top 500 gene for saving time.
object <- RunDimensionReduction(object, mat.source = "UMImatrix",reduction = "umap")
object <- RunClassification(object, k.param = 20, resolution = 0.8, algorithm = 1)
```

## 4.2 Predict cell clusters based on Markove clustering

The cell cluster prediction of IRIS-FGM is based on the biclustering results. 
In short, it will construct a weighted graph based on the biclusters and then do clustering on the weighted graph. To facilitate the process, we will use the pre-generated object to perform cell clustering result based on biclustering results.
```{r cell_type, eval=FALSE, echo =TRUE}
object <- FindClassBasedOnMC(object)
```
```{r load_example_object, eval= TRUE, echo = TRUE}
data("example_object")
example_object@MetaInfo[1:5,]
```
|          | ncount_RNA | nFeature | Cluster | Seurat_r_0.8_k_20 | MC_Label |
|----------|------------|----------|---------|-------------------|----------|
| OoCyte_1 | 28202.91   | 657      | 1       | 2                 | 1        |
| OoCyte_2 | 27029.84   | 660      | 1       | 2                 | 1        |
| OoCyte_3 | 26520.14   | 672      | 1       | 2                 | 1        |
| Zygote_1 | 24796.04   | 681      | 2       | 2                 | 1        |
| Zygote_2 | 23575.54   | 679      | 2       | 2                 | 1        |

# 5. Visualization and interpretation
## 5.1 Bicluster results
### Check gene overlapping of FGMs. The results show first 19 FMGs have overlapping genes, and there is no overlapping gene between first 1 FMGs and the 15th FGM.
```{r bicluster_network, eval=TRUE, echo =TRUE}
PlotNetwork(example_object,N.bicluster = c(1:20))
```
<img src="https://user-images.githubusercontent.com/26455910/103784903-84a07200-5008-11eb-9ce1-b63c47f6232e.png" alt="PlotNetwork">

### Perform pathway enrichment analysis on the selected the 1st FGM and visualize the results.
```{r bicluster_pathway, eval=TRUE, echo =TRUE}
example_object <- RunPathway(example_object, N.bicluster =4, species = "Human",database = "GO",genes.source = "Bicluster")
DotPlotPathway(example_object,genes.source = "Bicluster")
```
<img src="https://user-images.githubusercontent.com/26455910/103785140-d9dc8380-5008-11eb-89f7-dec559d5392e.png" alt="DotPlotPathway">

### Heatmap shows relations between any two biclusters (1th and 15th bicluster). 
```{r bicluster_heatmap, eval=TRUE, echo =TRUE}
PlotHeatmap(example_object,N.bicluster = c(1, 20),show.clusters = TRUE,show.annotation=TRUE)
```
<img src="https://user-images.githubusercontent.com/26455910/103785345-1314f380-5009-11eb-9eee-d44639f7ab45.png" alt="PlotHeatmap">

### Co-expression network shows gene correlation relations in one select FGM. 
```{r bicluster_network_module, eval=TRUE, echo =TRUE}
PlotModuleNetwork(example_object,N.bicluster = 3,method = "spearman",
                  cutoff.neg = -0.5,
                  cutoff.pos = 0.5,
                  layout = "circle",
                  node.label = TRUE,
                  node.col = "black",
                  node.label.cex = 10)

```
<img src="https://user-images.githubusercontent.com/26455910/103785450-32138580-5009-11eb-8f70-8310027b26b8.png" alt="PlotModuleNetwork">

## 5.1 cell clustering results
### Visualize cell clustering results
```{r cell_clustering_umap1, eval=TRUE, echo =TRUE}
# cell clustering results based on Seurat clustering method 
PlotDimension(example_object,  reduction = "umap",idents = "Seurat_r_0.8_k_20")
```
<img src="https://user-images.githubusercontent.com/26455910/103785696-86b70080-5009-11eb-9b58-aed94deb1a9a.png" alt="cell_clustering_umap1">

```{r cell_clustering_umap2, eval=TRUE, echo =TRUE}
# cell clustering results based on MCL clustering method 
PlotDimension(example_object, reduction = "umap",idents = "MC_Label")
```
<img src="https://user-images.githubusercontent.com/26455910/103785902-bf56da00-5009-11eb-8010-192eb7d73268.png" alt="cell_clustering_umap2">

### Find global marker based on Seurat method 
```{r cell_clustering_globalmarker, eval=TRUE, echo =TRUE}
global_marker <- FindGlobalMarkers(example_object,idents = "Seurat_r_0.8_k_20")
PlotMarkerHeatmap(Globalmarkers = global_marker,object = example_object,idents ="Seurat_r_0.8_k_20")
```
<img src="https://user-images.githubusercontent.com/26455910/103786020-ea412e00-5009-11eb-9311-050a95d2a625.png" alt="PlotMarkerHeatmap">
