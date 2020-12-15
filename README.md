# IRIS-FGM
---<span style="color:green;">I</span>ntegrative sc<span style="color: red;">R</span>NA-Seq <span style="color: red;">I</span>nterpretation <span style="color: red;">S</span>ystem for <span style="color: red;">F</span>unctional <span style="color: red;">G</span>ene <span style="color: red;">M</span>odule analysis

## Abstract
Single-cell RNA-Seq data is useful in discovering cell heterogeneity and signature genes in specific cell populations in cancer and other complex diseases. Specifically, the investigation of functional gene modules (FGM) can help to understand gene interactive networks and complex biological processes. QUBIC2 is recognized as one of the most efficient and effective tools for FGM identification from scRNA-Seq data. However, its availability is limited to a C implementation, and its applicative power is affected by only a few downstream analyses functionalities. We developed an R package named IRIS-FGM (integrative scRNA-Seq interpretation system for functional gene module analysis) to support the investigation of FGMs and cell clustering using scRNA-Seq data. Empowered by QUBIC2, IRIS-FGM can identify co-expressed and co-regulated FGMs, predict types/clusters, identify differentially expressed genes, and perform functional enrichment analysis. It is noteworthy that IRIS-FGM also applies Seurat objects that can be easily used in the Seurat vignettes.

## Workflow
<img src="https://user-images.githubusercontent.com/26455910/90200476-f64ff900-dda5-11ea-902c-c91726c22eac.jpg" alt="IRISFGM_working_flow">

## Environment

1. System environment:

We recommend users to install IRIS-FGM on large memory (32GB) based UNIX/Linux operation system if the user aims at analyzing bicluster-based co-expression analysis. If the user seeks to analyze data by quick mode, we recommend installing IRIS-FGM on a small memory (8GB) based Windows or Linux operation system; IRIS-FGM does not support MAC. 

2. R environment:
R (equal or greater than 3.5)

# Usage

## Pre-installation
1. Install required packages from CRAN: 

```
install.packages(c('BiocManager','devtools', 'AdaptGauss', "pheatmap", 'mixtools','MCL', 'anocva', 'qgraph','Rtools','ggpubr',"ggraph","Seurat"))
```
                   
2. Install required packages from Bioconductor:  

```
BiocManager::install(c('org.Mm.eg.db','multtest','org.Hs.eg.db','clusterProfiler','DEsingle', 'DrImpute','scater', 'scran'))
```

3. Install IRIS-FGM from github:

```
devtools::install_github("BMEngineeR/IRIS-FGM", force = T)
```

Now the MetaQUBIC is successfully installed and ready for use. 
***
## Data preparation

1. Download example data from [Yan's RPKM 90 cell embroyonic single cell RNA-Seq data](https://bmbl.bmi.osumc.edu/downloadFiles/Yan_expression.txt), [meta file (cell type information)](https://bmbl.bmi.osumc.edu/downloadFiles/Yan_cell_label.txt) and [paper](https://www.nature.com/articles/nsmb.2660). The alternative way to download is to use the step 2 method in this section.

2. Set the working directory where you put data and import the IRIS-FGM library.
```{r setwd, eval =FALSE, echo = TRUE}
setwd("./my_working_dir/")
library(IRISFGM)
# this step is for download Yan's expression data and meta file.
download.file(url = "https://bmbl.bmi.osumc.edu/downloadFiles/Yan_expression.txt",destfile = "./my_working_dir/Yan_expression.txt")
download.file(url = "https://bmbl.bmi.osumc.edu/downloadFiles/Yan_cell_label.txt",destfile = "./my_working_dir/Yan_cell_label.txt")
```
***

## Read in data.
```
InputMatrix <- read.table("./my_working_dir/Yan_expression.txt",header = T, row.names = 1,check.names = F)
Meta_info <- read.table("./my_working_dir/Yan_cell_label.txt",header = T,row.names = 1)
```

IRIS-FGM also provides the function to read in 10X scRNA-Seq data format.

Read HDF5 file```ReadFrom10X_h5()``` or read the folder which contain three files (barcode, sparse matrix, genes)```ReadFrom10X_folder()```

***

## Analysis data
### Data preprocessing and LTMG modeling
1. **Create IRIS-FGM object**:
```{r create_object, eval= FALSE, echo = TRUE,message=FALSE}
object <- CreateIRISCEMObject(InputMatrix)
```
2. **Adding meta information for your cell**:

This step can add the customized cell label by the user. The format of a file passing to `meta.info` is a data frame of which row name should be cell ID, and column name should be cell type.
```{r add_metadata, eval= FALSE, echo = TRUE}
object <- AddMeta(object, meta.info = Meta_info)
```

3. **Filtering out low quality data**:

Use `PlotMeta` and `SubsetData` together to filter out low-quality cells. 
```{r plot_metadata,eval= TRUE, echo = TRUE}
PlotMeta(object)
```
<img src="https://user-images.githubusercontent.com/26455910/90322202-4acab400-df1f-11ea-95a7-94241338c7ca.png" alt="metadata" width="400" height="400">

```
object <- SubsetData(object , nFeature.upper=15000,nFeature.lower=8000,
                         Counts.upper=700000,Counts.lower=400000)
```
After cells get filtered, we can check the distribution again.
```
PlotMeta(object)
```
<img src="https://user-images.githubusercontent.com/26455910/90202006-6eb8b900-ddaa-11ea-9bcd-38ec1401e88b.png" alt="metadata" width="400" height="400">

4. **Normalization**:

Users can choose to perform normalization based on their needs. The normalization method has two options, one is the simplest CPM normalization (default `normalization = 'LibrarySizeNormalization'`). The other is from package scran and can be opened by using parameter `normalization = 'scran'`. Compared to the CPM normalization method, scran will be more accurate but takes more time.
```
object <- ProcessData(object, normalization = "LibrarySizeNormalization")
```
5. **LTMG modeling**:

Here, we will use Left-truncated Mixture Gaussian distribution to model the regulatory signal of each gene. Parameter, 'Gene_use', decides number of top highly variant gene for LTMG modeling, and here we use top 2000 highly variant genes.
```{r run_LTMG, echo = TRUE,eval = FALSE}
object <- RunLTMG(object, Gene_use = "2000", seed = 123)
```
### Biclustering

IRIS-FGM can provide biclustering function, which is based on our in-house novel algorithm, [QUBIC2] (https://github.com/maqin2001/qubic2) to predict functional gene module. 

**Discretization & biclustering**  
If you have more cells to analyze functional gene module you can use LTMG or Quantile based discretization ([QUBIC](https://academic.oup.com/nar/article/37/15/e101/2409951)). In this step, IRIS-FGM will generate files in the local working directory and do not remove them before finishing your analysis (files include tmp_expression.txt, tmp_expression.txt.chars, tmp_expression.txt.chars.blocks, tmp_expression.txt.chars.chars, and tmp_expression.txt.rules).

1. LTMG based discretization: 

You need to run LTMG modeling then binarizing them. It might take a lot of time if you have a large size of cells.
```
object <- RunLTMG(object, Gene_use = "all", seed = 123)
object <- CalBinaryMultiSignal(object)
object <- RunBicluster(object, DiscretizationModel = "LTMG",OpenDual = FALSE,
                          NumBlockOutput = 100, BlockOverlap = 0.7, BlockCellMin = 15)
```
2. Quantile based discretization: 

It might require users to adjust parameter q for deciding cutoff of binarizing. 
```
object <- RunDiscretization(object, q = 0.06)
object <- RunBicluster(object, DiscretizationModel = "Quantile",OpenDual = FALSE, Extension = 0.90,
                          NumBlockOutput = 1000, BlockOverlap = 0.7, BlockCellMin = 15)
```
**Analyzing functional gene module**
This section is based on the quantile based discretization and biclustering results.

1. Visualize the gene module-to-module relationship globally:

The figure via the number of overlap genes (controlled `edge.by = "gene"`) shows that bicluster 20, 23, 24, and 25 have similar gene components. Bicluster 21, 22, 25, and until bicluster 35 have similar gene components. In this figure, we can see gene modules in bicluster 20, 23, 24, and 25 may have a similar component and similar functionality, whereas the gene components from this gene module may differ from the other gene modules from the other biclusters (i.e., bicluster 21, 22, 25, and until bicluster 35)
```
PlotNetwork(object,N.bicluster =c(20:30), edge.by = "gene")
```

<img src="https://user-images.githubusercontent.com/26455910/90253364-0819b680-de0f-11ea-91a6-a61e1df576e8.png" alt="metadata" width="600" height="300">

2. Visualize bicluster 20 & bicluster 35 on heatmap.

From step 1 in "Analysis functional gene module," we postulate gene components from bicluster 20, 23, 24, and 25 may differ from gene components from bicluster 21, 22, 25, and until bicluster 35. Therefore, in this section, we will focus on how the difference is. Therefore, we use bicluster 20 and bicluster 35 to generate heatmap and show such a difference.
```
PlotHeatmap(object,N.bicluster = c(20,35),show.annotation = F)
```
<img src="https://user-images.githubusercontent.com/26455910/90255235-01d90980-de12-11ea-8469-8992f578ee4b.png" alt="metadata" width="400" height="300">


3. Visualize local co-expression gene module network: 

Since we already know the bicluster 20 and bicluster 35 showing the difference at the global level. We then will focus on a local gene module and investigate the co-expression gene network with in the module. Yellow nodes represent the gene module network from bicluster #20. The nodes' size indicates the degree of presence (the number of connected edges for one node). The thickness of edges suggests the value of the correlation coefficient. 

* From this figure (bicluster 20, click the figure for zooming in), we can tell the EIFAD gene show a negative correlation (red color edge) to GOSR1 & BBS5, and also show a positive correlation (grey edge) to ZNF394 & POTEM in a group of cells from this bicluster.
```
PlotModuleNetwork(object, N.bicluster = 20, cutoff=0.6, Node.color = "#E8E504")
```
<img src="https://user-images.githubusercontent.com/26455910/90256029-1a95ef00-de13-11ea-87a4-a4302396df8e.png" alt="metadata" width="600" height="400">

* From this figure (bicluster 35, click the figure for zooming in), we can tell the ROBO1 gene shows a negative correlation (red color edge) to NLRP4 & BPGM. GNPDA2 gene shows a positive correlation to BPGM, KIT, and CCDC25 in a group of cells from this bicluster.
```
PlotModuleNetwork(object, N.bicluster = 35, cutoff=0.6, Node.color = "#E8E504")
```
<img src="https://user-images.githubusercontent.com/26455910/90259954-d6a5e880-de18-11ea-8ade-8d822fde3053.png" alt="metadata" width="600" height="400">

4. Functional enrichment analysis.

ISIR-FGM provide a functional enrichment analysis for a selected gene module. 

1. For gene module from bicluster 20: 

```
object <- RunPathway(object ,module.number =20, selected.gene.cutoff = 0.05, species = "Human", database = "GO", genes.source = "Bicluster")
```


|            | ONTOLOGY | ID         | Description                   | GeneRatio | BgRatio  | pvalue               | p.adjust           | qvalue             | geneID          | Count |
|------------|----------|------------|-------------------------------|-----------|----------|----------------------|--------------------|--------------------|-----------------|-------|
| GO:0008278 | CC       | GO:0008278 | cohesin complex               | 2/34      | 16/19717 | 0.000341142521037922 | 0.0313851119354888 | 0.0280095964641662 | STAG3L3/STAG3L2 | 2     |
| GO:0005689 | CC       | GO:0005689 | U12-type spliceosomal complex | 2/34      | 27/19717 | 0.000986050510033922 | 0.0453583234615604 | 0.0404799683066557 | ZMAT5/SNRNP35   | 2     |
| GO:0000062 | MF       | GO:0000062 | fatty-acyl-CoA binding        | 2/32      | 22/17697 | 0.000715377512795969 | 0.0470414266188569 | 0.038513448693801  | ACBD5/SOAT2     | 2     |
| GO:1901567 | MF       | GO:1901567 | fatty acid derivative binding | 2/32      | 28/17697 | 0.00116271921939764  | 0.0470414266188569 | 0.038513448693801  | ACBD5/SOAT2     | 2     |
| GO:0005484 | MF       | GO:0005484 | SNAP receptor activity        | 2/32      | 31/17697 | 0.001425497776329    | 0.0470414266188569 | 0.038513448693801  | STX11/GOSR1     | 2     |


2. For gene module from bicluster 35: 

```
object <- RunPathway(object ,module.number =35, selected.gene.cutoff = 0.05, species = "Human", database = "GO", genes.source = "Bicluster")
```

|            | ONTOLOGY | ID         | Description                                                  | GeneRatio | BgRatio   | pvalue               | p.adjust           | qvalue             | geneID                            | Count |
|------------|----------|------------|--------------------------------------------------------------|-----------|-----------|----------------------|--------------------|--------------------|-----------------------------------|-------|
| GO:0001667 | BP       | GO:0001667 | ameboidal-type cell migration                                | 6/29      | 461/18670 | 6.43045336163777e-05 | 0.0376874533879247 | 0.0301960094109077 | ROBO1/PTPRG/ERBB4/KIT/PTPRR/HDAC9 | 6     |
| GO:0021889 | BP       | GO:0021889 | olfactory bulb interneuron differentiation                   | 2/29      | 12/18670  | 0.000152281291259449 | 0.0376874533879247 | 0.0301960094109077 | ROBO1/ERBB4                       | 2     |
| GO:0010631 | BP       | GO:0010631 | epithelial cell migration                                    | 5/29      | 351/18670 | 0.000186973532900314 | 0.0376874533879247 | 0.0301960094109077 | ROBO1/PTPRG/KIT/PTPRR/HDAC9       | 5     |
| GO:0090132 | BP       | GO:0090132 | epithelium migration                                         | 5/29      | 354/18670 | 0.000194519550646714 | 0.0376874533879247 | 0.0301960094109077 | ROBO1/PTPRG/KIT/PTPRR/HDAC9       | 5     |
| GO:0090130 | BP       | GO:0090130 | tissue migration                                             | 5/29      | 360/18670 | 0.000210309449709401 | 0.0376874533879247 | 0.0301960094109077 | ROBO1/PTPRG/KIT/PTPRR/HDAC9       | 5     |
| GO:0042692 | BP       | GO:0042692 | muscle cell differentiation                                  | 5/29      | 385/18670 | 0.000286907770833753 | 0.0428448937778404 | 0.0343282631067754 | KIT/SYNE1/TMOD2/UCHL1/HDAC9       | 5     |
| GO:0005001 | MF       | GO:0005001 | transmembrane receptor protein tyrosine phosphatase activity | 2/32      | 17/17697  | 0.000423558866693072 | 0.0222368405013863 | 0.0164965032290986 | PTPRG/PTPRR                       | 2     |
| GO:0019198 | MF       | GO:0019198 | transmembrane receptor protein phosphatase activity          | 2/32      | 17/17697  | 0.000423558866693072 | 0.0222368405013863 | 0.0164965032290986 | PTPRG/PTPRR                       | 2     |

***

### Cell type identification based on MCL

IRIS-FGM provides the function to use the output biclusters to predict cell type based on the MCL algorithm, which predicted results outperformed the other nine bicluster and cluster methods (Details can be found from figure B in IRIS-FGM paper). 
In this case, we identify eight cell clusters.
```
PlotModuleNetwork(object, N.bicluster = 35, cutoff=0.6, Node.color = "#E8E504")
unique(object@MetaInfo$MC_Label)
[1] 1 2 3 4 5 6 7 8
```

### Visualization cell clusters on UMAP

IRIS-FGM integrates the Seurat package for carrying out dimension reduction. Furthermore, IRIS-FGM also integrates the Seurat package cell clustering method as an alternative approach. Before running dimension reduction, you should perform LTMG modeling if you would like to LTMG signal matrix as input rather than the processed expression matrix. The advantage of using the LTMG signal matrix as input is that it can improve dimension reduction and cell clustering. (Detail can be found in this [LTMG paper](https://pubmed.ncbi.nlm.nih.gov/31372654/))

1.  Use processed expression matrix as input to generate UMAP and cell clusters based on the Seurat package.
* runing program:
```
object <- IRISFGM::RunDimensionReduction(object,mat.source = "UMImatrix")
object <- IRISFGM::RunClassification(object,resolution = 0.8)
```
* Visualization by using different cell label:
"Cluster" is the benchmark label from paper. "MC_Label" is the label based on MCL prediction. "Seurat0.6" is the label based on the Seurat cell clustering method. 
```
# This function will ask you to choose an index (represent a cluster identity).
# Use benchmark cell clusters.
PlotDimension(object)
select condition to present
1 : ncount_RNA
2 : nFeature
3 : Cluster
4 : MC_Label
5 : Seurat0.8
select index of cell condition: 3
```
<img src="https://user-images.githubusercontent.com/26455910/90307176-fe8d5e80-dea1-11ea-9440-5d0d4e405bc1.png" alt="metadata" width="600" height="400">

```
# Use MCL predicted cell clusters.
PlotDimension(object)
select condition to present
1 : ncount_RNA
2 : nFeature
3 : Cluster
4 : MC_Label
5 : Seurat0.8
select index of cell condition: 4
```

<img src="https://user-images.githubusercontent.com/26455910/90307184-077e3000-dea2-11ea-9742-fb56b019d114.png" alt="metadata" width="600" height="400">

```
# Use Seurat predicted cell clusters.
PlotDimension(object)
select condition to present
1 : ncount_RNA
2 : nFeature
3 : Cluster
4 : MC_Label
5 : Seurat0.8
select index of cell condition: 5
```

<img src="https://user-images.githubusercontent.com/26455910/90307189-0cdb7a80-dea2-11ea-8849-3be6ece2e104.png" alt="metadata" width="600" height="400">

2.  Use LTMG signal matrix as input to generate UMAP and cell clusters based on the Seurat package.
* runing program:
```
object <- IRISFGM::RunDimensionReduction(object,mat.source = "LTMG")
object <- IRISFGM::RunClassification(object,resolution = 0.8)
```
* Visualization by using different cell label:
"Cluster" is the benchmark label from paper. "MC_Label" is the label based on MCL prediction. "Seurat0.6" is the label based on the Seurat cell clustering method. 
```
# This function will ask you to choose an index (represent a cluster identity).
# Use benchmark cell clusters.
PlotDimension(object)
select condition to present
1 : ncount_RNA
2 : nFeature
3 : Cluster
4 : MC_Label
5 : Seurat0.8
6 : Seurat0.8
select index of cell condition: 3
```
<img src="https://user-images.githubusercontent.com/26455910/90307193-15cc4c00-dea2-11ea-934e-d02fab9150d7.png" alt="metadata" width="600" height="400">

```
# Use MCL predicted cell clusters.
PlotDimension(object)
select condition to present
1 : ncount_RNA
2 : nFeature
3 : Cluster
4 : MC_Label
5 : Seurat0.8
6 : Seurat0.8
select index of cell condition: 4
```
<img src="https://user-images.githubusercontent.com/26455910/90307198-1bc22d00-dea2-11ea-9d2b-c45e0ef37f95.png" alt="metadata" width="600" height="400">

```
# Use Seurat predicted cell clusters.
PlotDimension(object)
select condition to present
1 : ncount_RNA
2 : nFeature
3 : Cluster
4 : MC_Label
5 : Seurat0.8
6 : Seurat0.8
select index of cell condition: 6
```
<img src="https://user-images.githubusercontent.com/26455910/90307201-22e93b00-dea2-11ea-8876-a73a3acec428.png" alt="metadata" width="600" height="400">


### DEG analysis and functional enrichment analysis.

Cell clusters have been identified by a different method. We can use these cell clusters to analyze DEGs and perform functional enrichment analysis.

1. Find DEG based on MCL prediction.

* Find DEGs. 

This function will ask you to provide an index of a cluster label, select a target cluster (group 1), and select a comparing cluster (group 2). Here we select "ML_Label" cell clusters (index is 4), choose the 8th cluster (index is 3, as group 1), and compare to the rest of all clusters (index is 9, as group 2). The DEGs' result can be found in ```object@BiCluster@MarkerGene ```.
```
object <- IRISFGM::FindMarkers(object)
select condition to compare
1 : ncount_RNA
2 : nFeature
3 : Cluster
4 : MC_Label
5 : Seurat0.8
6 : Seurat0.8

select index of cell condition: 4
select index (left) of first group to compare : 
1 : 1
2 : 2
3 : 3
4 : 4
5 : 5
6 : 6
7 : 7
8 : 8
input first group index : 3
select index (left) of second group to compare : 
1 : 1
2 : 2
3 : 3
4 : 4
5 : 5
6 : 6
7 : 7
9 : rest of all

select index of group 2: 9
# DEGs based on MCL predicted label was store in object@BiCluster@MarkerGene
object@BiCluster@MarkerGene[1:5,]
```

|          | LFC       | pval         | pvalue.adj.FDR | 
|----------|-----------|--------------|----------------|
| GNA14    | 1.775559  | 0.000000e+00 | 0.000000e+00   |
| SYTL3    | 1.753460  | 0.000000e+00 | 0.000000e+00   |
| OSBPL10  | 1.682444  | 3.330669e-16 | 2.220446e-13   |
| MIXL1    | 1.465656  | 7.771561e-16 | 3.885781e-13   |
| PRR15    | 1.778150  | 9.992007e-16 | 3.996803e-13   |

* Functional enrichment analysis.

The result can be found at ```object@BiCluster@PathwayFromMC```.

```
object <- RunPathway(object, species = "Human", database = "GO", genes.source = "Bicluster")
object@BiCluster@PathwayFromMC[1:5,]                 
```

| ONTOLOGY   | ID | Description | GeneRatio                | BgRatio  | pvalue    | p.adjust             | qvalue               | geneID               | Count                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |     |
|------------|----|-------------|--------------------------|----------|-----------|----------------------|----------------------|----------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----|
| GO:0006402 | BP | GO:0006402  | mRNA catabolic process   | 110/1417 | 364/18670 | 3.68880405623349e-38 | 2.00855380861913e-34 | 1.69646157070359e-34 | CNOT7/CSDE1/RPL7A/ALKBH5/PATL2/DDX5/PAIP1/IGF2BP3/LSM1/YTHDF2/PSMC6/HNRNPM/LSM6/PSMB1/PSMD2/PSMD8/PSMA3/RPL17/EIF4G1/LSM2/ETF1/PSMA5/RBM46/PSMD14/EIF4A3/PSMA6/HNRNPC/RPL34/PSMA2/HSPA1A/PSMD11/HNRNPU/YBX1/HSPA8/RPL21/EXOSC3/RPS4Y1/PSMD1/RPL14/EXOSC7/RNPS1/PSMB7/PSMD12/HSPA1B/HNRNPR/HSPB1/PPP2CA/RPL4/RPS24/RBM8A/PABPC4/HNRNPD/MAGOH/NCBP2/RPL19/RPL6/EXOSC5/EXOSC9/YWHAZ/RPL13/PSMC2/PSMB3/PSMA1/POLR2G/SERBP1/PSMD9/PCID2/RPS20/XPO1/PSMB6/RPL7/RPL22/RPS7/RPL9/PSMB4/RPS3/MRTO4/RPS6/PABPC1/RPS8/YWHAB/APEX1/RPS3A/RPL10/PSMC3/RPL36A/UBC/PNRC2/RPL29/PSMD4/RPL12/UBB/PSMD6/RPS21/RPL18A/RPS26/RPL15/RPL36/RPL18/RPL38/LSM4/RPL23A/RPL26/RPL10A/RPL37A/SET/EIF3E/RPL30/EXOSC8/RPS19                                                                              | 110 |
| GO:0006401 | BP | GO:0006401  | RNA catabolic process    | 113/1417 | 397/18670 | 1.90474189922193e-36 | 5.18565982063171e-33 | 4.37990387247401e-33 | CNOT7/CSDE1/RPL7A/ALKBH5/PATL2/DDX5/PAIP1/IGF2BP3/LSM1/YTHDF2/PSMC6/HNRNPM/LSM6/PSMB1/PSMD2/PSMD8/PSMA3/RPL17/EIF4G1/LSM2/ETF1/XRN2/PSMA5/RBM46/PSMD14/EIF4A3/PSMA6/HNRNPC/RPL34/PSMA2/HSPA1A/PSMD11/HNRNPU/YBX1/HSPA8/RPL21/DIS3L/EXOSC3/RPS4Y1/PSMD1/RPL14/EXOSC7/RNPS1/PSMB7/PSMD12/HSPA1B/HNRNPR/HSPB1/PPP2CA/RPL4/RPS24/RBM8A/PABPC4/HNRNPD/MAGOH/NCBP2/RPL19/RPL6/EXOSC5/EXOSC9/YWHAZ/RPL13/PSMC2/PSMB3/PSMA1/POLR2G/SERBP1/PSMD9/PCID2/RPS20/XPO1/PSMB6/LIN28A/RPL7/RPL22/RPS7/RPL9/PSMB4/RPS3/MRTO4/RPS6/PABPC1/RPS8/YWHAB/APEX1/RPS3A/RPL10/PSMC3/RPL36A/UBC/PNRC2/RPL29/PSMD4/RPL12/UBB/PSMD6/RPS21/RPL18A/RPS26/RPL15/RPL36/RPL18/RPL38/LSM4/RPL23A/RPL26/RPL10A/RPL37A/SET/EIF3E/RPL30/EXOSC8/RPS19                                                            | 113 |
| GO:0042254 | BP | GO:0042254  | ribosome biogenesis      | 96/1417  | 297/18670 | 4.35948140939203e-36 | 7.91245875804653e-33 | 6.68300851846799e-33 | NOP58/NOL11/RPL7A/SBDS/NOP56/RPF2/DDX47/NOP16/FBL/DDX21/MRPS7/PWP1/ISG20L2/LSM6/RPL7L1/NOLC1/GTPBP4/NHP2/DDX17/C1QBP/DDX10/XRN2/NSA2/EIF4A3/RIOK1/DDX3X/EBNA1BP2/C1D/PAK1IP1/IMP3/RIOK3/UTP18/NIP7/EXOSC3/MRPL20/IMP4/RPL14/EXOSC7/ZNF622/MRPS11/WDR43/EIF6/RPS24/RPP25/RPL6/EXOSC5/CEBPZ/EXOSC9/RPF1/RPL10L/DDX49/MRPS9/TRMT112/RSL1D1/RRN3/XPO1/RPL7/RPS7/MPHOSPH6/POP5/UTP3/NPM3/RSL24D1/ZNHIT3/ZNF593/KRR1/POP7/MRTO4/RPS6/UTP14A/GLUL/PA2G4/GTF3A/RPS8/LTV1/RPL10/WDR18/RPL12/RPS21/EMG1/RPL26L1/RPL38/RPS27L/DCAF13/RPL23A/TFB2M/RPL26/RPL10A/NOP10/MRPL36/MRPL11/MPHOSPH10/MRPL1/EXOSC8/MRPS2/RPS19                                                                                                                                                                 | 96  |
| GO:0008380 | BP | GO:0008380  | RNA splicing             | 120/1417 | 469/18670 | 1.02798676122236e-33 | 1.39934697871394e-30 | 1.18191425257381e-30 | ESRP1/TRA2B/IWS1/SRSF2/RBMX2/SRSF3/SF3B5/HNRNPA1L2/BCAS2/DDX5/SF3A3/DDX47/SNW1/CDC5L/LSM1/NONO/BUD31/SRPK1/HNRNPM/HNRNPA2B1/LSM6/SNRPC/SNRNP25/SNRPB/DDX17/C1QBP/SAP18/GEMIN5/LSM2/HNRNPK/EIF4A3/SRSF11/PCBP1/POLR2L/SRSF7/DHX15/HNRNPC/PSIP1/ISY1/BRDT/SCNM1/RBM7/HSPA1A/SFPQ/CWC25/SNRPA1/RBM3/HNRNPU/LSM10/YBX1/HSPA8/HNRNPA3/ZMAT5/THOC2/THOC7/SLU7/IK/FIP1L1/RNPS1/SYF2/PRPF38A/RBMX/HNRNPR/PPP2CA/U2AF1/RBM8A/DBR1/HNRNPD/PLRG1/MAGOH/PPIL1/PRPF40A/NCBP2/CIR1/SRSF10/POLR2C/POLR2F/HNRNPH1/PNN/PUF60/CLP1/POLR2G/RBM4/SNRPD2/ZNF830/POLR2E/SRSF9/CCAR1/PRDX6/CLK1/MFAP1/PAPOLA/RBM25/CWC15/PRPF31/SNRPD3/CTNNBL1/PPIG/PTBP1/CLK4/SNRPF/PABPC1/RNF113A/SF1/C9orf78/ZCRB1/TGS1/C2orf49/POLR2H/SNUPN/SRSF8/RPS26/STRAP/LSM4/PHF5A/GTF2F2/SRSF1/MPHOSPH10/KHDRBS3/ZMAT2 | 120 |
| GO:0006413 | BP | GO:0006413  | translational initiation | 71/1417  | 193/18670 | 5.74054371002534e-31 | 6.2514521002176e-28  | 5.28009167770542e-28 | RPL7A/EIF4A1/PAIP1/YTHDF2/EIF5/EIF4A2/EIF3I/RPL17/EIF4G1/EIF4E2/DDX3X/RPL34/RPL21/RPS4Y1/EIF3D/RPL14/HSPB1/EIF3B/RPL4/EIF6/RPS24/EIF3J/NCBP2/RPL19/RPL6/EIF1AD/RPL13/ATF4/POLR2G/RBM4/EIF3A/EIF3G/RPS20/RPL7/RPL22/PPP1CA/RPS7/EIF3K/RPL9/CDC123/RPS3/RPS6/PABPC1/EIF2S2/RPS8/RPS3A/EIF3M/RPL10/RPL36A/EIF3F/RPL29/RPL12/RPS21/RPL18A/RPS26/RPL15/RPL36/MTIF2/RPL18/RPL38/RPL23A/RPL26/RPL10A/EIF4G2/RPL37A/EIF3C/EIF1/EIF3E/RPL30/PPP1R15A/RPS19                                                                                                                                                                                                                                                                                                                          | 71  |



2. Find DEG based on Seurat prediction (LTMG as input matrix).
* Find DEGs. 

This function will ask you to provide an index of a cluster label, select a target cluster (group 1), and select a comparing cluster (group 2). Here we select "Seurat0.8" (result based on inputing LTMG signal matrix) cell clusters (index is 6), choose the 1th cluster (index is 2, as group 1), and compare to the rest of all clusters (index is 4, as group 2). The DEGs' result can be found in ```object@LTMG@MarkerGene ```.
```
object <- IRISFGM::FindMarkers(object)
select condition to compare
1 : ncount_RNA
2 : nFeature
3 : Cluster
4 : MC_Label
5 : Seurat0.8
6 : Seurat0.8

select index of cell condition: 6
select index of cell condition: 14
select index (left) of first group to compare : 
1 : 0
2 : 1
3 : 2

input first group index : 2
select index (left) of second group to compare : 
1 : 0
3 : 2
4 : rest of all

select index of group 2: 4

# DEGs based on MCL predicted label was store in object@LTMG@MarkerGene
object@LTMG@MarkerGene[1:5,]
```
|       | LFC              | pval | pvalue.adj.FDR |
|-------|------------------|------|----------------|
| PAIP1 | 2.72862490853613 | 0    | 0              |
| ARG2  | 2.50563425970776 | 0    | 0              |
| TPRXL | 4.14463501839634 | 0    | 0              |
| INTS8 | 2.15470215492115 | 0    | 0              |
| RNF38 | 2.41515880455291 | 0    | 0              |


* Functional enrichment analysis.

The result can be found at ```object@LTMG@Pathway```.

```
object <- RunPathway(object, species = "Human", database = "GO", genes.source = "CTS")
object@LTMG@Pathway[1:5,]            
```
| ONTOLOGY   | ID | Description | GeneRatio                | BgRatio  | pvalue    | p.adjust             | qvalue               | geneID               | Count                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |     |
|------------|----|-------------|--------------------------|----------|-----------|----------------------|----------------------|----------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----|
| GO:0006402 | BP | GO:0006402  | mRNA catabolic process   | 134/1549 | 364/18670 | 3.15782308102538e-53 | 1.7298554837857e-49  | 1.44794498325753e-49 | PAIP1/CSDE1/ETF1/ALKBH5/PSMA6/PSMD8/PSMB1/SERBP1/EIF4G1/HNRNPM/CNOT7/RPL17/PATL2/RPL7A/PSMA5/PSMD12/HSPA1A/YWHAB/RPL6/RPL34/RPL21/HSPB1/LSM2/MAGOH/RPS7/HSPA1B/PSMD1/LSM1/PSMD2/RPL4/RBM46/EXOSC3/RPL19/PCID2/YBX1/RBM8A/PSMA2/IGF2BP3/EXOSC5/HNRNPC/RPL14/RPS4Y1/PSMC6/PSMA1/EXOSC7/RPL22/HNRNPD/PSMC3/PSMD11/XPO1/PABPC1/HSPA8/HNRNPR/POLR2G/RPL13/HNRNPU/EXOSC9/PSMA3/PSMB7/RPL5/PSMD9/YTHDF2/YWHAZ/RNPS1/EXOSC8/SET/RPL9/PSMA7/RPL7/DDX5/DCP1A/PSMB2/RPS24/MRTO4/RPS3A/RPS3/EIF3E/RPS20/APEX1/EIF4A3/LSM4/RPL18/UBB/RPS8/RPL23/RPS21/NCBP2/LSM3/RPL36A/PSMB5/RPL36/RPL37A/PSMB3/PSMC2/RPL38/PSMD14/PSMB6/RPL15/RPL10/RPL29/RPL30/RPL18A/PPP2CA/PRPF18/LSM6/LSM7/RPS6/UBC/PSMD6/RPL26/RPS28/RPS11/RPL10A/PNRC2/NPM1/RPS9/PSMC4/RPLP1/RPL11/RPLP2/RPS29/RPL23A/RPL24/RPL3/RPLP0/SSB/PSMC1/RPSA/RPS19/PABPC4/PSMD7/RPS27/RPS27A/RPL35                                          | 134 |
| GO:0006401 | BP | GO:0006401  | RNA catabolic process    | 138/1549 | 397/18670 | 2.47116937223922e-51 | 6.76853291056321e-48 | 5.66548093972317e-48 | PAIP1/CSDE1/ETF1/ALKBH5/PSMA6/PSMD8/PSMB1/SERBP1/XRN2/EIF4G1/HNRNPM/CNOT7/RPL17/PATL2/RPL7A/PSMA5/PSMD12/HSPA1A/YWHAB/RPL6/RPL34/RPL21/HSPB1/LSM2/MAGOH/RPS7/HSPA1B/PSMD1/LSM1/PSMD2/RPL4/RBM46/EXOSC3/RPL19/PCID2/YBX1/RBM8A/PSMA2/DIS3L/IGF2BP3/EXOSC5/HNRNPC/RPL14/LIN28A/RPS4Y1/PSMC6/PSMA1/EXOSC7/RPL22/HNRNPD/PSMC3/PSMD11/XPO1/PABPC1/HSPA8/HNRNPR/POLR2G/RPL13/HNRNPU/EXOSC9/PSMA3/PSMB7/RPL5/PSMD9/YTHDF2/YWHAZ/RNPS1/EXOSC8/SET/RPL9/PSMA7/RPL7/DDX5/DCP1A/PSMB2/RPS24/MRTO4/RPS3A/RPS3/EIF3E/RPS20/APEX1/EIF4A3/LSM4/RPL18/UBB/RPS8/RPL23/RPS21/NCBP2/LSM3/RPL36A/PSMB5/RPL36/RPL37A/PSMB3/PSMC2/RPL38/PSMD14/PSMB6/RPL15/RPL10/RPL29/RPL30/RPL18A/PPP2CA/PRPF18/LSM6/LSM7/RPS6/UBC/PSMD6/RPL26/RPS28/RPS11/RPL10A/PNRC2/NPM1/RPS9/PSMC4/RPLP1/RPL11/RPLP2/RPS29/RPL23A/RPL24/RPL3/RPLP0/SSB/PSMC1/RPSA/RPS19/PABPC4/PSMD7/RPS27/RPS27A/RPL35/DKC1                   | 138 |
| GO:0006413 | BP | GO:0006413  | translational initiation | 94/1549  | 193/18670 | 3.81465668554073e-50 | 6.96556310779738e-47 | 5.83040158674226e-47 | PAIP1/EIF4A1/EIF4G1/EIF3B/RPL17/RPL7A/EIF3D/RPL6/RPL34/RPL21/HSPB1/EIF6/RPS7/EIF4A2/EIF3I/RPL4/RPL19/PPP1CA/RPL14/RPS4Y1/EIF3M/RPL22/EIF3J/PABPC1/ATF4/POLR2G/RPL13/EIF3A/EIF2S2/EIF3G/EIF4E2/EIF5B/RPL5/EIF3K/YTHDF2/EIF4EBP1/PPP1R15A/RPL9/EIF1AD/EIF1/RPL7/RPS24/RPS3A/MTIF2/RPS3/EIF4G2/EIF3E/EIF3F/RPS20/RPL18/RPS8/EIF1B/RPL23/RPS21/NCBP2/RPL36A/EIF4H/EIF5/RPL36/DDX3X/RPL37A/RPL38/EIF3C/RPL15/RPL10/RPL29/RPL30/RPL18A/EIF2S3/RPS6/RPL26/RPS28/EIF3H/RPS11/RPL10A/NPM1/RPS9/RPLP1/RPL11/DDX1/RPLP2/RPS29/RPL23A/RBM4/RPL24/RPL3/RPLP0/EIF1AX/RPSA/RPS19/RPS27/RPS27A/CDC123/RPL35                                                                                                                                                                                                                                                                                     | 94  |
| GO:0042254 | BP | GO:0042254  | ribosome biogenesis      | 111/1549 | 297/18670 | 6.79251406976065e-45 | 9.3023480185372e-42  | 7.78636612838878e-42 | C1QBP/NSA2/RPL7L1/FBL/EBNA1BP2/PWP1/RIOK1/DDX21/XRN2/RPL7A/RPF2/MRPS7/DDX17/RPL6/UTP18/TRMT112/ISG20L2/NHP2/EIF6/RPS7/WDR43/NOP16/MRPL20/MRPS11/NOLC1/EXOSC3/DDX10/NIP7/SBDS/RIOK3/GTPBP4/EXOSC5/RPL14/NPM3/POP7/GLUL/NOP58/EXOSC7/DDX49/NOL11/GTF3A/DDX47/XPO1/RSL1D1/EMG1/WDR18/ZNF622/KRR1/IMP4/RPP25/EXOSC9/CEBPZ/C1D/NOP56/RPL5/PAK1IP1/IMP3/UTP3/MRPL1/RRS1/EXOSC8/MPHOSPH6/TFB2M/MRPS2/UTP14A/RPL7/RAN/MRPL36/RPS24/MRTO4/RRN3/LTV1/RPL10L/DCAF13/EIF4A3/MRPL10/RPS8/RPS21/ZNHIT3/RPL26L1/POP5/NGDN/MRPL11/DDX3X/RPS27L/LYAR/ZNF593/RPL38/RPF1/PA2G4/RPL10/LSM6/RPS6/RPL26/RPS28/NOB1/WBP11/BYSL/RPL10A/NPM1/RPS9/RPL11/RPL23A/RPL24/RPL3/RPLP0/RPSA/RPS19/RPS27/RPL35/DKC1                                                                                                                                                                                              | 111 |
| GO:0008380 | BP | GO:0008380  | RNA splicing             | 137/1549 | 469/18670 | 4.15036127638313e-41 | 4.54713581440535e-38 | 3.8060997305105e-38  | GEMIN5/RBMX2/C1QBP/IWS1/SRPK1/ESRP1/TRA2B/PSIP1/HNRNPA3/HNRNPM/POLR2E/SRSF2/SF3A3/DDX17/SF3B5/SRSF3/CDC5L/HSPA1A/SNRPC/SRSF7/SNRPB/NONO/LSM2/MAGOH/THOC2/HNRNPK/LSM1/BUD31/PUF60/CIR1/POLR2L/PCBP1/SNRPA1/PAPOLA/RRAGC/HNRNPA2B1/YBX1/RBM8A/DHX15/SNW1/HNRNPC/CWC25/PTBP1/RBM3/PRPF40A/HNRNPD/SFPQ/DDX47/PABPC1/HSPA8/SLU7/BCAS2/FIP1L1/ZCRB1/HNRNPR/SRSF9/POLR2G/THOC7/HNRNPU/SRSF11/SNRPD2/SNRNP25/PPIL1/SCNM1/HNRNPA1L2/LSM10/POLR2F/RBM7/ISY1/SAP18/TGS1/SRSF10/RNPS1/POLR2H/ZNF830/ZMAT5/PRDX6/U2AF1/C2orf49/PRPF38A/DDX5/PQBP1/HNRNPH1/SYF2/POLR2C/CLP1/MFAP1/BRDT/TMBIM6/PRPF31/CWC15/CTNNBL1/RBM25/PNN/CLNS1A/RNF113A/EIF4A3/C9orf78/LSM4/SNRPF/CLK4/NCBP2/LSM3/CCAR1/UBL5/PLRG1/PPIG/SF1/IK/PHF5A/SNRNP27/PCBP2/GTF2F2/SRSF1/ZMAT2/SNRPE/SRSF8/PPP2CA/SRSF5/PRPF18/DBR1/PPIH/LSM6/LSM7/SNRPD3/WBP11/CLK1/LUC7L3/DDX1/POLR2K/RBM4/STRAP/RBMX/SNRPG/PRMT5/DNAJC8/HNRNPA1 | 137 |
