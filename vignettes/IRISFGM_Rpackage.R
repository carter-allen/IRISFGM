## ----setwd, eval =TRUE, echo = TRUE-------------------------------------------
# dir.create("your working directory",recursive = TRUE)
# setwd("your working directory")
library(IRISFGM)

## ----txt, eval= TRUE, echo = TRUE---------------------------------------------
InputMatrix <- read.table(url("https://bmbl.bmi.osumc.edu/downloadFiles/Yan_expression.txt"),
                          header = TRUE, 
                          row.names = 1,
                          check.names = FALSE)

## ----create_object, eval= TRUE, echo = TRUE,message=TRUE----------------------
set.seed(123)
seed_idx <- sample(1:nrow(InputMatrix),3000)
InputMatrix_sub <- InputMatrix[seed_idx,]
object <- CreateIRISFGMObject(InputMatrix_sub)

## ----add_metadata, eval= TRUE, echo = TRUE------------------------------------
my_meta <- read.table(url("https://bmbl.bmi.osumc.edu/downloadFiles/Yan_cell_label.txt"),header = TRUE,row.names = 1)
object <- AddMeta(object, meta.info = my_meta)

## ----plot_metadata,eval= TRUE, echo = TRUE------------------------------------
PlotMeta(object)

## ----subset_data,eval= TRUE, echo =  TRUE-------------------------------------
object <- SubsetData(object , nFeature.upper=2000,nFeature.lower=250)

## ----ProcessData,echo = TRUE, eval= TRUE--------------------------------------
object <- ProcessData(object, normalization = "cpm", IsImputation = FALSE)

## ----run_LTMG, echo = TRUE,eval = FALSE---------------------------------------
#  # do not show progress bar
#  quiet <- function(x) {
#    sink(tempfile())
#    on.exit(sink())
#    invisible(force(x))
#  }
#  # demo only run top 500 gene for saving time.
#  object <- quiet(RunLTMG(object, Gene_use = "500"))
#  # you can get LTMG signal matrix
#  LTMG_Matrix <- GetLTMGmatrix(object)
#  LTMG_Matrix[1:5,1:5]

## ----biclustering_basedLTMG,eval= FALSE,echo = TRUE---------------------------
#  object <- CalBinaryMultiSignal(object)
#  object <- RunBicluster(object, DiscretizationModel = "LTMG",OpenDual = FALSE,
#                         NumBlockOutput = 100, BlockOverlap = 0.5, BlockCellMin = 25)
#  

## ----biclustering_basedQuantile,eval=FALSE,echo = TRUE------------------------
#  object <- RunDiscretization(object)
#  object <- RunBicluster(object, DiscretizationModel = "Quantile",OpenDual = FALSE,
#                         NumBlockOutput = 100, BlockOverlap = 0.5, BlockCellMin = 25)

## ----Run_dimReduce, eval= FALSE, echo = TRUE----------------------------------
#  # demo only run top 500 gene for saving time.
#  object <- RunDimensionReduction(object, mat.source = "UMImatrix",reduction = "umap")
#  object <- RunClassification(object, k.param = 20, resolution = 0.8, algorithm = 1)

## ----cell_type, eval=FALSE, echo =TRUE----------------------------------------
#  object <- FindClassBasedOnMC(object)

## ----load_example_object, eval= TRUE, echo = TRUE-----------------------------
data("example_object")
example_object@MetaInfo[1:5,]

## ----bicluster_network, eval=TRUE, echo =TRUE---------------------------------
PlotNetwork(example_object,N.bicluster = c(1:20))

## ----bicluster_pathway, eval=TRUE, echo =TRUE---------------------------------
example_object <- RunPathway(example_object, N.bicluster =4, species = "Human",database = "GO",genes.source = "Bicluster")
DotPlotPathway(example_object,genes.source = "Bicluster")

## ----bicluster_heatmap, eval=TRUE, echo =TRUE---------------------------------
PlotHeatmap(example_object,N.bicluster = c(1, 20),show.clusters = TRUE,show.annotation=TRUE)

## ----bicluster_network_module, eval=TRUE, echo =TRUE--------------------------
PlotModuleNetwork(example_object,N.bicluster = 3,method = "spearman",
                  cutoff.neg = -0.5,
                  cutoff.pos = 0.5,
                  layout = "circle",
                  node.label = TRUE,
                  node.col = "black",
                  node.label.cex = 10)


## ----cell_clustering_umap, eval=TRUE, echo =TRUE------------------------------
# cell clustering results based on Seurat clustering method 
PlotDimension(example_object,  reduction = "umap",idents = "Seurat_r_0.8_k_20")
# cell clustering results based on MCL clustering method 
PlotDimension(example_object, reduction = "umap",idents = "MC_Label")

## ----cell_clustering_globalmarker, eval=TRUE, echo =TRUE----------------------
global_marker <- FindGlobalMarkers(example_object,idents = "Seurat_r_0.8_k_20")
PlotMarkerHeatmap(Globalmarkers = global_marker,object = example_object,idents ="Seurat_r_0.8_k_20")

## -----------------------------------------------------------------------------
sessionInfo()

