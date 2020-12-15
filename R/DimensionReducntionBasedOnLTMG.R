#' @include generics.R
#' @include Classes.R
NULL

#' Run dimension reduction based on seurat method 
#' 
#' This function is based on the Seurat package to perform dimension reduction. The input matirx is the LTMG signaling matrix.
#' @param object input IRIS-FGM
#'
#' @param reduction select a method for dimension reduction, including umap, tsne, and pca.
#' @param dims select the number of PCs from PCA results to perform the following dimension reduction and cell clustering.
#' @param mat.source choose source data for running this function either from LTMG signal matrix or from processed data. Values of this parameter are 'LTMG' and 'UMImatrix' 
#' @name RundimensionReduction
#' @importFrom Seurat CreateSeuratObject ScaleData RunPCA RunTSNE RunUMAP FindVariableFeatures
#' @return This function will generate pca, tsne, or umap dimension reduction results.
#' @examples \dontrun{
#' obejct <- RunDimensionReduction(object,
#'    mat.source= 'LTMG',
#'    reduction = 'umap', 
#'    dims = 1:15 ,
#'    perplexity = 15, 
#'    seed = 1)}
.runDimensionReduction <- function(object, mat.source = c("LTMG", "UMImatrix"), reduction = "umap", dims = 1:15, perplexity = 15, seed = 1) {
    if (mat.source == "LTMG") {
        Tmp.seurat <- CreateSeuratObject(object@LTMG@LTMG_discrete)
        Tmp.seurat <- ScaleData(Tmp.seurat)
    } else if (mat.source == "UMImatrix") {
        Tmp.seurat <- CreateSeuratObject(object@Processed_count)
        Tmp.seurat <- FindVariableFeatures(Tmp.seurat, selection.method = "vst", nfeatures = 2000)
        Tmp.seurat <- ScaleData(Tmp.seurat)
    } else {
        stop("please select either LTMG or UMImatrix as input matrix")
    }
    Tmp.seurat <- suppressMessages(RunPCA(Tmp.seurat, features = rownames(Tmp.seurat@assays$RNA)))
    object@LTMG@DimReduce@PCA <- Tmp.seurat@reductions$pca@cell.embeddings
    if (grepl("tsne", reduction, ignore.case = T) || grepl("umap", reduction, ignore.case = T)) {
        if (grepl("tsne", reduction, ignore.case = T)) {
            Tmp.seurat <- RunTSNE(Tmp.seurat, dims = dims, seed.use = seed, perplexity = perplexity)
            object@LTMG@DimReduce@TSNE <- Tmp.seurat@reductions$tsne@cell.embeddings
        }
        if (grepl("umap", reduction, ignore.case = T)) {
            Tmp.seurat <- suppressMessages(RunUMAP(Tmp.seurat, dims = dims))
            object@LTMG@DimReduce@UMAP <- Tmp.seurat@reductions$umap@cell.embeddings
        }
    } else {
        stop("choose a dimension reduction method between umap or tsne")
    }
    object@LTMG@Tmp.seurat <- Tmp.seurat
    return(object)
}


#' @rdname RunDimensionRecution
#' @export
setMethod("RunDimensionReduction", "IRISFGM", .runDimensionReduction)



#' Classify cell type prediction
#' 
#' This function is based on Seurat package.
#' @param object input IRIS-FGM object.
#' @param k.param Defines k for the k-nearest neighbor algorithm.
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param algorithm Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). Leiden requires the leidenalg python.
#' @param dims \tDimensions of reduction to use as input.
#' @name RunClassification
#' @importFrom Seurat FindNeighbors FindClusters
#' @import ggplot2
#' @return It will generate cell type inforamtion.
#' @examples \dontrun{object <- RunClassification(object, dims = 1:15, k.param = 20, resolution = 0.6, algorithm = 1)}
.runClassification <- function(object, dims = 1:15, k.param = 20, resolution = 0.6, algorithm = 1) {
    if (is.null(object@LTMG@Tmp.seurat)) {
        stop("There is no temporary seurat obejct getting detected. \n Try to run RundimensionRuduction first.")
    }
    Tmp.seurat <- object@LTMG@Tmp.seurat
    Tmp.seurat <- FindNeighbors(Tmp.seurat, dims = dims, k.param = k.param)
    Tmp.seurat <- FindClusters(Tmp.seurat, resolution = resolution, algorithm = algorithm)
    tmp.meta <- object@MetaInfo
    tmp.colname <- colnames(tmp.meta)
    res.index <- grep(paste0("res.", resolution), colnames(Tmp.seurat@meta.data))
    tmp.colname <- c(tmp.colname, paste0("Seurat", resolution))
    tmp.meta <- cbind(tmp.meta, Tmp.seurat@meta.data[, res.index])
    colnames(tmp.meta) <- tmp.colname
    object@MetaInfo <- tmp.meta
    object@LTMG@Tmp.seurat <- Tmp.seurat
    return(object)
}

#' @rdname RunDimensionRecution
#' @export
setMethod("RunClassification", "IRISFGM", .runClassification)




#' Visualize dimension reduction results
#'
#' Generate Umap and it requires user to input cell label index on console window.
#' @param object Input IRIS-FGM Object
#' @param reduction Choose one of approaches for dimension reduction, including 'pca', 'tsne', 'umap'.
#' @param pt_size Point size, default is 0.
#'
#' @return generate plot on umap space.
#' @export
#' @name PlotDimension
#' @examples \dontrun{PlotDimension(object)}
.plotDimension <- function(object, reduction = "umap", pt_size = 1) {
    
    if (grepl("tsne", reduction, ignore.case = T) || grepl("umap", reduction, ignore.case = T) || grepl("pca", reduction, ignore.case = T)) {
        
        if (grepl("tsne", reduction, ignore.case = T)) {
            tmp.plot.table <- object@LTMG@DimReduce@TSNE[, c(1, 2)]
        } else if (grepl("umap", reduction, ignore.case = T)) {
            tmp.plot.table <- object@LTMG@DimReduce@UMAP[, c(1, 2)]
        } else if (grepl("tsne", reduction, ignore.case = T)) {
            tmp.plot.table <- object@LTMG@DimReduce@PCA[, c(1, 2)]
        }
    } else {
        stop("choose a dimension reduction method from pca, umap, or tsne")
    }
    message("select condition to present")
    message(paste0(c(1:ncol(object@MetaInfo)), " : ", c(colnames(object@MetaInfo)), "\n"))
    ident.index <- readline(prompt = "select index of cell condition: ")
    ident.index <- as.numeric(ident.index)
    tmp.ident <- object@MetaInfo[, ident.index]
    names(tmp.ident) <- rownames(object@MetaInfo)
    # check name later add
    tmp.plot.table <- cbind.data.frame(tmp.plot.table, Cell_type = as.character(tmp.ident))
    p.cluster <- ggplot(tmp.plot.table, aes(x = tmp.plot.table[, 1], y = tmp.plot.table[, 2], col = tmp.plot.table[, "Cell_type"]))
    
    p.cluster <- p.cluster + geom_point(stroke = pt_size, size = pt_size)
    
    p.cluster <- p.cluster + guides(colour = guide_legend(override.aes = list(size = 5))) + labs(color = colnames(object@MetaInfo)[ident.index]) + xlab("Dimension 1") + 
        ylab("Dimentsion 2")
    p.cluster <- p.cluster + theme_classic() + coord_fixed()
    print(p.cluster)
}

#' @rdname PlotDimension
#' @export
setMethod("PlotDimension", "IRISFGM", .plotDimension)





