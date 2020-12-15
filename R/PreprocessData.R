#' @include generics.R
NULL

#' Process data
#'
#' @description  Process data via normalization and imputation.
#' @param object
#'
#' @param normalization two options: (1)library size noralization by using library size factor: 1e6, equal to CPM, or (2) using 'scran' normalization method.
#' @param seed set seed for reproducibility.
#' @param IsImputation imputation method is provided by DrImpute. Default is FALSE.
#'
#' @name ProcessData
#'
#' @importFrom scater normalize logNormCounts
#' @importFrom SingleCellExperiment SingleCellExperiment normcounts
#' @importFrom scran quickCluster computeSumFactors
#' @importFrom Seurat as.sparse
#' @importFrom DrImpute DrImpute

.processData <- function(object = NULL, normalization = "cpm", library.size = 1e+05, IsImputation = FALSE, seed = 123) {
    if (is.null(object@MetaInfo)) {
        stop("Can not find meta data, please run AddMeta")
    }
    Input <- object@Raw_count[, rownames(object@MetaInfo)]
    set.seed(seed)
    random.number <- sample(c(1:nrow(Input)), 100)
    if (all(as.numeric(unlist(Input[random.number, ]))%%1 == 0)) {
        ## normalization##############################
        if (grepl("cpm", ignore.case = T, normalization)) {
            my.normalized.data <- (Input/colSums(Input)) * library.size
        } else if (grepl("scran", ignore.case = T, normalization)) {
            sce <- SingleCellExperiment(assays = list(counts = Input))
            clusters <- quickCluster(sce, min.size = floor(ncol(Input)/3))
            sce <- computeSumFactors(sce, clusters = clusters)
            sce <- scater::logNormCounts(sce, log = FALSE)
            my.normalized.data <- normcounts(sce)
        }
    } else {
        my.normalized.data <- Input
        
    }
    
    ## imputation#################################
    if (IsImputation == TRUE) {
        my.imputated.data <- DrImpute(as.matrix(my.normalized.data), dists = "spearman")
    } else {
        my.imputated.data <- my.normalized.data
    }
    colnames(my.imputated.data) <- colnames(Input)
    rownames(my.imputated.data) <- rownames(Input)
    object@Processed_count <- my.imputated.data
    return(object)
}


#' @export
#' @rdname ProcessData
setMethod("ProcessData", "IRISFGM", .processData)





