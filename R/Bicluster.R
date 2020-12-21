#' @include generics.R
#' @include Classes.R
NULL

#' getblock function
#' 
#' get block from generated temporal files
#'
#' @param object IRISFGM object
#' @param keyword 'Conds' for co-regulatory or 'Genes' for co-expression gene
.getBlock <- function(object = NULL, keyword = "Conds") {
    tmp.block <- readLines(paste0(getwd(), "/tmp_expression.txt.chars.blocks"))
    tmp.bc <- grep(keyword, tmp.block, value = T)
    tmp.cel.module <- sapply(strsplit(tmp.bc, ":", 2), "[", 2)
    CONDS <- as.character()  # store the conditions
    label_C <- as.numeric()  # store the occurence of one condistions
    
    for (j in 1:length(tmp.cel.module)) {
        BCcond <- unlist(strsplit(tmp.cel.module[j], split = " "))
        BCcond <- BCcond[BCcond != ""]  # exclude the blank string
        CONDS <- c(CONDS, BCcond)
        label_C <- c(label_C, rep(j, length(BCcond)))
    }
    df_C <- data.frame(cell_name = CONDS, Condition = label_C)
    if (keyword == "Conds") {
        df_C$cell_name <- as.character(df_C$cell_name)
        object@BiCluster@CoCond_cell <- df_C
        return(object)
    } else if (keyword == "Genes") {
        tmp.df_C <- df_C
        tmp.gene.list <- as.character(tmp.df_C$cell_name)
        tmp.gene.list <- gsub("_[0-9]$", "", tmp.gene.list)
        tmp.df_C$cell_name <- tmp.gene.list
        object@BiCluster@CoReg_gene <- tmp.df_C
        colnames(object@BiCluster@CoReg_gene) <- c("Gene", "Condition")
        return(object)
    }
    
}

#' 
#' @name RunDiscretization
#' @title RunDiscretization
#'
#' @description Run discretization based on Quantile method
#'
#' @param object input IRIS-FGM object
#' @param q quantile number which is used as discretized cutoff. The bigger q means more cells will be categorized into 1 in terms of binarizing one gene.
#' 
#' 
#' @return It will generate quantile based binary matrix.
.runDiscretization <- function(object = NULL, q = 0.06) {
    message("writing temporary expression file ...")
    tmp.dir <- paste0(getwd(), "/tmp_expression.txt")
    tmp.count <- object@Processed_count
    tmp.count <- cbind(ID = rownames(tmp.count), tmp.count)
    write.table(tmp.count, file = tmp.dir, row.names = F, quote = F, sep = "\t")
    message("writing quantile discretization file ...")
    qubic(i = tmp.dir, Fa = TRUE, q = q, R = FALSE)
    tmp.chars <- paste0(getwd(), "/tmp_expression.txt.chars")
    tmp.readin <- read.table(tmp.chars, row.names = 1, header = T)
    object@Discretization <- as.matrix(tmp.readin)
    return(object)
}


#' @export
#' @rdname RunDiscretization
setMethod("RunDiscretization", "IRISFGM", .runDiscretization)

#' RunBicusterBaseOnLTMG
#'
#' Run bicluster based on LTMG
#'
#' @param object input IRISFMG object
#' @param OpenDual the flag using the lower bound of condition number. Default: 5 percent of the gene number in current bicluster.
#' @param Extension consistency level of the block (0.5-1.0], the minimum ratio between the number of identical valid symbols in a column and the total number of rows in the output. Default: 1.0.
#' @param NumBlockOutput number of blocks to report. Default: 100.
#' @param BlockOverlap filtering overlapping blocks. Default: 0.7.
#' @param BlockCellMin minimum column width of the block. Default: 15 columns.
.runBiclusterBaseOnLTMG <- function(object = NULL, OpenDual = FALSE, Extension = 1, NumBlockOutput = 100, BlockOverlap = 0.9, BlockCellMin = 15) {
    print("writing LTMG discretization file ...")
    tmp.dir <- paste0(getwd(), "/tmp_expression.txt.chars")
    tmp.multi <- object@LTMG@LTMG_BinaryMultisignal
    tmp.multi <- cbind(ID = rownames(tmp.multi), tmp.multi)
    write.table(tmp.multi, file = tmp.dir, row.names = F, quote = F, sep = "\t")
    print("finished!")
    print("running biclustering . . .")
    qubic(i = tmp.dir, d = TRUE, C = OpenDual, c = Extension, o = NumBlockOutput, f = BlockOverlap, k = BlockCellMin)
}

#' Run Discretization
#' Generate temporal discretized file
#'
#' @param object input IRISFMG object
#' @param OpenDual the flag using the lower bound of condition number. Default: 5 percent of the gene number in current bicluster.
#' @param Extension consistency level of the block (0.5-1.0], the minimum ratio between the number of identical valid symbols in a column and the total number of rows in the output. Default: 1.0.
#' @param NumBlockOutput number of blocks to report. Default: 100.
#' @param BlockOverlap filtering overlapping blocks. Default: 0.7.
#' @param BlockCellMin minimum column width of the block. Default: 15 columns.
.runBiclusterBaseOnDiscretization <- function(object = NULL, OpenDual = TRUE, Extension = 1, NumBlockOutput = 100, BlockOverlap = 1, BlockCellMin = 15) {
    tmp.dir <- paste0(getwd(), "/tmp_expression.txt.chars")
    if (file.exists(tmp.dir)) {
        qubic(i = tmp.dir, d = TRUE, C = OpenDual, c = Extension, o = NumBlockOutput, f = BlockOverlap, k = BlockCellMin)
    } else {
        print("please use `RunDiscretization` first and then execute this command")
    }
    
}

#' Run cluster
#' 
#' This function will identify the Biclusters based on LTMG or Quantile normalization
#' 
#' @param object input IRIS-FGM object
#' @param DiscretizationModel use different discretization method, including 'Quantile' and 'LTMG.'
#' @param OpenDual the flag using the lower bound of condition number. Default: 5 percent of the gene number in current bicluster.
#' @param NumBlockOutput number of blocks to report. Default: 100.
#' @param BlockOverlap filtering overlapping blocks. Default: 0.7.
#' @param Extension 
#' @param BlockCellMin minimum column width of the block. Default: 15 columns.
#'
#' @name RunBicluster
#' @return It will generate a temporal file on local directory for processed data named 'tmp_expression.txt', discretized file named 
#' 'tmp_expression.txt.chars', and biclsuter block named 'tmp_expression.txt.chars.block'.
#' @examples
#' # based on LTMG discretization
#' \dontrun{
#' object <- RunBicluster(object, 
#' DiscretizationModel = 'LTMG',
#' OpenDual = F,
#' NumBlockOutput = 1000, 
#' BlockOverlap = 0.7, 
#' BlockCellMin = 15)
#' }
#' # based on quantile discretization
#' \dontrun{
#' object <- RunBicluster(object, 
#' DiscretizationModel = 'Quantile',
#' OpenDual = F, 
#' NumBlockOutput = 1000, 
#' BlockOverlap = 0.7, 
#' BlockCellMin = 15)
#' }
.runBicluster <- function(object = NULL, DiscretizationModel = "Quantile", OpenDual = FALSE, Extension = 0.9, NumBlockOutput = 100, BlockOverlap = 0.7, BlockCellMin = 15) {
    if (DiscretizationModel != "LTMG" && DiscretizationModel != "Quantile") {
        stop("please select either LTMG or Quantile")
    }
    if (DiscretizationModel == "LTMG") {
        .runBiclusterBaseOnLTMG(object = object, OpenDual = OpenDual, Extension = Extension, NumBlockOutput = NumBlockOutput, BlockOverlap = BlockOverlap, 
                                BlockCellMin = BlockCellMin)
    }
    if (DiscretizationModel == "Quantile") {
        .runBiclusterBaseOnDiscretization(object = object, OpenDual = OpenDual, Extension = Extension, NumBlockOutput = NumBlockOutput, BlockOverlap = BlockOverlap, 
                                          BlockCellMin = BlockCellMin)
    }
    object <- .getBlock(object = object, keyword = "Conds")
    object <- .getBlock(object = object, keyword = "Genes")
    return(object)
    
}
#' @export
#' @rdname RunBicluster
setMethod("RunBicluster", "IRISFGM", .runBicluster)












