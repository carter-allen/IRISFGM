#' Create DimReduce object
#'
#' @slot PCA ANY.
#' @slot UMAP ANY.
#' @slot TSNE ANY.
setClass("DimReduce", slots = c(PCA = "ANY", UMAP = "ANY", TSNE = "ANY"))
# create Bicluster object
#'
#' @slot block ANY.
#' @slot CoReg_gene ANY.
#' @slot CoCond_cell ANY.
setClass("Bicluster", slots = c(CoReg_gene = "ANY", CoCond_cell = "ANY", MarkerGene = "ANY", PathwayFromMC = "ANY", PathwayFromModule = "ANY"))
# create LTMG object
#'
#' @slot LTMG_discrete matrix.
#' @slot LTMG_BinarySingleSignal matrix.
#' @slot DimReduce DimReduce
#' @slot Cluster ANY
#' @slot MarkerGene ANY
#' @slot Pathway ANY
#' @slot LTMG_BinaryMultisignal matrix.
setClass("LTMGr", slots = c(LTMG_discrete = "matrix", LTMG_BinarySingleSignal = "matrix", LTMG_BinaryMultisignal = "matrix", DimReduce = "DimReduce", MarkerGene = "ANY", 
                            Pathway = "ANY", Tmp.seurat = "ANY"))
# set IRISFGM class
#' @title IRISFGM class
#'
#' @slot MetaInfo ANY.
#' @slot Discretization matrix.
#' @slot Raw_count ANY
#' @slot Processed_count ANY
#' @slot BiCluster Bicluster
#' @slot LTMG LTMGr.
#'
#' @return IRIS-FGM object
#' @name IRISFGM
#' @rdname IRISFGM
#' @exportClass IRISFGM
setClass("IRISFGM", slots = c(Raw_count = "ANY", 
                              Processed_count = "ANY", 
                              MetaInfo = "ANY", 
                              Discretization = "matrix", 
                              LTMG = "LTMGr", 
                              BiCluster = "Bicluster"))

setMethod(f = "show", signature = "IRISFGM", definition = function(object) {
  cat("The data contain", "\n", nrow(object@Raw_count), "genes.\n", ncol(object@Raw_count), "cells.")
  invisible(x = NULL)
})

