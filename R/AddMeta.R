#' @include generics.R
#' @include Classes.R
NULL

#' Create addmeta function for adding meta information
#' 
#' This function can import cell annotation information to IRIS-FGM object.
#' @param object input IRIS-FGM object
#' @param meta.info meta information table should be a data frame with rows representing cell and coloumn representing different group condition
#'
#' @return It will add meta informatoin to IRISFGM.
#' @name AddMeta
#' @examples \dontrun{
#' AddMeta(object, 
#' meta.info = my.meta)
#' }  
.addMeta <- function(object = NULL, meta.info = NULL) {
    if (is.null(meta.info)) {
        message("Do not provide meta info table for the object.")
        message("using original cell identity")
        meta.info <- data.frame(row.names = as.character(colnames(object@Raw_count)), 
                                Original = as.character(colnames(object@Raw_count)), 
                                ncount_RNA = as.numeric(colSums(object@Raw_count)), 
                                nFeature = as.numeric(colSums(object@Raw_count > 0)))
        object@MetaInfo <- meta.info
    } else {
        if (!identical(colnames(object@Raw_count),rownames(meta.info))) {
            stop("\n There is inconsisten cell names between meta info and original raw count.
           \n Please check rownames of meta info and colnames of original raw count.")
        } else {
            tmp.meta <- cbind(ncount_RNA = as.numeric(colSums(object@Raw_count)), nFeature = as.numeric(colSums(object@Raw_count > 0)), meta.info)
            object@MetaInfo <- tmp.meta
        }
    }
    return(object)
}

#' @rdname AddMeta
#' @export
setMethod("AddMeta", "IRISFGM", .addMeta)
globalVariables(c("RNA_count.name","RNA_count","Feature_number.name","Feature_number"))
#' Generate violin plot based on meta data.
#' 
#' This function can plot figure based on numebr total count information and this step is for the quality control. we shoud exclude extreme value in data.
#' @import ggplot2
#' @importFrom  ggpubr ggarrange
#' @param object input IRIS-FGM object
#' @name PlotMeta
#' @return It will generate two violin plots regarding number of RNA count and identified gene number.
#'
#' @examples \dontrun{PlotMeta(object)}
.plotMeta <- function(object = NULL) {
    if (is.null(object@MetaInfo)) {
        stop("Can not find meta data, please run AddMeta")
    }
    my.data <- data.frame(row.names = rownames(object@MetaInfo), RNA_count = object@MetaInfo$ncount_RNA, Feature_number = object@MetaInfo$nFeature, RNA_count.name = rep("RNA_count", 
                                                                                                                                                                         nrow(object@MetaInfo)), Feature_number.name = rep("Feature_number", nrow(object@MetaInfo)))
    p.1 <- ggplot(data = my.data, aes(x = RNA_count.name, y = RNA_count)) + geom_violin(trim = FALSE, color = "#E69F00", fill = "#E69F00") + geom_jitter(shape = 16) + 
        labs(x = NULL) + ggtitle("RNA_count")
    p.2 <- ggplot(data = my.data, aes(x = Feature_number.name, y = Feature_number)) + geom_violin(trim = FALSE, color = "#E69F00", fill = "#E69F00") + geom_jitter(shape = 16) + 
        labs(x = NULL) + ggtitle("Feature_number")
    combine.p <- ggpubr::ggarrange(p.1, p.2, ncol = 2)
    print(combine.p)
}

#' @rdname PlotMeta
#' @export
setMethod("PlotMeta", "IRISFGM", .plotMeta)

#' Subset data by number of count and numebr of gene.
#' 
#' This function is used for filtering out low-quality data based on previous result generated from \code{\link{PlotMeta}}
#'
#' @param object input IRIS-FGM object
#' @param nFeature.upper select upper limit for number of feature
#' @param nFeature.lower select lower limit for number of feature
#' @param Counts.upper select upper limit for number of UMI counts
#' @param Counts.lower select lower limit for number of UMI counts
#' @name SubsetData 
#' @return it will filter out some cell regarding threshhold.
#'
#' @examples # Use Yan's data which posts on github tutorial
#' \dontrun{ 
#' object <- SubsetData(object, 
#'     nFeature.upper=15000,
#'     nFeature.lower=8000, 
#'     Counts.upper=700000,
#'     Counts.lower=400000)}
#' 
.subset_data <- function(object, nFeature.upper = Inf, nFeature.lower = -Inf, Counts.upper = Inf, Counts.lower = -Inf) {
    Meta <- object@MetaInfo
    if (is.null(Meta)) {
        stop("Can not find meta data, please run AddMeta")
    } else {
        Meta <- Meta[Meta$nFeature < nFeature.upper & Meta$nFeature > nFeature.lower & Meta$ncount_RNA < Counts.upper & Meta$ncount_RNA > Counts.lower, ]
    }
    object@MetaInfo <- Meta
    return(object)
}

#' @rdname SubsetData
#' @export
setMethod("SubsetData", "IRISFGM", .subset_data)
