
#' CreateIRISFGMObject
#' @description Create IRIS-FGM object
#' @param x Input expression matrix which should be a matrix or dataframe.
#' @param min.cell each gene should be expressed by at least this many cell.
#' @param min.gene each cell should express this many gene at least.
#' @param LTMGr Automatically create LTMG object.
#' @param Bicluster Automatically create Bicluster object.
#'
#' @return it should return a IRISFGM object of which structure can be found in tutorial.
#' @export
#'
#' @examples \dontrun{object <- CreateIRISCEMObject(x = input_matrix, min.cell = 0, min.gene =0}
globalVariables(c("input_matrix"))
CreateIRISFGMObject <- function(x = input_matrix, min.cell = 0, 
                                min.gene = 0, 
                                LTMGr = new(Class = "LTMGr"), 
                                Bicluster = new(Class = "Bicluster")) {
    raw.matrix <- as.matrix(x)
    raw.matrix.filterbycell <- raw.matrix[(rowSums(raw.matrix > 0) > min.cell), ]
    raw.matrix.filterbygene <- raw.matrix.filterbycell[, (colSums(raw.matrix.filterbycell > 0) > min.gene)]
    message("Creating IRISCEM object. \n", "The original input file contains ", dim(raw.matrix)[2], " cells and ", dim(raw.matrix)[1], " genes \n", "Removed ", 
            dim(raw.matrix)[1] - dim(raw.matrix.filterbycell)[1], " genes that total expression value is equal or less than ", min.cell, "\n", "Removed ", dim(raw.matrix.filterbycell)[2] - 
                dim(raw.matrix.filterbygene)[2], " cells that number of expressed gene is equal or less than ", min.gene)
    my_Object <- new(Class = "IRISFGM", Raw_count = raw.matrix.filterbygene)
    my_Object <- suppressMessages(AddMeta(my_Object))
    return(my_Object)
}
