#' @include generics.R
#' @include Classes.R
NULL

#' @title plot heatmap based on bicluster
#'
#'
#' @param object Input IRISFGM object
#' @param show.overlap Parameter (logic) indicates whether the figure shows the overlap part between two selected biclusters.
#' @param N.bicluster Number of biclusters.
#' @param show.annotation Parameter (logic) indicates whether to show annotation (biclusters number and cell cluster labels).
#' @param show.clusters Parameter (logic) indicates whether to show the cell cluster label.
#' @return It will generate a heatmap based on selected two FGMs.
#' @name PlotHeatmap
#' @import pheatmap
#' @import Polychrome RColorBrewer
#' @examples 
#' data(example_object)
#' PlotHeatmap(example_object,
#' N.bicluster = c(1,20),
#' show.annotation = TRUE, 
#' show.cluster =TRUE)
#' 
.plotHeatmap <- function(object = object, 
                         N.bicluster = c(1, 5), 
                         show.overlap = FALSE, 
                         show.annotation = FALSE, 
                         show.clusters = FALSE) {
    vec.boolean <- vector(mode = "logical")
    for (i in seq_along(N.bicluster)) {
        vec.boolean[i] <- is.double(N.bicluster[i])
    }
    if (!all(vec.boolean)) {
        stop("please type two bicluster numbers")
    }
    if (length(N.bicluster) != 2) {
        stop("Only plot two bicluster; please type in two numbers")
    }
    condition.index <- N.bicluster
    gene.sub <- c()
    cell.sub <- c()
    for (i in condition.index) {
        gene.sub <- rbind(gene.sub, object@BiCluster@CoReg_gene[object@BiCluster@CoReg_gene$Condition == i, ])
        cell.sub <- rbind(cell.sub, object@BiCluster@CoCond_cell[object@BiCluster@CoCond_cell$Condition == i, ])
    }
    cell.bicluster.1 <- cell.sub$cell_name[cell.sub$Condition == N.bicluster[1]]
    cell.bicluster.2 <- cell.sub$cell_name[cell.sub$Condition == N.bicluster[2]]
    cell.bicluster1.diff <- setdiff(cell.bicluster.1, cell.bicluster.2)
    cell.bicluster2.diff <- setdiff(cell.bicluster.2, cell.bicluster.1)
    cell.overlap <- intersect(cell.bicluster.1, cell.bicluster.2)
    cell.bicluster.vec <- c(rep(N.bicluster[1], 
                                length(cell.bicluster1.diff)), 
                            rep("overlap", length(cell.overlap)), 
                            rep(N.bicluster[2], 
                                length(cell.bicluster2.diff)))
    annotation_col <- data.frame(row.names = c(cell.bicluster1.diff, 
                                               cell.overlap, 
                                               cell.bicluster2.diff), 
                                 bicluster_cell = cell.bicluster.vec)
    if (show.clusters == TRUE) {
        unwanted_lab <- which(colnames(object@MetaInfo) %in% c("Original","ncount_RNA","nFeature"))
        annotation_col <- cbind(annotation_col, 
                                object@MetaInfo[rownames(annotation_col), 
                                                c(-unwanted_lab)])
        colnames(annotation_col)[2:ncol(annotation_col)] <- colnames(object@MetaInfo)[c(-unwanted_lab)]
    }
    gene.bicluster.1 <- gene.sub[, 1][gene.sub$Condition == N.bicluster[1]]
    gene.bicluster.2 <- gene.sub[, 1][gene.sub$Condition == N.bicluster[2]]
    gene.bicluster1.diff <- setdiff(gene.bicluster.1, gene.bicluster.2)
    gene.bicluster2.diff <- setdiff(gene.bicluster.2, gene.bicluster.1)
    gene.overlap <- intersect(gene.bicluster.1, gene.bicluster.2)
    gene.bicluster.vec <- c(rep(N.bicluster[1], 
                                length(gene.bicluster1.diff)), 
                            rep("overlap", length(gene.overlap)), 
                            rep(N.bicluster[2], 
                                length(gene.bicluster2.diff)))
    annotation_row <- data.frame(row.names = c(gene.bicluster1.diff, 
                                               gene.overlap, 
                                               gene.bicluster2.diff), 
                                 bicluster_gene = gene.bicluster.vec)
    if (show.overlap == TRUE) {
        heatmap.matrix <- object@Processed_count[rownames(annotation_row), rownames(annotation_col)]
    } else {
        gene.bicluster.vec <- c(rep(N.bicluster[1], 
                                    length(gene.bicluster1.diff)), 
                                rep(N.bicluster[2], 
                                    length(gene.bicluster2.diff)))
        cell.bicluster.vec <- c(rep(N.bicluster[1], 
                                    length(cell.bicluster1.diff)), 
                                rep(N.bicluster[2], 
                                    length(cell.bicluster2.diff)))
        annotation_row <- data.frame(row.names = c(gene.bicluster1.diff, 
                                                   gene.bicluster2.diff), 
                                     bicluster_gene = as.factor(gene.bicluster.vec))
        annotation_col <- data.frame(row.names = c(cell.bicluster1.diff, 
                                                   cell.bicluster2.diff), 
                                     bicluster_cell = as.factor(cell.bicluster.vec))
        if (show.clusters == TRUE) {
            unwanted_lab <- which(colnames(object@MetaInfo) %in% c("Original","ncount_RNA","nFeature"))
            annotation_col <- cbind(annotation_col, object@MetaInfo[rownames(annotation_col), c(-unwanted_lab)])
            colnames(annotation_col)[2:ncol(annotation_col)] <- colnames(object@MetaInfo)[c(-unwanted_lab)]
        }
        heatmap.matrix <- object@Processed_count[rownames(annotation_row), rownames(annotation_col)]
        
    }
    for (i in 1:ncol(annotation_col)) {
        annotation_col[, i] <- as.factor(annotation_col[, i])
    }
    ann_colors <- list()
    for (i in 1:ncol(annotation_col)) {
        if (colnames(annotation_col)[i] == "bicluster_cell"){
            tmp.col.name <- as.character(unique(annotation_col[, i]))
            if(length(tmp.col.name) == 2 ){
                tmp.color <- brewer.pal(n = 3, name = "Dark2")[1:2]
            } else if (length(tmp.col.name) == 3){
                tmp.color <- brewer.pal(n = length(tmp.col.name), name = "Dark2")  
            }  
        } else {
            tmp.col.name <- as.character(unique(annotation_col[, i]))
            if (length(tmp.col.name) <= 8 ){
                tmp.color <- brewer.pal(n = length(tmp.col.name), name = "Dark2")  
            }  else {
                tmp.color <- Polychrome::palette36.colors(n =length(tmp.col.name))
            } 
        }
       
        names(tmp.color) <- tmp.col.name
        ann_colors[[i]] <- tmp.color
    }
    names(ann_colors) <- colnames(annotation_col)
    ann_colors[["bicluster_gene"]] <- ann_colors$bicluster_cell[unique(annotation_row$bicluster_gene)]
    if (show.annotation == FALSE) {
        annotation_col <- NULL
        annotation_row <- NULL
        ann_colors <- NULL
    }
    if(max(heatmap.matrix) > 100){
        heatmap.matrix <- log1p(heatmap.matrix)
    }
    pheatmap(heatmap.matrix, 
             color = colorRampPalette(c("#0018FF", "#F6F9BE", "#FF6106"))(100), 
             scale = "row", 
             border_color = NA, 
             cluster_rows = FALSE, 
             cluster_cols = FALSE, 
             show_rownames = FALSE, 
             show_colnames = FALSE, 
             fontsize_number = 1, 
             main = paste0("bicluster ", 
                           N.bicluster[1], 
                           " and bicluster ", 
                           N.bicluster[2]), 
             annotation_col = annotation_col, 
             annotation_row = annotation_row, 
             annotation_colors = ann_colors)
}

#' @export
#' @rdname PlotHeatmap
setMethod("PlotHeatmap", "IRISFGM", .plotHeatmap)


