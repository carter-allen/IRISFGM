#' @include generics.R
#' @include Classes.R
NULL

#' @title plot heatmap based on bicluster
#'
#'
#' @param object IRISFGM object
#' @param show.overlap parameter indicates whether figure show the overlap part betwween two selected biclusters.
#' @param N.bicluster number of biclusters
#'
#' @name PlotHeatmap
#' @import pheatmap
#' @import Polychrome RColorBrewer
.plotHeatmap <- function(object = object, N.bicluster = c(1, 5), show.overlap = F, seed = 123, show.annotation = F, show.clusters = F) {
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
    cell.bicluster.vec <- c(rep(N.bicluster[1], length(cell.bicluster1.diff)), rep("overlap", length(cell.overlap)), rep(N.bicluster[2], length(cell.bicluster2.diff)))
    annotation_col <- data.frame(row.names = c(cell.bicluster1.diff, cell.overlap, cell.bicluster2.diff), bicluster_cell = cell.bicluster.vec)
    if (show.clusters == T) {
        annotation_col <- cbind(annotation_col, object@MetaInfo[rownames(annotation_col), c(-1, -2)])
    }
    gene.bicluster.1 <- gene.sub[, 1][gene.sub$Condition == N.bicluster[1]]
    gene.bicluster.2 <- gene.sub[, 1][gene.sub$Condition == N.bicluster[2]]
    gene.bicluster1.diff <- setdiff(gene.bicluster.1, gene.bicluster.2)
    gene.bicluster2.diff <- setdiff(gene.bicluster.2, gene.bicluster.1)
    gene.overlap <- intersect(gene.bicluster.1, gene.bicluster.2)
    gene.bicluster.vec <- c(rep(N.bicluster[1], length(gene.bicluster1.diff)), rep("overlap", length(gene.overlap)), rep(N.bicluster[2], length(gene.bicluster2.diff)))
    annotation_row <- data.frame(row.names = c(gene.bicluster1.diff, gene.overlap, gene.bicluster2.diff), bicluster_gene = gene.bicluster.vec)
    if (show.overlap == T) {
        heatmap.matrix <- object@Processed_count[rownames(annotation_row), rownames(annotation_col)]
    } else {
        gene.bicluster.vec <- c(rep(N.bicluster[1], length(gene.bicluster1.diff)), rep(N.bicluster[2], length(gene.bicluster2.diff)))
        cell.bicluster.vec <- c(rep(N.bicluster[1], length(cell.bicluster1.diff)), rep(N.bicluster[2], length(cell.bicluster2.diff)))
        annotation_row <- data.frame(row.names = c(gene.bicluster1.diff, gene.bicluster2.diff), bicluster_gene = as.factor(gene.bicluster.vec))
        annotation_col <- data.frame(row.names = c(cell.bicluster1.diff, cell.bicluster2.diff), bicluster_cell = as.factor(cell.bicluster.vec))
        if (show.clusters == T) {
            annotation_col <- cbind(annotation_col, object@MetaInfo[rownames(annotation_col), c(-1, -2)])
        }
        heatmap.matrix <- object@Processed_count[rownames(annotation_row), rownames(annotation_col)]
        
    }
    for (i in 1:ncol(annotation_col)) {
        annotation_col[, i] <- as.factor(annotation_col[, i])
    }
    ann_colors <- list()
    for (i in 1:ncol(annotation_col)) {
        tmp.col.name <- as.character(unique(annotation_col[, i]))
        tmp.color <- brewer.pal(n = length(tmp.col.name), name = "Dark2")
        names(tmp.color) <- tmp.col.name
        ann_colors[[i]] <- tmp.color
    }
    names(ann_colors) <- colnames(annotation_col)
    ann_colors[["bicluster_gene"]] <- ann_colors$bicluster_cell[unique(annotation_row$bicluster_gene)]
    if (show.annotation == F) {
        annotation_col <- NULL
        annotation_row <- NULL
        ann_colors <- NULL
    }
    pheatmap(heatmap.matrix, color = colorRampPalette(c("#0018FF", "#F6F9BE", "#FF6106"))(100), scale = "row", border_color = NA, cluster_rows = F, cluster_cols = F, 
        show_rownames = F, show_colnames = F, fontsize_number = 1, main = paste0("bicluster ", N.bicluster[1], " and bicluster ", N.bicluster[2]), annotation_col = annotation_col, 
        annotation_row = annotation_row, annotation_colors = ann_colors)
}

#' @export
#' @rdname PlotHeatmap
setMethod("PlotHeatmap", "IRISFGM", .plotHeatmap)


