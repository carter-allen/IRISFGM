#' @include generics.R
#' @include Classes.R
NULL

#' @param object input IRISFGM object
#'
#' @param SimpleResult marker gene only output log fold change (LFC), p-value, and adjusted p-value.
#' @param FDR a number to specify the threshold of FDR, default by 0.05
#' @name FindMarker
#' @importFrom DEsingle DEsingle DEtype
#'
.findMarker <- function(object, SimpleResult = T, FDR = 0.05) {
    # two group number as factor.
    message("select condition to compare")
    message(paste0(c(1:ncol(object@MetaInfo)), " : ", c(colnames(object@MetaInfo)), "\n"))
    ident.index <- readline(prompt = "select index of cell condition: ")
    ident.index <- as.numeric(ident.index)
    tmp.ident <- object@MetaInfo[, ident.index]
    names(tmp.ident) <- rownames(object@MetaInfo)
    label.used <- colnames(object@MetaInfo)[ident.index]
    # create index table
    tmp.group.table <- data.frame(index = 1:length(unique(tmp.ident)), condition = as.character(sort(unique(tmp.ident))), stringsAsFactors = F)
    tmp.group.table <- rbind(tmp.group.table, c(nrow(tmp.group.table) + 1, "rest of all"))
    # select groups to compare
    message("select index (left) of first group to compare : ")
    message(paste0(tmp.group.table$index[1:nrow(tmp.group.table) - 1], " : ", tmp.group.table$condition[1:nrow(tmp.group.table) - 1], "\n"))
    group.1.idx <- readline("input first group index : ")
    message("select index (left) of second group to compare : ")
    message(paste0(tmp.group.table$index[tmp.group.table$index != group.1.idx], " : ", tmp.group.table$condition[tmp.group.table$index != group.1.idx], "\n"))
    group.2.idx <- readline(prompt = "select index of group 2: ")
    group.1 <- tmp.group.table[tmp.group.table$index == group.1.idx, 2]
    group.2 <- tmp.group.table[tmp.group.table$index == group.2.idx, 2]
    tmp.expression.table <- object@LTMG@LTMG_discrete
    if (group.2 == "rest of all") {
        new.condition <- as.character(tmp.ident)
        # set group in new condition
        new.condition[tmp.ident == group.1] <- group.1
        new.condition[tmp.ident != group.1] <- "all"
        results <- DEsingle(counts = tmp.expression.table, group = as.factor(new.condition))
    } else {
        new.condition <- vector(mode = "logical", length = length(tmp.ident))
        new.condition[tmp.ident == group.1] <- "TRUE"
        new.condition[tmp.ident == group.2] <- "TRUE"
        tmp.new.condition <- as.factor(as.character(tmp.ident)[new.condition == "TRUE"])
        tmp.expression.table <- tmp.expression.table[, new.condition == "TRUE"]
        results <- DEsingle(counts = tmp.expression.table, group = as.factor(tmp.new.condition))
    }
    results.classified <- DEtype(results = results, threshold = FDR)
    if (SimpleResult == TRUE) {
        gene.name <- rownames(results.classified)
        results.classified <- cbind(log(results.classified$foldChange), results.classified$pvalue, results.classified$pvalue.adj.FDR)
        colnames(results.classified) <- c("LFC", "pval", "pvalue.adj.FDR")
        rownames(results.classified) <- gene.name
        results.classified <- as.data.frame(results.classified)
        results.classified <- results.classified[results.classified$pvalue.adj.FDR < FDR, ]
    }
    if (grepl("MC", label.used, ignore.case = T)) {
        object@BiCluster@MarkerGene <- results.classified
    } else {
        object@LTMG@MarkerGene <- results.classified
    }
    
    return(object)
}

#' @rdname FindMarker
#' @export
setMethod("FindMarker", "IRISFGM", .findMarker)

#' This function is for finding global marker
#' FindGlobalMarkers is based on Seurat FindAllMarkers
#'
#' @param object input IRIS-FGM object
#' @param idents choose an idents for labelling cells
#' @param logfc.threshold Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.25 Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' @param test.use same as \code{\link{FindAllMarkers}}
#' @param only.pos
#' @param random.seed
#' @param min.pct
#' @name FindGlobalMarkers
#' @return Output is a differentially expressed gene list
#' @importFrom Seurat Idents FindAllMarkers AddMetaData Idents<-
#'
#' @examples \dontrun{Global_marker <- FindGlobalMarkers(object)}
.findglobalMarkers <- function(object = NULL, idents = "Seurat0.6", logfc.threshold = 0.25, test.use = "wilcox", only.pos = TRUE, random.seed = 1, min.pct = 0.1) {
    tmp.seurat <- object@LTMG@Tmp.seurat
    tmp.seurat <- AddMetaData(tmp.seurat, metadata = object@MetaInfo)
    if (!idents %in% colnames(object@MetaInfo)) {
        stop(paste0("Cannot find idents: ", idents, ".", "\nPlease choose from here:\n", paste0(colnames(object@MetaInfo)[c(-1, -2)], collapse = ",")))
    }
    Idents(tmp.seurat) <- idents
    tmp.all.marker <- FindAllMarkers(tmp.seurat, logfc.threshold = logfc.threshold, test.use = test.use, only.pos = only.pos, random.seed = random.seed, min.pct = min.pct)
    return(tmp.all.marker)
}

#' @rdname FindGlobalMarkers
#' @export
setMethod("FindGlobalMarkers", "IRISFGM", .findglobalMarkers)


#' PlotMarkerHeatmap will visualize global marker
#' This function will generate global marker gene heatmap
#'
#' @param Globalmarkers output from FindGlobalMarkers
#' @param top.gene this number of genes will be used for generating heatmap
#' @param p.adj adjusted pvalue cutoff for gene selection threshold
#' @param idents set current idents
#' @param scale 
#'
#' @return heatmap
#' @name PlotMarkerHeatmap
#' @export
#' @import colorspace
#' @examples 
PlotMarkerHeatmap <- function(Globalmarkers = NULL, idents = "Seurat0.6", top.gene = 50, p.adj = 0.05, scale = "row") {
    marker.list <- Globalmarkers
    marker.list <- marker.list[marker.list$p_val_adj < p.adj, ]
    marker.cluster.index <- as.data.frame(cbind(index = 1:length(unique(marker.list$cluster)), cluster = unique(as.character(marker.list$cluster))), stringsAsFactors = F)
    sub.marker.list <- c()
    for (i in 1:length(unique(marker.list$cluster))) {
        tmp.cluster <- marker.cluster.index$cluster[marker.cluster.index$index == i]
        tmp.marker.list <- marker.list[marker.list$cluster == tmp.cluster, ]
        tmp.marker.list <- tmp.marker.list[1:top.gene, ]
        sub.marker.list <- rbind(sub.marker.list, tmp.marker.list)
    }
    if (!idents %in% colnames(object@MetaInfo)) {
        stop(paste0("Cannot find idents: ", idents, ".", "\nPlease choose from here:\n", paste0(colnames(object@MetaInfo)[c(-1, -2)], collapse = ",")))
    }
    my.cluster <- data.frame(cell = rownames(object@MetaInfo), cluster = object@MetaInfo[, idents])
    my.cluster <- my.cluster[order(my.cluster$cluster), ]
    if (!is.factor(my.cluster$cluster)) {
        col.data <- as.factor(my.cluster$cluster)
    } else {
        col.data <- my.cluster$cluster
    }
    heatmap.matrix <- object@Processed_count[sub.marker.list$gene, as.character(my.cluster$cell)]
    heatmap.matrix <- na.omit(heatmap.matrix)
    colSide <- qualitative_hcl(length(unique(col.data)), palette = "Dynamic")[col.data]
    colMain <- colorRampPalette(c("#2D80BD", "#DCEAF5", "#B0390C"))(100)
    heatmap(as.matrix(heatmap.matrix), Colv = NA, Rowv = NA, scale = scale, ColSideColors = colSide, col = colMain, labCol = "0")
    legend("topright", legend = c(unique(as.character(my.cluster$cluster))), bg = "transparent", pch = rep(19, length(unique(col.data))), col = qualitative_hcl(length(unique(col.data)), 
        palette = "Dynamic")[unique(my.cluster$cluster)], cex = 2, box.lty = 0)
    
}



