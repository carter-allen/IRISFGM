#' @include generics.R
#' @include Classes.R
NULL

#' separateBic
#' separate biclusters 
#'
#' @param object Input IRIS-FGM object.
#'
#' @return It will reture a list for shoring gene and cell. 
#'
.separateBic <- function(object = NULL) {
    tmp.expression <- object@Processed_count
    bic.number <- length(unique(object@BiCluster@CoCond_cell$Condition))
    Bic <- c()
    Bic.name <- c()
    for (i in seq_len(bic.number)) {
        gene.name <- object@BiCluster@CoReg_gene$Gene[object@BiCluster@CoReg_gene$Condition == i]
        cell.name <- object@BiCluster@CoCond_cell$cell_name[object@BiCluster@CoCond_cell$Condition == i]
        tmp.Bic <- tmp.expression[gene.name, cell.name]
        Bic.name <- c(Bic.name, paste0("Bicluster", i))
        Bic <- c(Bic, list(tmp.Bic))
    }
    names(Bic) <- Bic.name
    return(Bic)
}



#' PlotNetwork
#' This function is to plot the network for selected biclusters to show the genes (or cells) overlapping relations.
#'
#' @param object Input IRIS-FGM object.
#' @param edge.by Should be "cell" or by "gene," indicating nodes label. The default value is by "gene."
#' @param lay.out Should be one type of layouts to show nodes' arrangement, including 'linear'(default), 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 
#' 'randomly', 'fr', 'kk', 'drl', 'lgl'.
#' @param N.bicluster Should be two integers indicating the number of two biclusters.
#'
#' @name PlotNetwork
#' @return It will generate a global network regarding overlapping genes or cells.
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom ggraph ggraph geom_node_point geom_edge_arc scale_edge_width geom_node_text theme_graph
#' @import ggplot2
#' @examples \dontrun{object <- PlotNetwork(object,edge.by = 'gene',N.bicluster =c(1:20) )}
globalVariables(c("name", "weight", "status","weight"))
.plotnetwork <- function(object, edge.by = "gene", lay.out = "linear", N.bicluster = seq_len(20)) {
    Bic.list <- .separateBic(object)
    Bic.list.select <- Bic.list[N.bicluster]
    ntwork.adjacency.mtx <- matrix(seq_len(length(N.bicluster) * length(N.bicluster)), nrow = length(N.bicluster))
    rownames(ntwork.adjacency.mtx) <- colnames(ntwork.adjacency.mtx) <- names(Bic.list.select)
    for (i in seq_len(nrow(ntwork.adjacency.mtx))) {
        for (j in seq_len(ncol(ntwork.adjacency.mtx))) {
            tmp.block.1 <- Bic.list.select[[rownames(ntwork.adjacency.mtx)[i]]]
            tmp.block.2 <- Bic.list.select[[colnames(ntwork.adjacency.mtx)[j]]]
            if (edge.by == "gene") {
                n.intersect <- length(intersect(rownames(tmp.block.1), rownames(tmp.block.2)))
                n.intersect.proportion <- n.intersect/length(union(rownames(tmp.block.1), rownames(tmp.block.2)))
            } else if (edge.by == "cell") {
                n.intersect <- length(intersect(colnames(tmp.block.1), colnames(tmp.block.2)))
                n.intersect.proportion <- n.intersect/length(union(colnames(tmp.block.1), colnames(tmp.block.2)))
            } else {
                stop(paste("Please select 'gene' or 'cell' to edge.by parameter"))
            }
            ntwork.adjacency.mtx[i, j] <- n.intersect.proportion
        }
    }
    edge.list <- graph_from_adjacency_matrix(adjmatrix = ntwork.adjacency.mtx, mode = "undirected", weighted = TRUE, diag = FALSE)
    label = rownames(ntwork.adjacency.mtx)
    if ( lay.out == "linear"){
        p <- ggraph(edge.list, layout = lay.out) + 
            geom_node_point(size = 3, color = "gray50") + 
            geom_edge_arc(aes(width = weight), alpha = 0.8, color = "orange") + 
            scale_edge_width(range = c(0.5, 2)) + 
            geom_node_text(aes(label = label), repel = TRUE) + 
            labs(edge_width = paste0(edge.by," percentage\n between\n two nodes")) + 
            theme_void() 
    } else if( lay.out!= "linear"){
        p <- ggraph(edge.list, layout = lay.out) + 
            geom_node_point(size = 3, color = "gray50") + 
            geom_edge_link(aes(width = weight),alpha = 0.8, color = "orange") +
            #geom_edge_arc(aes(width = weight), alpha = 0.8, color = "orange") + 
            scale_edge_width(range = c(0.5, 2)) + 
            geom_node_text(aes(label = label), repel = TRUE) + 
            labs(edge_width = paste0(edge.by," percentage\n between\n two nodes")) + 
            theme_void()
    }
    print(p)
}


#' @rdname PlotNetwork
#' @export
#'
setMethod("PlotNetwork", "IRISFGM", .plotnetwork)



#' GenerateNetObject
#' Generate Net Object for the network use.
#' @param object Input IRIS-FGM object.
#'
#' @param N.bicluster Should be two integers indicating the number of two biclusters.
#' @param method Should be a statistical method to calculate edge weight based on expression data. It can be either "Spearman" (default) or "Pearson."
#'
#' @importFrom  stats cor
.generateNetObject <- function(object, N.bicluster = c(1, 5), method = "spearman") {
    groups <- N.bicluster
    tmp.data <- object@Processed_count
    if (length(N.bicluster) < 1) {
        stop("please choose number of biclusters")
    }
    rownamelist <- list()
    colnamelist <- list()
    for (i in seq_along(N.bicluster)) {
        idx <- i
        rownamelist[[paste0("Bicluster_", N.bicluster[idx])]] <- object@BiCluster@CoReg_gene$Gene[object@BiCluster@CoReg_gene$Condition == N.bicluster[idx]]
        colnamelist[[paste0("Bicluster_", N.bicluster[idx])]] <- object@BiCluster@CoCond_cell$cell_name[object@BiCluster@CoCond_cell$Condition == N.bicluster[idx]]
    }
    allrownames <- Reduce(union, rownamelist)
    allcolnames <- Reduce(union, colnamelist)
    un <- tmp.data[allrownames, allcolnames]
    rowidlist <- list()
    index <- paste0("Bicluster_", N.bicluster)
    # Bic.list <- .separateBic(object) bics <- Bic.list
    if (length(groups) > 2) 
        stop("length(group) > 2")
    if (length(groups) == 1) {
        rowidlist[[index[[1]]]] <- match(rownamelist[[index[[1]]]], rownames(un))
    } else if (length(groups) == 2) {
        rowidlist[[paste0(index[[1]], " & ", index[[2]])]] <- match(intersect(rownamelist[[index[[1]]]], rownamelist[[index[[2]]]]), rownames(un))
        rowidlist[[index[[1]]]] <- match(setdiff(rownamelist[[index[[1]]]], rownamelist[[index[[2]]]]), rownames(un))
        rowidlist[[index[[2]]]] <- match(setdiff(rownamelist[[index[[2]]]], rownamelist[[index[[1]]]]), rownames(un))
        rowidlist[["Others"]] <- match(setdiff(allrownames, union(rownamelist[[index[[1]]]], rownamelist[[index[[2]]]])), rownames(un))
    }
    cort <- stats::cor(t(un), method = method)
    return(list(cort, rowidlist))
}


#' @title PlotModuleNetwork
#' @description This function will visualize co-expression gene network based on selected two biclusters. 
#' The nodes represent the gene module network from the selected bicluster. 
#' The size of the nodes indicates the degree of presence.
#'  The thickness of edges indicates the value of the correlation coefficient.
#'
#' @param object Input IRIS-FGM object.
#' @param method Should be a statistical method to calculate edge weight based on expression data. It can be either "Spearman" (default) or "Pearson."
#' @param cutoff.neg Should be a cutoff to show a negative correlation between two nodes (default: -0.8).
#' @param cutoff.pos Should be a cutoff to show a positive correlation between two nodes (default: 0.8).
#' @param layout Should be one type of layouts to show nodes' arrangement, including 'linear', 'star', 'circle'(default), 'gem', 'dh', 'graphopt', 'grid', 'mds', 
#' 'randomly', 'fr', 'kk', 'drl', 'lgl'.
#' @param node.label Should be logic to show the nodes' label (default: TRUE).
#' @param node.label.cex Should be a number to control the label size.
#' @param N.bicluster Should be the two numbers of biclsuters.
#' @param node.col Should a color name (or color code) for nodes
#'
#' @importFrom igraph graph_from_adjacency_matrix degree vertex.attributes edge.attributes vertex.attributes<- edge.attributes<-
#' @import ggraph
#' @return It will generate co-expression network based on selected bicluster (can be one or multiple.)
#' @name PlotModuleNetwork
#'
#' @examples \dontrun{
#' object <- PlotModuleNetwork(object = NULL, 
#' N.bicluster = c(1,5), 
#' Node.color = '#E8E504', 
#' cutoff=0.7, 
#' node.label.cex = 1) 
#' }
.plotmodulenetwork <- function(object = NULL, 
                               method = "spearman",
                               node.col = "orange",
                               N.bicluster = c(1, 5),
                               cutoff.neg = -0.8, 
                               cutoff.pos = 0.8, 
                               layout = "circle", 
                               node.label = TRUE, 
                               node.label.cex = 1) 
{
    my.list <- .generateNetObject(object = object, N.bicluster = N.bicluster, method = method)
    cort <- my.list[[1]]
    my.adjacency <- ifelse(cort < cutoff.neg | cort > cutoff.pos, cort, 0)
    g <- graph_from_adjacency_matrix(my.adjacency, weighted = TRUE, diag = FALSE, mode = "undirected")
    if (length(N.bicluster) > 1) {
        vertex.attributes(g)$Bicluster <- c(rep(names(my.list[[2]][1]), length(my.list[[2]][[1]])), rep(names(my.list[[2]][2]), length(my.list[[2]][[2]])), 
                                            rep(names(my.list[[2]][3]), length(my.list[[2]][[3]])), rep(names(my.list[[2]][4]), length(my.list[[2]][[4]])))
    }
    if (length(N.bicluster) == 1) {
        vertex.attributes(g)$Bicluster <- rep(names(my.list[[2]][1]), length(my.list[[2]][[1]]))
    }
    edge.attributes(g)$status <- ifelse(edge.attributes(g)$weight < 0, "negative", "positive")
    degree_number <- 10*(degree(g)-min(degree(g))/(max(degree(g)-min(degree(g)))))
    layout <- create_layout(g, layout = layout)
    # create color panel
    p.base <- ggraph(layout)
    p.base <- p.base + geom_node_point(aes(size = degree_number), color = node.col)
    if (node.label == TRUE) {
        p.base <- p.base + geom_node_text(aes(label = name, size = node.label.cex), repel = TRUE)
    }
    p.base <- p.base + geom_edge_link0(aes(width = abs(weight), color = status), alpha = 1) + scale_edge_width(range = c(0.5, 2)) + coord_fixed() + theme_void()
    print(p.base)
}


#' @rdname PlotModuleNetwork
#' @export
#'
setMethod("PlotModuleNetwork", "IRISFGM", .plotmodulenetwork)



