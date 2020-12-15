#' @include generics.R
#' @include Classes.R
NULL

.RunGO <- function(genes.use = NULL, species = "mouse") {
    if (grepl("mouse", species, ignore.case = T)) {
        pathway <- invisible(suppressMessages(enrichGO(gene = genes.use, OrgDb = org.Mm.eg.db, ont = "ALL", keyType = "SYMBOL", pAdjustMethod = "BH", pvalueCutoff = 0.05, 
            qvalueCutoff = 0.05)))
    } else if (grepl("human", species, ignore.case = T)) {
        pathway <- invisible(suppressMessages(enrichGO(gene = genes.use, OrgDb = org.Hs.eg.db, ont = "ALL", keyType = "SYMBOL", pAdjustMethod = "BH", pvalueCutoff = 0.05, 
            qvalueCutoff = 0.05)))
    }
    return(pathway)
}

.IDConvert <- function(genes.use = NULL, species = NULL) {
    if (grepl("mouse", species, ignore.case = T)) {
        gene.convert <- bitr(genes.use, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db)
    } else if (grepl("human", species, ignore.case = T)) {
        gene.convert <- bitr(genes.use, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Hs.eg.db)
    }
    return(gene.convert)
}

.RunKEGG <- function(genes.use = NULL, species = "mouse") {
    if (grepl("mouse", species, ignore.case = T)) {
        pathway <- invisible(suppressMessages(enrichKEGG(gene = genes.use, organism = "mmu", keyType = "kegg", pAdjustMethod = "BH", pvalueCutoff = 0.05, 
            qvalueCutoff = 0.05)))
    } else if (grepl("human", species, ignore.case = T)) {
        pathway <- invisible(suppressMessages(enrichKEGG(gene = genes.use, organism = "hsa", keyType = "kegg", pAdjustMethod = "BH", pvalueCutoff = 0.05, 
            qvalueCutoff = 0.05)))
    }
    return(pathway)
}
#' Functional enrichment analysis
#' 
#' This function will perform enrichment analysis based on a gene module or identified differentially expressed genes (DEG).
#' This function is also depended on clusterProfiler, AnnotationDbi, org.Mm.eg.db, and org.Hs.eg.db package.
#' @param object Input IRIS-FGM object
#'
#' @param N.bicluster Select the numebr of bicluster to perform this function.
#' @param selected.gene.cutoff Set up a statistical significance cutoff for all identified DEGs.
#' @param species You can choose either 'Human' or 'Mouse'
#' @param database You can choose either 'GO' or 'KEGG' database
#' @param genes.source You can choose a gene list source, either 'CTS' or 'Bicluster.' 'CTS' means from cell-type-specific DEGs, 
#' and 'Bicluster means using gene module from the selected bicluster.'
#'
#' @importFrom clusterProfiler enrichKEGG enrichGO bitr
#' @import org.Mm.eg.db org.Hs.eg.db
#' @importFrom AnnotationDbi select
#' @name RunPathway
#' @return It will reture a function enrichment analysis.
#' @examples # If you want to perform this function based on identified DEGs, you should use: 
#' \dontrun{object <- RunPathway(object = NULL,N.bicluster = NULL, selected.gene.cutoff = 0.05,
#' species = 'Human', database = 'GO', genes.source = 'CTS' }
#' # If you want to perform this function based on the gene module from an identified bicluster, you should use: 
#' \dontrun{object <- RunPathway(object = NULL,N.bicluster = NULL, selected.gene.cutoff = 0.05,
#' species = 'Human', database = 'GO', genes.source = 'Bicluster' }
.runPathway <- function(object = NULL, N.bicluster = c(10, 11, 12, 13), selected.gene.cutoff = 0.05, species = "Human", database = "GO", genes.source = c("CTS", 
    "Bicluster")) {
    
    if (genes.source == "CTS") {
        tmp.table <- as.data.frame(object@LTMG@MarkerGene)
        genes.use.LTMG <- rownames(tmp.table)[tmp.table$pvalue.adj.FDR < selected.gene.cutoff]
        if (database == "GO") {
            pathway <- .RunGO(genes.use = genes.use.LTMG, species = species)
        } else if (database == "KEGG") {
            ID.covert <- .IDConvert(genes.use.LTMG, species = species)
            pathway <- .RunKEGG(genes.use = ID.covert$ENTREZID, species = species)
        }
        object@LTMG@Pathway <- pathway
    } else if (genes.source == "Bicluster") {
        if (length(N.bicluster) < 1) {
            stop("please choose number of biclusters")
        }
        rownamelist <- list()
        for (i in seq_along(N.bicluster)) {
            idx <- i
            rownamelist[[paste0("Bicluster_", N.bicluster[idx])]] <- object@BiCluster@CoReg_gene$Gene[object@BiCluster@CoReg_gene$Condition == N.bicluster[idx]]
        }
        allrownames <- Reduce(union, rownamelist)
        genes.use.module <- allrownames
        # run on Bicluster marker gene
        if (is.null(object@BiCluster@MarkerGene)) {
            message("There is no gene in MarkerGene slot.
                                                     \n Ignore pathway analysis based on marker gene derived from MC defined cell type. ")
            genes.use.MC <- NULL
        } else {
            tmp.table <- object@BiCluster@MarkerGene
            genes.use.MC <- rownames(tmp.table)[tmp.table$pvalue.adj.FDR < selected.gene.cutoff]
        }
        gene.list <- c(genes.use.module = list(genes.use.module), genes.use.MC = list(genes.use.MC))
        if (database == "GO") {
            pathway <- lapply(gene.list, function(x) .RunGO(genes.use = x, species = species))
        } else if (database == "KEGG") {
            ID.covert <- lapply(gene.list, function(x) .IDConvert(genes.use = x, species = species))
            pathway <- lapply(ID.covert, function(x) .RunKEGG(genes.use = x$ENTREZID, species = species))
        }
        object@BiCluster@PathwayFromModule <- pathway$genes.use.module
        if (is.null(genes.use.MC)) {
            message("not find MC based marker genes, ignore pathway enrichment.......")
        } else {
            object@BiCluster@PathwayFromMC <- pathway$genes.use.MC
        }
    }
    
    return(object)
    
}

#' @rdname RunPathway
#' @export
setMethod("RunPathway", "IRISFGM", .runPathway)

#' Plot dotplot for enrichment pathway
#'
#' @param object Input IRIS-FGM object
#' @param genes.source decide the plot source either 'CTS' or 'Bicluster.' 'CTS' means from cell-type-specific DEGs, 
#' and 'Bicluster means using gene module from the selected bicluster.'
#' @importfrom clusterProfiler dotplot
#'
#' @return dot plot
#' @name DotPlotPathway
#' 
#' @examples
.dotplotpathway <- function(object = object, genes.source = c("CTS", "MC", "module")) {
    if (genes.source == "CTS") {
        print(dotplot(object@LTMG@Pathway, showCategory = 20))
    }
    if (genes.source == "MC") {
        print(dotplot(object@BiCluster@PathwayFromMC, showCategory = 20))
    }
    if (genes.source == "Bicluster") {
        print(dotplot(object@BiCluster@PathwayFromModule, showCategory = 20))
    }
}

#' @rdname DotPlotPathway
#' @export 
setMethod("DotPlotPathway", "IRISFGM", .dotplotpathway)
