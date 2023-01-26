#' Wrapper for sctype by Lanevski et al.
#'
#' @param object An object of class Seurat
#' @param genes.list A named list
#' @param assay The assay from Seurat to use. Default is RNA
#' @param is.scaled Is the input data scaled. Should be scaled
#' @param tissue Tissue type supplied if gene list of markers is not
#'
#' @return
#' @export
#' @import data.table Seurat dplyr
#' @importFrom methods is
#' @importFrom utils head
#' @importFrom openxlsx read.xlsx
#'
#' @examples
#'
#'
AnnotateCellTypes <- function(object, genes.list=NULL, assay='RNA', is.scaled=TRUE, tissue=NULL) {

    source('https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R')
    source('https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R')
    db_ <- 'https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx'

    if (!is(object, 'Seurat')) stop('please provide a Seurat object...')

    tissue <- ifelse(!is.null(gene.list), NULL, tissue)

    gene.list <- ifelse(is.null(gene.list), gene_sets_prepare(db_, tissue)$gs_positive, gene.list)

    # get cell-type by cell matrix
    counts <- object %>% GetAssayData(assay=assay, slot=ifelse(is.scaled, 'scale.data', 'data'))
    es.max <- sctype_score(scRNAseqData=counts, scaled=is.scaled, gs=genes.list)

    # merge by cluster
    cL_resutls <- do.call('rbind', lapply(unique(object@meta.data$seurat_clusters), function(cl) {

        es.max.cl <- sort(rowSums(es.max[, rownames(object@meta.data[object@meta.data$seurat_clusters == cl, ])]), decreasing=!0)
        head(data.frame(cluster=cl, type=names(es.max.cl), scores=es.max.cl, ncells=sum(object@meta.data$seurat_clusters == cl)), 10)

        }
    ))

    sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n=1, wt=scores)

    # set low-confident (low ScType score) clusters to 'unknown'
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells / 4] <- 'Unknown'

    setDT(sctype_scores)
    setnames(sctype_scores, old='cluster', new='seurat_clusters')

    return(sctype_scores)
}

