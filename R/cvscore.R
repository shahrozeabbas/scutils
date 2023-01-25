#' CV Score Wrapper for Seurat object as described by Lakkis et. al (2021)
#'
#' @param object An object of class Seurat
#' @param groupby Clustering metadata variable
#' @param batch Batch metadata variable
#' @param weighted Boolean to determine if weighted score should be calculated for all clusters
#'
#' @return A vector of scores
#' @export
#' @import Seurat data.table
#' @importFrom methods is
#' @importFrom stats sd
#'
#' @examples
#'
CVScore <- function(object, groupby='seurat_clusters', batch='tumorID', weighted=TRUE) {

    if (!is(object, 'Seurat')) stop('please provide a Seurat object...')

    lambda <- 10^-12
    m <- copy(object@meta.data)
    setDT(m, keep.rownames='cells')

    group_table <- m[, .N, keyby=groupby][, prop := N / sum(N)]


    groups <- as.character(group_table[, get(groupby)])

    scores <- sapply(X=groups, FUN=function(group) {
        subcells <- m[get(groupby) %in% group, cells]
        avg <- as.data.frame(object %>% subset(cells=subcells) %>%
                                AverageExpression(assays='RNA', group.by=batch, slot='data'))

        sdVec <- apply(X=avg, FUN=sd, MARGIN=1)
        meanVec <- apply(X=avg, FUN=mean, MARGIN=1)

        if (weighted) {
            (sdVec / (meanVec + lambda)) * group_table[get(groupby) %in% group, prop]
        } else {
            sdVec / (meanVec + lambda)
        }
    })

    return(scores)
}
