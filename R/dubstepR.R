#' Wrapper for DUBstepR to be used with Seurat objects
#'
#' @param object An object of class Seurat
#' @param k Number of nearest neighbors for DUBstepR
#' @param pcs Number of dimensions to represent single-cell data
#' @param optimize Determine optimal feature set using density index
#'
#' @return A Seurat object with variable features calculated using DUBstepR
#' @export
#' @import DUBStepR Seurat
#' @importFrom methods is
#'
#' @examples
#'
dubstepR <- function(object, k=NULL, pcs=NULL, optimize=TRUE) {

    if (!is(object, 'Seurat')) stop('please provide a Seurat object...')

    nCells <- ncol(object)
    k <- ifelse(is.null(k), round(sqrt(nCells) * 0.5), k)
    pcs <- ifelse(is.null(pcs), ifelse(nCells < 50, nCells, 50), pcs)

    message(paste0('k: ', k))
    message(paste0('pcs: ', pcs))

    e <- ifelse(nCells > 10000, 1, 0)
    input <- GetAssayData(object, slot='data', assay='RNA')
    dubstepR.out <- DUBStepR(input.data=input, min.cells=0.05*nCells, optimise.features=optimize, k=ifelse(k < 5, 5, k), num.pcs=pcs, error=e)

    VariableFeatures(object) <- dubstepR.out$optimal.feature.genes

    return(object)
}
