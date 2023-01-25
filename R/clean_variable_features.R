
#' Remove Mitochondrial and Ribosomal Genes from Variable Features in Seurat
#'
#' @param object An object of class Seurat
#' @param nHVG Number of variable features to clean
#'
#' @return A Seurat object with variable features that do not contain mito or rb genes
#' @export
#' @import Seurat dplyr data.table
#' @importFrom methods is
#'
#' @examples
#'
CleanVarGenes <- function(object, nHVG=2000) {

    if (!is(object, 'Seurat')) stop('please provide a Seurat object...')

    g <- data.table(genes=rownames(object))

    genes <- object %>%
        subset(features=g[!(genes %like% '^MT' | genes %like% '^RP[LS]'), genes]) %>%
        NormalizeData() %>% FindVariableFeatures(nfeatures=nHVG) %>% VariableFeatures()

    VariableFeatures(object) <- genes

    return(object)
}
