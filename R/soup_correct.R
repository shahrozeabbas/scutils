#' Wrapper for Ambient RNA removal using SoupX
#'
#' @param raw.matrix Raw expression matrix from 10x Chromium
#' @param filt.matrix Filtered expression matrix 10x Chromium
#' @param contamination_rate Contamination rate to correct for soup
#'
#' @return Adjusted count matrix to be used downstream
#' @export
#' @import Seurat SoupX
#' @importFrom stats setNames
#'
#' @examples
#'
SoupCorrect <- function(raw.matrix, filt.matrix, contamination_rate=NULL) {

  srat  <- CreateSeuratObject(counts=filt.matrix[['Gene Expression']])
  soup.channel  <- SoupChannel(raw.matrix[['Gene Expression']], filt.matrix[['Gene Expression']])

  srat <- srat %>%
    SCTransform(verbose=FALSE) %>%
    RunPCA(verbose=FALSE) %>%
    RunUMAP(dims=1:30, verbose=FALSE) %>%
    FindNeighbors(dims=1:30, verbose=FALSE) %>%
    FindClusters(verbose=FALSE, algorithm=3)

  meta <- srat@meta.data
  umap <- srat@reductions$umap@cell.embeddings
  soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel <- setDR(soup.channel, umap)

  if (is.null(contamination_rate)) {
    soup.channel  <- autoEstCont(soup.channel, forceAccept=TRUE, doPlot=FALSE)
  } else {
    soup.channel <- setContaminationFraction(soup.channel, contamination_rate, forceAccept=TRUE)
  }

  adj.matrix  <- adjustCounts(soup.channel)

  return(adj.matrix)

}

