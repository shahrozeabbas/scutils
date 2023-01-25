
#' Wrapper function to run Kmeans Clustering on a Seurat object
#'
#' @param object An object of class Seurat
#' @param k Parameter for number of clusters in kmeans
#' @param reduction Dimension reduction to use from Seurat object for clustering
#' @param ... Any other parameters for Base R kmeans function
#'
#' @return A Seurat object with kmean clusters added as metadata
#' @export
#' @import Seurat
#' @importFrom stats kmeans
#'
#' @examples
#'
FindClustersKmeans <- function(object, k=2, reduction='pca', ...) {

  pca <- Embeddings(object=object, reduction=reduction)
  groups <- stats::kmeans(x=pca, centers=k, ...)$cluster

  object <- AddMetaData(object=object, metadata=factor(groups), col.name=paste0('kmeans_', k))

  return(object)
}




