

#' Wrapper for using scVI for list of Seurat objects to correct for batch effects
#'
#' @param object.list A list of Seurat objects
#' @param batch.key Metadata variable or variables to indicate batch
#' @param features Either a list of variable features or numeric to calculate number of variable features used for integration
#' @param assay Assay to calculate scVI. Default is RNA from Seurat
#' @param merge.data If list of objects are normalized using the same method, values are carried over when merging objects
#' @param vars.to.regress List of features to regress such as cell cycle, if calculated on each Seurat object
#'
#' @return Integrated Seurat object using scVI
#' @export
#' @import Seurat
#' @importFrom reticulate import
#' @importFrom sceasy convertFormat
#'
#' @examples
#'
RunSCVI <- function(object.list, batch.key=NULL, features=2000, assay='RNA', merge.data=TRUE, vars.to.regress=NULL) {

    if (is.character(features)) {
        features <- features
    } else if(is.numeric(features)) {
        features <- SelectIntegrationFeatures(object.list=object.list, nfeatures=features)
    } else {
        stop('features not provided...try again')
    }

    object <- merge(x=object.list[[1]], y=object.list[-1], merge.data=merge.data)

    sc <- reticulate::import('scanpy', convert=FALSE)
    scvi <- reticulate::import('scvi', convert=FALSE)

    adata <- sceasy::convertFormat(object[features], from='seurat', to='anndata', main_layer='counts', drop_single_values=FALSE)

    # run setup_anndata, use column for batch if applicable
    if (is.character(vars.to.regress)) {
        scvi$model$SCVI$setup_anndata(adata, batch_key=batch.key, continuous_covariate_keys=as.list(vars.to.regress))
    } else {
        scvi$model$SCVI$setup_anndata(adata, batch_key=batch.key)
    }

    # create the model
    model <- scvi$model$SCVI(adata)

    # train the model
    model$train()

    # get the latent represenation
    latent <- model$get_latent_representation()

    # put it back in our original Seurat object
    latent <- as.matrix(latent)
    rownames(latent) <- colnames(object)
    object[['scvi']] <- CreateDimReducObject(embeddings=latent, key='scvi_', assay=assay)

    VariableFeatures(object) <- features

    return(object)
}
