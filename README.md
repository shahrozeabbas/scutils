# scutils - Wrapper Functions for Analyzing Single-Cell Data with Seurat

### Package Overview
Over the years, I've found useful packages for analyzing single-cell data. As a result, there are many different approaches to handle problems one might have when analyzing these data. The `Seurat` package has become a popular choice for researchers to perform end to end analysis. However, many tools may not be ready to accept `Seurat` objects out of the box. I designed this package to keep track of the different tools I have used and found useful, and to allow users to easily apply these tools to data stored in `Seurat` objects.      

### Dependencies
R packages from CRAN that are required to run `scutils`

```
install.packages(c('dplyr', 'Seurat', 'data.table', 'reticulate', 'SoupX', 'DUBStepR'))
```

The `sceasy` package is required to run `scVI`

```
devtools::install_github('cellgeni/sceasy')
```


### Installation
This package can be installed using the following code in R. Please note that some functions will require external dependencies to be downloaded and installed. 

```
devtools::install_github('shahrozeabbas/scutils')
```

<br />
<br />

# Functions Currently Supported by `scutils`

## Removing Mitochondrial and Ribosomal Genes from Variable Feature List
Many machine learning algorithms rely on a set of variable or informative features to create models. This becomes especially necessary for single-cell data which are sparse. Many genes do not provide enough information to be used for downstream clustering. `Seurat` calculates a set a variable features for this purpose, which can sometimes retain mitochondrial or ribosomal genes. Depending on what biological question is, this may inform downstream clustering in a negative way. The `CleanVarGenes` function can be used to remove these genes from the variable feature list. 


## Calculate Variable Features Using DUBstepR
*Ranjan et al.* created a novel way to calculate variable features for single-cell data. The authors use pairwise gene-gene correlations to produce a list of genes that are robust across different conditions for the same cell types. It can be a useful tool in your arsenal and can be run using the `dupstepR` function. For more information about the method please refer to their [paper](https://www.nature.com/articles/s41467-021-26085-2). To install DUBstepR, please use [CRAN](https://cran.r-project.org/web/packages/DUBStepR/index.html).

## Use K-means Clustering on Reduced Dimensions in Seurat
K-means clustering is a classic machine learning approach to group data into similar categories or clusters. The `FindClustersKmeans` function in this package is applied typically to PCA dimensions. There are many resources available online to learn more about k-means, but you can start [here](https://en.wikipedia.org/wiki/K-means_clustering). 

## Correct for Ambient RNA contamination using SoupX
Droplet-based Single-cell experiments often capture mRNA that is present in the assay solution and encapsulated along with a given cell. Before sequencing, these cells are lysed and 'cell-free' RNA is added to the mixture, contaminating the downstream expression counts. The `SoupX` method attempts to correct this problem and the `SoupCorrect` wrapper function included in this package helps to facilitate this process. Please read more about this awesome package [here](https://academic.oup.com/gigascience/article/9/12/giaa151/6049831).

## Using the Coefficient of Variation Score to Assess Integration of Multiple Datasets
Batch effect is a major issue when combining multiple datasets in single-cell analyses. There are currently many packages for addressing this problem along with ways to assess how well the integration performed. *Lakkis et al.* created the [CarDEC](https://genome.cshlp.org/content/early/2021/05/25/gr.271874.120) method for this purpose. The authors also define a metric known as the *CV Score* to assess how well their method and others are able to correct for batch effect. This package implements the metric in an easy-to-use `CVScore` function that can be included in any single-cell analysis. 

## Integrate Data using scVI and Seurat
This is another method that is part of a toolkit and stands for [Single-Cell Variational Inference](https://scvi-tools.org) and is a personal favorite. `scVI` uses probabilistic models to integrate and remove batch effect. Similar to `Seurat`, the toolkit provides methods for many different types of data such as scATAC-seq and spatial sequencing data. The function in this package is written to accept a `Seurat` object and allows the user to supply batch labels present in the dataset. It's extremely fast and can be used to build atlases with reasonable runtimes. The website linked above contains tutorials for all their methods. The `scVI` package must be installed seperately and can be run using the `RunSCVI` function.

