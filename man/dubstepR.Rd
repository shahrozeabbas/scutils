% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dubstepR.R
\name{dubstepR}
\alias{dubstepR}
\title{Wrapper for DUBstepR to be used with Seurat objects}
\usage{
dubstepR(object, k = NULL, pcs = NULL, optimize = TRUE)
}
\arguments{
\item{object}{An object of class Seurat}

\item{k}{Number of nearest neighbors for DUBstepR}

\item{pcs}{Number of dimensions to represent single-cell data}

\item{optimize}{Determine optimal feature set using density index}
}
\value{
A Seurat object with variable features calculated using DUBstepR
}
\description{
Wrapper for DUBstepR to be used with Seurat objects
}
