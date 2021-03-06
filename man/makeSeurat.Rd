% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapper_run_scfeatures.R
\name{makeSeurat}
\alias{makeSeurat}
\title{Format SingleCellExperiment or SpatialExperiment object into Seurat object for scFeatures functions}
\usage{
makeSeurat(data, sample, celltype, assay, spatialCoords, spotProbability)
}
\arguments{
\item{data}{input data, either a SingleCellExperiment or SpatialExperiment object.
The object needs to contain a column named "sample" and a column named "celltype".
Alternatively, user can provide the name of the column containing sample and celltype into the \code{sample} and \code{celltype} argument.
When passing as SingleCellExperiment or SpatialExperiment, by default we use the assay stored in "logcount".
Alternatively, user can specify the assay to use in the \code{assay} argument.
If users want to construct features from the spatial category, by default we need columns called "x_cord" and "y_cord".
Alternatively, please specify the relevant column in the \code{spatialCoords} argument.
For spot-based spatial transcriptomics, we also requires a matrix containing cell type prediction probability of each spot, in the format of celltype x spot}

\item{sample}{the name of column containing sample identifier. If not provided, we use the column called "sample".}

\item{celltype}{the name of column containing sample identifier. If not provided, we use the column called "sample".}

\item{assay}{the assay identifier if using a SingleCellExperiment or SpatialExperiment object.}

\item{spatialCoords}{the spatialCoords identifiers if users want to construct features from the spatial category.
If not provided, we use the columns called "x_cord" and "y_cord".}

\item{spotProbability}{a matrix in the format of celltype x spot, where each entry is the prediction probability of that cell type for each spot.}
}
\description{
Format SingleCellExperiment or SpatialExperiment object into Seurat object for scFeatures functions
}
