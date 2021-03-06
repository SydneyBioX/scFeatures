% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_scfeatures.R
\name{run_L_function}
\alias{run_L_function}
\title{generate L stats}
\usage{
run_L_function(data, type = "spatial_p", ncores = 1)
}
\arguments{
\item{data}{input data, a Seurat object containing \code{celltype} and \code{sample} label}

\item{type}{input data type, either scrna, spatial_p, or spatial_t}

\item{ncores}{number of cores}
}
\value{
a matrix of samples x features
}
\description{
generate L stats
}
\examples{

data <- readRDS(system.file("extdata", "example_spatial_proteomics.rds", package = "scFeatures"))
feature_L_function <- run_L_function(data, type = "spatial_p", ncores = 1 )

}
