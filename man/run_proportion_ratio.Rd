% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_scfeatures.R
\name{run_proportion_ratio}
\alias{run_proportion_ratio}
\title{generate cell type proportion ratio}
\usage{
run_proportion_ratio(data, type = "scrna", ncores = 1)
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
generate cell type proportion ratio
}
\examples{

data <- readRDS(system.file("extdata", "example_scrnaseq.rds", package = "scFeatures"))
feature_proportion_ratio <- run_proportion_ratio(data, type = "scnrna", ncores = 1)

}
