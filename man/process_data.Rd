% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{process_data}
\alias{process_data}
\title{perform pre-processing}
\usage{
process_data(data, normalise = T)
}
\arguments{
\item{data}{input data}

\item{normalise}{whether to normalise the data. Note if the data has already been normalised (eg, log2CPM), there is no need to normalise again.}
}
\value{
a matrix of samples x features
}
\description{
perform pre-processing
}
