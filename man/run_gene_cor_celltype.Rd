% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_scfeatures.R
\name{run_gene_cor_celltype}
\alias{run_gene_cor_celltype}
\title{Generate cell type specific gene expression correlation}
\usage{
run_gene_cor_celltype(
  data,
  type = "scrna",
  genes = NULL,
  num_top_gene = NULL,
  ncores = 1
)
}
\arguments{
\item{data}{A list object containing \code{data} matrix and \code{celltype} and \code{sample} vector.}

\item{type}{The type of dataset, either "scrna", "spatial_t", or "spatial_p".}

\item{genes}{Optional dataframe with 2 columns: 'marker' and 'celltype'.
The 'marker' column should contain the genes of interest (e.g. 'S100A11', 'CCL4'),
and the 'celltype' column should contain the celltype that the gene expression
is to be computed from (e.g. 'CD8', 'B cells').
If not provided, the top variable genes will be used based on the
num_top_gene parameter.}

\item{num_top_gene}{Number of top genes to use when genes is not provided.
Defaults to 5.}

\item{ncores}{Number of cores for parallel processing.}
}
\value{
a dataframe of samples x features.
The features are in the form of gene 1 vs gene 2 cell type a ,
gene 1 vs gene 3 cell type b ... etc, with the numbers representing the
correlation of the two given genes in the given cell type.
}
\description{
This function computes the correlation of expression of a set of genes
for each cell type in the input data. The input data can be of three types:
'scrna', 'spatial_p' or 'spatial_t'. If the genes parameter is not provided
by the user, the top variable genes will be selected based on the
num_top_gene parameter (defaults to 100).
}
\examples{

utils::data("example_scrnaseq" , package = "scFeatures")
data <- example_scrnaseq[1:50, 1:20]
celltype <- data$celltype
sample <- data$sample
data <- data@assays$RNA@data

alldata <- scFeatures:::formatData(data = data, celltype = celltype, sample = sample )

feature_gene_cor_celltype <- run_gene_cor_celltype(
   alldata,
   type = "scrna", num_top_gene = 5, ncores = 1
 )

}
