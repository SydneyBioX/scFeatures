#' Example of scRNA-seq data
#'
#' This is a subsampled version of the melanoma patients dataset as used in our
#' manuscript.
#' The original dataset is available at:
#'      <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575>.
#' 
#' @usage data("example_scrnaseq")
#' @format ## `example_scrnaseq`
#' A Seurat object with 3523 genes and 550 cells. 
#' Some of the key metadata columns are: 
#' \describe{
#'   \item{celltype}{cell type of the cell}
#'   \item{sample}{patient ID of the cell}
#'   \item{condition}{whether the patient is a responder or non-responder}
#' }
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575>
"example_scrnaseq"


#' Example of scFeatures() output
#'
#' This is an example output of the scFeatures() function for example_scrnaseq. 
#'
#' @usage data("scfeatures_result")
#' @format ## `scfeatures_result`
#' A list with two dataframes. In each dataframe the columns are each patient and the rows are the feature values.
#' The first dataframe contains the feature type "proportion_raw". 
#' The second dataframe contains the feature type "proportion_logit". 
"scfeatures_result"