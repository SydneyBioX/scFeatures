#' A toy scRNA-seq data 
#'
#' This data is included in the package for demonstration purposes. 
#' It is a subsampled version of a melanoma data, containing treatment
#' response status (ie, responder and non-responder) of melanoma patients. 
#'
#' @docType data
#' 
#' @usage data(example_scrnaseq)
#' 
#' @format A \code{Seurat} object with 3,523 rows and 550 columns.
#' Each row is a gene and each column is a cell. 
#' It include metadata that describe each cell. 
#' Some important metdata columns are:
#' \describe{
#'   \item{celltype}{Cell type label of the cell}
#'   \item{sample}{Which sample the cell is from}
#'   \item{patient}{The patient ID}
#'   \item{respond}{Response under therapy}
#'   \item{pre_post}{Pre or post therapy}
#' }
#' 
#' @references Sade-Feldman, M., Yizhak, K., Bjorgaard, S. L., Ray, J. P., 
#' de Boer, C. G., Jenkins, R. W., ... & Hacohen, N.(2018). Defining T cell 
#' states associated with response to checkpoint immunotherapy in melanoma. 
#' Cell, 175(4), 998-1013.
#' 
#' @source <https://www.sciencedirect.com/science/article/pii/S0092867418313941>
#' 
#' @examples
#' data <- data(example_scrnaseq)
#' str(data)
#' 
"example_scrnaseq"