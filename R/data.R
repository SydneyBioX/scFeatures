#' A toy scRNA-seq data 
#'
#' This data is included in the package for demonstration purposes. 
#' 
#' It is a subsampled version of the melanoma pre-treatment dataset, from 
#' the publication Sade-Feldman, M., Yizhak, K., Bjorgaard, S. L., Ray, J. P., 
#' de Boer, C. G., Jenkins, R. W., ... & Hacohen, N.(2018). Defining T cell 
#' states associated with response to checkpoint immunotherapy in melanoma. 
#' Cell, 175(4), 998-1013.
#'
#' @format ## `example_scrnaseq`
#' A Seurat object with 7,240 rows and 60 columns. 
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
#' @source <https://www.sciencedirect.com/science/article/pii/S0092867418313941>
"example_scrnaseq"