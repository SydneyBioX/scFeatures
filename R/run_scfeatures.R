# this file contains the 17 functions that generate the 17 feature classes



#' generate cell type proportion raw
#'
#' @param data the input data, a Seurat object containing `celltype` and
#'             `sample` label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores
#'
#' @return a matrix of samples x features
#'
#' @examples
#'
#' data <- readRDS(
#'     system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#' )
#' feature_proportion_raw <- run_proportion_raw(
#'     data,
#'     type = "scrna", ncores = 1
#' )
#'
#' @importFrom gtools logit
#' @importFrom tidyr pivot_wider
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom cli cli_abort cli_warn
#'
#' @export
run_proportion_raw <- function(data, type = "scrna", ncores = 1) {
    check_data(data, type)

    if (type %in% c("scrna", "spatial_p")) {
        if (ncores > 1) {
            cli::cli_warn(
                "{type} does not currently support parallel computation ",
                "for feature `proportion_raw`"
            )
        }
        X <- helper_proportion_raw(data, logit = FALSE)
    } else if (type == "spatial_t") {
        X <- helper_proportion_raw_st(data, logit = FALSE, ncores)
    } else {
        cli::cli_abort(c(
            "Parameter {.var type} must be 'scrna', 'spatial_p' or 'spatial_t'",
            "x" = "'{type}' is not a valid input data type"
        ))
    }

    X <- as.data.frame(X)
    return(X)
}




#' generate cell type proportion raw and cell type proportion logit transformed
#'
#' @param data the input data, a Seurat object containing `celltype` and
#'             `sample` label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores
#'
#' @return a matrix of samples x features
#'
#' @examples
#'
#' data <- readRDS(
#'     system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#' )
#' feature_proportion_logit <- run_proportion_logit(
#'     data,
#'     type = "scrna", ncores = 1
#' )
#'
#' @importFrom gtools logit
#' @importFrom tidyr pivot_wider
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom cli cli_abort cli_warn
#'
#' @export
run_proportion_logit <- function(data, type = "scrna", ncores = 1) {
    check_data(data, type)

    if (type %in% c("scrna", "spatial_p")) {
        if (ncores > 1) {
            cli::cli_warn(
                "{type} does not currently support parallel computation ",
                "for feature proportion_logit."
            )
        }
        X <- helper_proportion_raw(data, logit = TRUE)
    } else if (type == "spatial_t") {
        X <- helper_proportion_raw_st(data, logit = TRUE, ncores)
    } else {
        cli::cli_abort(c(
            "Parameter {.var type} must be 'scrna', 'spatial_p' or 'spatial_t'",
            "x" = "'{type}' is not a valid input data type"
        ))
    }


    X <- as.data.frame(X)
    return(X)
}



#' generate cell type proportion ratio
#'
#' @param data input data, a Seurat object containing `celltype` and
#'             `sample` label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores
#'
#' @return a matrix of samples x features
#'
#' @examples
#'
#' data <- readRDS(
#'     system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#' )
#' feature_proportion_ratio <- run_proportion_ratio(
#'     data,
#'     type = "scrna", ncores = 1
#' )
#'
#' @importFrom tidyr pivot_wider
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom cli cli_abort
#'
#' @export
run_proportion_ratio <- function(data, type = "scrna", ncores = 1) {
    check_data(data, type)

    if (type %in% c("scrna", "spatial_p")) {
        X <- helper_proportion_ratio(data, ncores)
    } else if (type == "spatial_t") {
        X <- helper_proportion_ratio_st(data, ncores)
    } else {
        cli::cli_abort(c(
            "Parameter {.var type} must be 'scrna', 'spatial_p' or 'spatial_t'",
            "x" = "'{type}' is not a valid input data type"
        ))
    }

    X <- as.data.frame(X)
    return(X)
}










#' generate cell type specific gene mean expression
#'
#' @param data input data, a Seurat object containing `celltype` and
#'             `sample` label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param genes defaults to NULL, in which case the top variable genes are used
#'  If provided by user, need to be in the format of a dataframe with 2 columns
#'  'marker' and 'celltype'.
#'  The marker column contains the genes of interest, eg: S100A11 , CCL4 ,
#'  the celltype column contains the celltype that the gene expression is to
#'  be computed from, eg: CD8, B cells
#' @param num_top_gene when the genes is not provided by the user, the top
#'  variable genes will be used
#'  The number of genes is set by this number.
#'  Defaults to NULL, in which case top 100 genes from each cell type will be
#'  selected
#' @param ncores number of cores
#'
#' @return a matrix of samples x features
#'
#' @examples
#'
#' data <- readRDS(
#'     system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#' )
#' # optional step, if mito and ribo genes are not of interest
#' data_remove_mito <- remove_mito(data)
#' feature_gene_mean_celltype <- run_gene_mean_celltype(
#'     data_remove_mito,
#'     type = "scrna", num_top_gene = 100, ncores = 1
#' )
#'
#' @importFrom proxyC simil
#' @importFrom DelayedMatrixStats rowVars rowMeans2
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom DelayedArray DelayedArray
#' @importFrom cli cli_abort cli_warn
#'
#' @export
run_gene_mean_celltype <- function(data,
    type = "scrna",
    genes = NULL,
    num_top_gene = NULL,
    ncores = 1) {
    # check_data(data, type)

    if (type %in% c("scrna", "spatial_p")) {
        X <- helper_gene_mean_celltype(data, genes, num_top_gene, ncores)
    } else if (type == "spatial_t") {
        if (!is.null(genes) || !is.null(num_top_gene) || ncores > 1) {
            cli::cli_warn(c(
                "{.var gene}, {.var num_top_gene} and {.var ncores} are not ",
                "implemented for {type} data.\n",
                "i" = "Defaults will be used instead."
            ))
        }
        X <- helper_gene_mean_celltype_st(data)
    } else {
        cli::cli_abort(c(
            "Parameter {.var type} must be 'scrna', 'spatial_p' or 'spatial_t'",
            "x" = "'{type}' is not a valid input data type"
        ))
    }

    X <- as.data.frame(X)
    return(X)
}


# TODO: Several warnings below should probably be errors.


#' generate cell type specific gene proportion expression
#'
#' @param data input data, a Seurat object containing `celltype` and `sample`
#'  label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param genes default to NULL, in which case the top variable genes will be
#'  used
#'  If provided by user, need to be in the format of a dataframe with 2 columns
#'  'marker' and 'celltype'.
#'  The marker column contains the genes of interest, eg: S100A11 , CCL4 ,
#'  the celltype column contains the celltype that the gene expression is to be
#'  computed from, eg: CD8, B cells
#' @param num_top_gene when the genes is not provided by the user, the top
#'  variable genes will be used
#'  The number of genes is set by this number.
#'  default to NULL, in which case top 100 genes from each cell type will be
#'  selected
#' @param ncores number of cores
#'
#' @return a matrix of samples x features
#'
#' @examples
#'
#' data <- readRDS(
#'     system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#' )
#' # optional step, if mito and ribo genes are not of interest
#' data_remove_mito <- remove_mito(data)
#' feature_gene_prop_celltype <- run_gene_prop_celltype(
#'     data_remove_mito,
#'     type = "scrna", num_top_gene = 100, ncores = 1
#' )
#'
#' @importFrom proxyC simil
#' @importFrom DelayedMatrixStats rowVars rowMeans2
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom DelayedArray DelayedArray
#' @importFrom cli cli_abort cli_warn
#'
#' @export
run_gene_prop_celltype <- function(data,
    type = "scrna",
    genes = NULL,
    num_top_gene = NULL,
    ncores = 1) {
    check_data(data, type)

    if (type %in% c("scrna", "spatial_p")) {
        X <- helper_gene_prop_celltype(data, genes, num_top_gene, ncores)
    } else if (type == "spatial_t") {
        cli::cli_warn(c(
            "`gene_prop_celltype` currently does not support spatial transcriptomics"
        ))
        return(NULL)
    } else {
        cli::cli_abort(c(
            "Parameter {.var type} must be 'scrna', 'spatial_p' or 'spatial_t'",
            "x" = "'{type}' is not a valid input data type"
        ))
    }

    X <- as.data.frame(X)
    return(X)
}





#' generate cell type specific gene correlation
#'
#' @param data input data, a Seurat object containing `celltype` and `sample`
#'  label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param genes default to NULL, in which case the top variable genes will be
#'  used
#'  If provided by user, need to be in the format of a dataframe with 2 columns
#'  'marker' and 'celltype'.
#'  The marker column contains the genes of interest (e.g., S100A11 , CCL4),
#'  the celltype column contains the celltype that the gene expression is to be
#'  computed from, eg: CD8, B cells
#' @param num_top_gene when the genes is not provided by the user, the top
#'  variable genes will be used
#'  The number of genes is set by this number.
#'  default to NULL, in which case top 5 genes from each cell type will be
#'  selected
#' @param ncores number of cores
#'
#' @return a matrix of samples x features
#'
#' @examples
#'
#' data <- readRDS(system.file("extdata",
#'     "example_scrnaseq.rds",
#'     package = "scFeatures"
#' ))
#' # optional step, if mito and ribo genes are not of interest
#' data_remove_mito <- remove_mito(data)
#' feature_gene_cor_celltype <- run_gene_cor_celltype(
#'     data_remove_mito,
#'     type = "scrna", num_top_gene = 100, ncores = 1
#' )
#'
#' @importFrom proxyC simil
#' @importFrom DelayedMatrixStats rowVars rowMeans2
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom DelayedArray DelayedArray
#' @importFrom cli cli_abort cli_warn
#'
#' @export
run_gene_cor_celltype <- function(data,
    type = "scrna",
    genes = NULL,
    num_top_gene = NULL,
    ncores = 1) {
    check_data(data, type)

    if (type %in% c("scrna", "spatial_p")) {
        X <- helper_gene_cor_celltype(data, genes, num_top_gene, ncores)
    } else if (type == "spatial_t") {
        cli::cli_warn(c(
            "`gene_cor_celltype` currently does not support spatial transcriptomics"
        ))
        return(NULL)
    } else {
        cli::cli_abort(c(
            "Parameter {.var type} must be 'scrna', 'spatial_p' or 'spatial_t'",
            "x" = "'{type}' is not a valid input data type"
        ))
    }

    X <- as.data.frame(X)
    return(X)
}






#' generate pathway score using gene set enrichement analysis
#'
#' @param data input data, a Seurat object containing `celltype` and `sample`
#'  label
#' @param method type of pathway analysis method, currently support `ssgsea`
#'  and `aucell`
#' @param geneset By default (when the `geneset` argument is not specified),
#'  we use the 50 hallmark gene set from msigdb.
#'  The users can also provide their geneset of interest in a list format, with
#'  each list entry containing a vector of the names of genes in a gene set.
#'  eg, geneset <- list("pathway_a" = c("CAPN1", ...), "pathway_b" = c("PEX6"))
#' @param species whether the species is "Homo sapiens" or "Mus musculus".
#'  Default is "Homo sapiens".
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param subsample whether to subsample, either T or F. For larger datasets
#'  (eg, over 30,000 cells), the subsample function can be used to increase
#'  speed.
#' @param ncores number of cores
#'
#' @return a matrix of samples x features
#'
#' @examples
#'
#' data <- readRDS(
#'     system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#' )
#' feature_pathway_gsva <- run_pathway_gsva(
#'     data,
#'     geneset = NULL, species = "Homo sapiens",
#'     type = "scrna", subsample = FALSE, ncores = 1
#' )
#'
#' @importFrom msigdbr msigdbr
#' @importFrom GSVA gsva
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC getAUC
#' @import EnsDb.Hsapiens.v79
#' @import EnsDb.Mmusculus.v79
#' @importFrom dplyr select
#' @importFrom reshape2 melt
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colMeans2
#' @importFrom cli cli_abort cli_warn
#'
#' @export
run_pathway_gsva <- function(data, method = "ssgsea", geneset = NULL,
    species = "Homo sapiens",
    type = "scrna", subsample = TRUE, ncores = 1) {
    check_data(data, type)

    if (is.null(geneset)) {
        geneset <- get_geneset(species = species)
    }

    if (subsample && ncol(data) > 90000) {
        data <- data[, sample(seq_len(ncol(data)), 90000)]
    }

    if (type == "scrna") {
        # TODO: if the user does not provide geneset, need to get the geneset from msigdb
        X <- helper_pathway_gsva(
            data,
            method = method, geneset = geneset, ncores = ncores
        )
    } else if (type == "spatial_p") {
        cli::cli_warn(c(
            "Pathway GSVA currently does not support spatial proteomics"
        ))
        return(NULL)
    } else if (type == "spatial_t") {
        cli::cli_warn(c(
            "Pathway GSVA currently does not support spatial transcriptiomics"
        ))
        return(NULL)
    } else {
        cli::cli_abort(c(
            "Parameter {.var type} must be 'scrna', 'spatial_p' or 'spatial_t'",
            "x" = "'{type}' is not a valid input data type"
        ))
    }

    X <- as.data.frame(X)
    return(X)
}






#' generate pathway score using expression level
#'
#' @param data input data, a Seurat object containing `celltype` and `sample`
#'  label
#' @param geneset By default (when the `geneset` argument is not specified),
#'  we use the 50 hallmark gene set from msigdb.
#'  The users can also provide their geneset of interest in a list format, with
#'  each list entry containing a vector of the names of genes in a gene set.
#'  eg, geneset <- list("pathway_a" = c("CANS1", ...), "pathway_b" = c("PEX6"))
#' @param species whether the species is "Homo sapiens" or "Mus musculus".
#'  Default is "Homo sapiens".
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores
#'
#' @return a matrix of samples x features
#'
#' @examples
#'
#' data <- readRDS(
#'     system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#' )
#' feature_pathway_mean <- run_pathway_mean(
#'     data,
#'     geneset = NULL, species = "Homo sapiens",
#'     type = "scrna", ncores = 1
#' )
#'
#' @importFrom msigdbr msigdbr
#' @importFrom GSVA gsva
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC getAUC
#' @import EnsDb.Hsapiens.v79
#' @import EnsDb.Mmusculus.v79
#' @importFrom dplyr select
#' @importFrom reshape2 melt
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colMeans2
#' @importFrom cli cli_abort cli_warn
#'
#' @export
run_pathway_mean <- function(data, geneset = NULL,
    species = "Homo sapiens",
    type = "scrna",
    ncores = 1) {
    check_data(data, type)

    if (is.null(geneset)) {
        geneset <- get_geneset(species = species)
    }

    if (type == "scrna") {
        X <- helper_pathway_mean(data, geneset = geneset, ncores = ncores)
    } else if (type == "spatial_p") {
        cli::cli_warn(c(
            "`pathway_mean` currently does not support spatial proteomics"
        ))
        return(NULL)
    } else if (type == "spatial_t") {
        X <- helper_pathway_mean_st(data, geneset = geneset, ncores = ncores)
    } else {
        cli::cli_abort(c(
            "Parameter {.var type} must be 'scrna', 'spatial_p' or 'spatial_t'",
            "x" = "'{type}' is not a valid input data type"
        ))
    }

    X <- as.data.frame(X)
    return(X)
}





#' generate pathway score using proportion of expression
#'
#' @param data input data, a Seurat object containing `celltype` and `sample`
#'  label
#' @param geneset By default (when the `geneset` argument is not specified),
#'  we use the 50 hallmark gene set from msigdb.
#'  The users can also provide their geneset of interest in a list format, with
#'  each list entry containing a vector of the names of genes in a gene set.
#'  eg, geneset <- list("pathway_a" = c("CANS1", ...), "pathway_b" = c("PEX6"))
#' @param species whether the species is "Homo sapiens" or "Mus musculus".
#'  Default is "Homo sapiens".
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores
#'
#' @return a matrix of samples x features
#'
#' @examples
#'
#' data <- readRDS(
#'     system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#' )
#' feature_pathway_prop <- run_pathway_prop(
#'     data,
#'     geneset = NULL, species = "Homo sapiens",
#'     type = "scrna", ncores = 1
#' )
#'
#' @importFrom msigdbr msigdbr
#' @importFrom GSVA gsva
#' @importFrom AUCell AUCell_buildRankings AUCell_calcAUC getAUC
#' @import EnsDb.Hsapiens.v79
#' @import EnsDb.Mmusculus.v79
#' @importFrom dplyr select
#' @importFrom reshape2 melt
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom DelayedArray DelayedArray
#' @importFrom DelayedMatrixStats colMeans2
#' @importFrom cli cli_abort cli_warn
#'
#' @export
run_pathway_prop <- function(data, geneset = NULL,
    species = "Homo sapiens",
    type = "scrna",
    ncores = 1) {
    check_data(data, type)

    if (is.null(geneset)) {
        geneset <- get_geneset(species = species)
    }

    if (type == "scrna") {
        X <- helper_pathway_prop(data, geneset = geneset, ncores = ncores)
    } else if (type == "spatial_p") {
        cli::cli_warn(c(
            "`pathway_prop` currently does not support spatial proteomics"
        ))
        return(NULL)
    } else if (type == "spatial_t") {
        cli::cli_warn(c(
            "`pathway_prop` currently does not support spatial transcriptomics"
        ))
        return(NULL)
    } else {
        cli::cli_abort(c(
            "Parameter {.var type} must be 'scrna', 'spatial_p' or 'spatial_t'",
            "x" = "'{type}' is not a valid input data type"
        ))
    }

    X <- as.data.frame(X)
    return(X)
}








#' generate overall aggregated mean expression
#'
#' @param data input data, a Seurat object containing `celltype` and `sample`
#'  label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param genes default to NULL, in which case the top variable genes will be
#'  used
#'  If provided by user, need to be in the format of a list containing the
#'  genes of interest,
#'  eg, genes <- c(GZMA", "GZMK", "CCR7", "RPL38" )
#' @param num_top_gene when the genes is not provided by the user, the top
#'  variable genes will be used
#'  The number of genes is set by this number.
#'  default to NULL, in which case top 1500 variable genes will be selected
#' @param ncores number of cores
#'
#' @return a matrix of samples x features
#'
#' @examples
#'
#' data <- readRDS(
#'     system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#' )
#' feature_gene_mean <- run_gene_mean(
#'     data,
#'     type = "scrna", num_top_gene = 1500, ncores = 1
#' )
#'
#' @importFrom proxyC simil
#' @importFrom DelayedMatrixStats rowMeans2
#' @importFrom DelayedArray DelayedArray
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom reshape2 melt
#' @importFrom cli cli_abort
#'
#' @export
run_gene_mean <- function(data,
    type = "scrna",
    genes = NULL,
    num_top_gene = NULL,
    ncores = 1) {
    check_data(data, type)

    if (type %in% c("scrna", "spatial_p", "spatial_t")) {
        X <- helper_gene_mean(data, genes, num_top_gene, ncores)
    } else {
        cli::cli_abort(c(
            "Parameter {.var type} must be 'scrna', 'spatial_p' or 'spatial_t'",
            "x" = "'{type}' is not a valid input data type"
        ))
    }

    X <- as.data.frame(X)
    return(X)
}





#' generate overall aggregated gene proportion expression
#'
#' @param data input data, a Seurat object containing `celltype` and `sample` label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param genes default to NULL, in which case the top variable genes will be used
#' If provided by user, need to be in the format of a list containing the genes of interest,
#' eg, genes <- c(GZMA", "GZMK", "CCR7", "RPL38" )
#' @param num_top_gene when the genes is not provided by the user, the top variable genes will be used
#' The number of genes is set by this number.
#' default to NULL, in which case top 1500 variable genes will be selected
#' @param ncores number of cores
#'
#' @return a matrix of samples x features
#'
#' @examples
#'
#' data <- readRDS(system.file("extdata", "example_scrnaseq.rds", package = "scFeatures"))
#' feature_gene_prop <- run_gene_prop(data, type = "scrna", num_top_gene = 1500, ncores = 1)
#'
#' @importFrom proxyC simil
#' @importFrom DelayedMatrixStats rowMeans2
#' @importFrom DelayedArray DelayedArray
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom reshape2 melt
#' @importFrom cli cli_abort
#'
#' @export
run_gene_prop <- function(data, type = "scrna", genes = NULL, num_top_gene = NULL, ncores = 1) {
    check_data(data, type)

    if (type %in% c("scrna", "spatial_p", "spatial_t")) {
        X <- helper_gene_prop(data, genes, num_top_gene, ncores)
    } else {
        cli::cli_abort(c(
            "Parameter {.var type} must be 'scrna', 'spatial_p' or 'spatial_t'",
            "x" = "'{type}' is not a valid input data type"
        ))
    }

    X <- as.data.frame(X)
    return(X)
}





#' generate overall aggregated gene correlation
#'
#' @param data input data, a Seurat object containing `celltype` and `sample` label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param genes default to NULL, in which case the top variable genes will be used
#' If provided by user, need to be in the format of a list containing the feature of interest,
#' eg, genes <- c(GZMA", "GZMK", "CCR7", "RPL38" )
#' @param num_top_gene when the genes is not provided by the user, the top variable genes will be used
#' The number of genes is set by this number.
#' default to NULL, in which case top 50 variable genes will be selected
#' @param ncores number of cores
#'
#' @return a matrix of samples x features
#'
#' @examples
#'
#' data <- readRDS(system.file("extdata", "example_scrnaseq.rds", package = "scFeatures"))
#' feature_gene_cor <- run_gene_cor(data, type = "scrna", num_top_gene = 1500, ncores = 1)
#'
#' @importFrom proxyC simil
#' @importFrom DelayedMatrixStats rowMeans2
#' @importFrom DelayedArray DelayedArray
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom reshape2 melt
#' @importFrom cli cli_abort
#'
#' @export
run_gene_cor <- function(data, type = "scrna", genes = NULL, num_top_gene = NULL, ncores = 1) {
    check_data(data, type)

    if (type %in% c("scrna", "spatial_p", "spatial_t")) {
        X <- helper_gene_cor(data, genes, num_top_gene, ncores)
    } else {
        cli::cli_abort(c(
            "Parameter {.var type} must be 'scrna', 'spatial_p' or 'spatial_t'",
            "x" = "'{type}' is not a valid input data type"
        ))
    }

    X <- as.data.frame(X)
    return(X)
}




#' generate L stats
#'
#' @param data input data, a Seurat object containing `celltype` and `sample`
#'  label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores
#'
#' @return a matrix of samples x features
#'
#' @examples
#'
#' data <- readRDS(system.file(
#'     "extdata", "example_spatial_proteomics.rds",
#'     package = "scFeatures"
#' ))
#' feature_L_function <- run_L_function(data, type = "spatial_p", ncores = 1)
#'
#' @importFrom spatstat.geom ppp pairdist
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom reshape2 melt
#' @importFrom ape Moran.I
#' @importFrom spatstat.explore nncorr
#' @importFrom cli cli_abort cli_warn
#'
#' @export
run_L_function <- function(data, type = "spatial_p", ncores = 1) {
    check_data(data, type)

    if (type == "scrna") {
        cli::cli_warn(c(
            "`L_function` currently does not support 'scrna'"
        ))
        return(NULL)
    } else if (type == "spatial_p") {
        X <- helper_L_stat_sp(data)
    } else if (type == "spatial_t") {
        X <- helper_L_stat_st(data)
    } else {
        cli::cli_abort(c(
            "Parameter {.var type} must be 'scrna', 'spatial_p' or 'spatial_t'",
            "x" = "'{type}' is not a valid input data type"
        ))
    }

    X <- as.data.frame(X)
    return(X)
}





#' generate cell type interaction
#'
#' @param data input data, a Seurat object containing `celltype` and `sample`
#'  label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores
#'
#' @return a matrix of samples x features
#'
#' @examples
#'
#' data <- readRDS(system.file(
#'     "extdata", "example_spatial_proteomics.rds",
#'     package = "scFeatures"
#' ))
#' feature_celltype_interaction <- run_celltype_interaction(
#'     data,
#'     type = "spatial_p", ncores = 1
#' )
#'
#' @importFrom spatstat.geom ppp pairdist
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom reshape2 melt
#' @importFrom ape Moran.I
#' @importFrom spatstat.explore nncorr
#' @importFrom cli cli_abort cli_warn
#'
#' @export
run_celltype_interaction <- function(data, type = "spatial_p", ncores = 1) {
    check_data(data, type)

    if (type == "scrna") {
        cli::cli_warn(c(
            "`celltype_interaction` currently does not support 'scrna'"
        ))
        return(NULL)
    } else if (type == "spatial_p") {
        X <- helper_celltype_interaction_sp(data)
    } else if (type == "spatial_t") {
        X <- helper_celltype_interaction_st(data)
    }

    X <- as.data.frame(X)
    return(X)
}








#' generate Moran's I
#'
#' @param data input data, a Seurat object containing `celltype` and `sample`
#'  label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores
#'
#' @return a matrix of samples x features
#'
#' @examples
#'
#' data <- readRDS(system.file(
#'     "extdata", "example_spatial_proteomics.rds",
#'     package = "scFeatures"
#' ))
#' feature_Morans_I <- run_Morans_I(data, type = "spatial_p", ncores = 1)
#'
#' @importFrom spatstat.geom ppp pairdist
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom reshape2 melt
#' @importFrom ape Moran.I
#' @importFrom spatstat.explore nncorr
#' @importFrom cli cli_abort cli_warn
#'
#' @export
run_Morans_I <- function(data, type = "spatial_p", ncores = 1) {
    check_data(data, type)

    if (type == "scrna") {
        cli::cli_warn(c(
            "`Morans_I` currently does not support 'scrna'"
        ))
        return(NULL)
    } else if (type %in% c("spatial_p", "spatial_t")) {
        X <- helper_moran(data, num_top_gene = NULL, ncores = 1)
    } else {
        cli::cli_abort(c(
            "Parameter {.var type} must be 'scrna', 'spatial_p' or 'spatial_t'",
            "x" = "'{type}' is not a valid input data type"
        ))
    }

    X <- as.data.frame(X)
    return(X)
}









#' generate nearest neighbour correlation
#'
#' @param data input data, a Seurat object containing `celltype` and `sample`
#'  label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores
#'
#' @return a matrix of samples x features
#'
#' @examples
#'
#' data <- readRDS(system.file(
#'     "extdata", "example_spatial_proteomics.rds",
#'     package = "scFeatures"
#' ))
#' feature_nn_correlation <- run_nn_correlation(
#'     data,
#'     type = "spatial_p", ncores = 1
#' )
#'
#' @importFrom spatstat.geom ppp pairdist
#' @importFrom BiocParallel SerialParam bplapply
#' @importFrom reshape2 melt
#' @importFrom ape Moran.I
#' @importFrom spatstat.explore nncorr
#' @importFrom cli cli_abort cli_warn
#'
#' @export
run_nn_correlation <- function(data, type = "spatial_p", ncores = 1) {
    check_data(data, type)

    if (type == "scrna") {
        cli::cli_warn(c(
            "`nn_correlation` currently does not support 'scrna'"
        ))
        return(NULL)
    } else if (type %in% c("spatial_p", "spatial_t")) {
        X <- helper_nncorr_protein(data, num_top_gene = NULL, ncores = ncores)
    } else {
        cli::cli_abort(c(
            "Parameter {.var type} must be 'scrna', 'spatial_p' or 'spatial_t'",
            "x" = "'{type}' is not a valid input data type"
        ))
    }

    X <- as.data.frame(X)
    return(X)
}
