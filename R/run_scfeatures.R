# this file contains the 17 functions that generate the 17 feature classes



#' Generate cell type proportion raw
#' 
#' @description 
#' This function calculates the proportions of cells belonging to each cell type. 
#' The input data must contain `sample` and `celltype` metadata column. 
#' The function supports scRNA-seq and spatial proteomics.
#' The function returns a dataframe with samples as rows and cell types as columns.
#' 
#' @param data A Seurat object containing `celltype` and `sample` label
#' @param type The type of dataset, either "scrna", "spatial_t", or "spatial_p".
#' @param ncores Number of cores for parallel processing.
#'
#' @return a dataframe of samples x features
#'
#' @examples
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




#' Generate cell type proportions, with logit transformation
#' 
#' @description
#' This function calculates the proportions of cells belonging to each cell type,
#' and applies a logit transformation to the proportions. 
#' The input data must contain `sample` and `celltype` metadata column. 
#' The function supports scRNA-seq and spatial proteomics.
#' The function returns a dataframe with samples as rows and cell types as columns.
#' 
#' @param data A Seurat object containing `celltype` and `sample` label
#' @param type The type of dataset, either "scrna", "spatial_t", or "spatial_p".
#' @param ncores Number of cores for parallel processing.
#'
#' @return a dataframe of samples x features
#'
#' @examples
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



#' Generate cell type proportion ratio
#' 
#' @description
#' This function calculates pairwise cell type proportion ratio for each sample. 
#' and applies a logit transformation to the proportions. 
#' The input data must contain `sample` and `celltype` metadata column. 
#' The function supports scRNA-seq and spatial proteomics.
#' The function returns a dataframe with samples as rows and cell types as columns.
#' 
#' @param data A Seurat object containing `celltype` and `sample` label
#' @param type The type of dataset, either "scrna", "spatial_t", or "spatial_p".
#' @param ncores Number of cores for parallel processing.
#'
#' @return a dataframe of samples x features
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


 

#' Generate cell type specific gene mean expression
#' 
#' @description 
#' This function computes the mean expression of a set of genes for
#' each cell type in the input data. The input data can be of three types: 
#' 'scrna', 'spatial_p' or 'spatial_t'. If the genes parameter is not p
#' rovided by the user, the top variable genes will be selected based on 
#' the num_top_gene parameter (defaults to 100). 
#' 
#' @param data A Seurat object containing `celltype` and `sample` label
#' @param type The type of dataset, either "scrna", "spatial_t", or "spatial_p".
#' @param genes Optional dataframe with 2 columns: 'marker' and 'celltype'.
#' The 'marker' column should contain the genes of interest (e.g. 'S100A11', 'CCL4'), 
#' and the 'celltype' column should contain the celltype that the gene expression 
#' is to be computed from (e.g. 'CD8', 'B cells').
#' If not provided, the top variable genes will be used based on the 
#' num_top_gene parameter.
#' @param num_top_gene Number of top genes to use when genes is not provided.
#' Defaults to 100.
#' @param ncores Number of cores for parallel processing.
#'
#' @return a dataframe of samples x features
#'
#' @examples
#'
#'   data <- readRDS(
#'     system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#'   )
#' 
#'   # optional step, if mito and ribo genes are not of interest
#'   data_remove_mito <- remove_mito(data)
#' 
#'   feature_gene_mean_celltype <- run_gene_mean_celltype(
#'     data_remove_mito,
#'     type = "scrna", num_top_gene = 100, ncores = 1
#'   )
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
        X <- helper_gene_mean_celltype_st(data, genes, num_top_gene, ncores)
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

#' Generate cell type specific gene proportion expression
#'
#' This function computes the proportion of expression of a set of genes 
#' for each cell type in the input data. The input data can be of three types: 
#' 'scrna', 'spatial_p' or 'spatial_t'. If the genes parameter is not provided 
#' by the user, the top variable genes will be selected based on the 
#' num_top_gene parameter (defaults to 100). 
#' 
#' @param data A Seurat object containing `celltype` and `sample` label
#' @param type The type of dataset, either "scrna", "spatial_t", or "spatial_p".
#' @param genes Optional dataframe with 2 columns: 'marker' and 'celltype'.
#' The 'marker' column should contain the genes of interest (e.g. 'S100A11', 'CCL4'), 
#' and the 'celltype' column should contain the celltype that the gene expression 
#' is to be computed from (e.g. 'CD8', 'B cells').
#' If not provided, the top variable genes will be used based on the 
#' num_top_gene parameter.
#' @param num_top_gene Number of top genes to use when genes is not provided.
#' Defaults to 100.
#' @param ncores Number of cores for parallel processing.
#'
#' @return a dataframe of samples x features
#'
#' @examples
#'
#'  data <- readRDS(
#'     system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#'  )
#' 
#' # optional step, if mito and ribo genes are not of interest
#' data_remove_mito <- remove_mito(data)
#' 
#'  feature_gene_prop_celltype <- run_gene_prop_celltype(
#'     data_remove_mito,
#'     type = "scrna", num_top_gene = 100, ncores = 1
#'  )
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






#' Generate cell type specific gene expression correlation
#'
#' This function computes the correlation of expression of a set of genes 
#' for each cell type in the input data. The input data can be of three types: 
#' 'scrna', 'spatial_p' or 'spatial_t'. If the genes parameter is not provided 
#' by the user, the top variable genes will be selected based on the 
#' num_top_gene parameter (defaults to 100). 
#' 
#' @param data A Seurat object containing `celltype` and `sample` label
#' @param type The type of dataset, either "scrna", "spatial_t", or "spatial_p".
#' @param genes Optional dataframe with 2 columns: 'marker' and 'celltype'.
#' The 'marker' column should contain the genes of interest (e.g. 'S100A11', 'CCL4'), 
#' and the 'celltype' column should contain the celltype that the gene expression 
#' is to be computed from (e.g. 'CD8', 'B cells').
#' If not provided, the top variable genes will be used based on the 
#' num_top_gene parameter.
#' @param num_top_gene Number of top genes to use when genes is not provided.
#' Defaults to 5.
#' @param ncores Number of cores for parallel processing.
#'
#' @return a dataframe of samples x features
#'
#' @examples
#'
#'  data <- readRDS(system.file("extdata",
#'    "example_scrnaseq.rds",
#'     package = "scFeatures"
#'  ))
#' 
#'  # optional step, if mito and ribo genes are not of interest
#' 
#'  data_remove_mito <- remove_mito(data)
#' 
#'  feature_gene_cor_celltype <- run_gene_cor_celltype(
#'    data_remove_mito,
#'    type = "scrna", num_top_gene = 100, ncores = 1
#'  )
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






#' Generate pathway score using gene set enrichement analysis
#'
#' @description 
#' This function calculates pathway scores for a given input 
#' dataset and gene set using gene set enrichment analysis (GSVA). 
#' It supports scRNA-seq, spatial proteomics and spatial transcriptomics. 
#' It currently supports two pathway analysis methods: ssgsea and aucell. 
#' By default, it uses the 50 hallmark gene sets from msigdb. 
#' Alternatively, users can provide their own gene sets of interest 
#' in a list format.
#' 
#' @param data A Seurat object containing `celltype` and `sample` label
#' @param type The type of dataset, either "scrna", "spatial_t", or "spatial_p".
#' @param method Type of pathway analysis method, currently support `ssgsea`
#'  and `aucell`
#' @param geneset By default (when the `geneset` argument is not specified),
#'  we use the 50 hallmark gene set from msigdb.
#'  The users can also provide their geneset of interest in a list format, with
#'  each list entry containing a vector of the names of genes in a gene set.
#'  eg, geneset <- list("pathway_a" = c("CAPN1", ...), "pathway_b" = c("PEX6"))
#' @param species Whether the species is "Homo sapiens" or "Mus musculus".
#'  Default is "Homo sapiens".
#' @param subsample Whether to subsample, either T or F. For larger datasets
#'  (eg, over 30,000 cells), the subsample function can be used to increase
#'  speed.
#' @param ncores Number of cores for parallel processing.
#'
#' @return a dataframe of samples x features
#'
#' @examples
#' 
#'  data <- readRDS(
#'     system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#'  )
#'  feature_pathway_gsva <- run_pathway_gsva(
#'     data,
#'     geneset = NULL, species = "Homo sapiens",
#'     type = "scrna", subsample = FALSE, ncores = 1
#'  )
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
        capture.output( suppressMessages( X <- helper_pathway_gsva(
            data,
            method = method, geneset = geneset, ncores = ncores
        ) ))
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






#' Generate pathway score using expression level
#' 
#' @description 
#' This function calculates pathway scores for a given dataset and gene set 
#' using gene expression levels. It supports scRNA-seq, spatial transcriptomics 
#' and spatial proteomics and spatial transcriptomics). 
#' By default, it uses the 50 hallmark gene sets from msigdb. 
#' Alternatively, users can provide their own gene sets of interest in a list format.
#' 
#' @param data A Seurat object containing `celltype` and `sample` label
#' @param type The type of dataset, either "scrna", "spatial_t", or "spatial_p".
#' @param geneset By default (when the `geneset` argument is not specified),
#'  we use the 50 hallmark gene set from msigdb.
#'  The users can also provide their geneset of interest in a list format, with
#'  each list entry containing a vector of the names of genes in a gene set.
#'  eg, geneset <- list("pathway_a" = c("CANS1", ...), "pathway_b" = c("PEX6"))
#' @param species Whether the species is "Homo sapiens" or "Mus musculus".
#'  Default is "Homo sapiens".
#' @param ncores Number of cores for parallel processing.
#'
#' @return a dataframe of samples x features
#'
#' @examples
#' 
#'  data <- readRDS(
#'     system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#'  )
#'  feature_pathway_mean <- run_pathway_mean(
#'     data,
#'     geneset = NULL, species = "Homo sapiens",
#'     type = "scrna", ncores = 1
#'  )
#' 
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





#' Generate pathway score using proportion of expression
#' 
#' @description 
#' This function calculates pathway scores for a given input dataset and gene set 
#' using the proportion of gene expression levels. It supports scRNA-seq, spatial transcriptomics 
#' and spatial proteomics and spatial transcriptomics). 
#' By default, it uses the 50 hallmark gene sets from msigdb. 
#' Alternatively, users can provide their own gene sets of interest in a list format.
#' 
#' @param data A Seurat object containing `celltype` and `sample` label
#' @param type The type of dataset, either "scrna", "spatial_t", or "spatial_p".
#' @param geneset By default (when the `geneset` argument is not specified),
#'  we use the 50 hallmark gene set from msigdb.
#'  The users can also provide their geneset of interest in a list format, with
#'  each list entry containing a vector of the names of genes in a gene set.
#'  eg, geneset <- list("pathway_a" = c("CANS1", ...), "pathway_b" = c("PEX6"))
#' @param species Whether the species is "Homo sapiens" or "Mus musculus".
#'  Default is "Homo sapiens".
#' @param ncores Number of cores for parallel processing.
#'
#' @return a dataframe of samples x features
#'
#' @examples
#' 
#'  data <- readRDS(
#'     system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#'  )
#'  feature_pathway_prop <- run_pathway_prop(
#'     data,
#'     geneset = NULL, species = "Homo sapiens",
#'     type = "scrna", ncores = 1
#'  )
#' 
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






#' Generate cell cell communication score
#' 
#' @description 
#' This function calculates the ligand receptor interaction score using SingleCellSignalR. 
#' The output features are in the form of celltype a -> celltype b -- ligand 1 -> receptor 2 ,
#' which indicates the interaction between ligand 1 in celltype a and receptor 2 in celltype b. 
#' 
#' It supports scRNA-seq. 
#' 
#' @param data input data, a Seurat object containing `celltype` and `sample` label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @examples
#'
#' data <- readRDS(system.file("extdata", "example_scrnaseq.rds", package = "scFeatures"))
#' feature_CCI <- run_CCI(data, type = "scrna" ,  ncores = 1 )
#' 
#' @import dplyr
#' @import DelayedArray
#' @import SingleCellSignalR
#' 
#' @export
run_CCI <- function( data, type = "scrna" , ncores = 1  ){
  
  check_data(data, type)
  
  if ( type == "scrna" )  {
      X <- helper_CCI(data, ncores =  ncores )
  }
  
  if ( type == "spatial_p" )  {
    print("This feature class currently does not support spatial proteomics")
    return(NULL)
  }
  
  if ( type == "spatial_t" ) {
    print("This feature class currently does not support spatial transcriptomics")
    return(NULL)
  }
  
  return (X)
  
}







#' Generate overall aggregated mean expression
#'
#' @description 
#' This function computes the mean expression of genes across samples. The user
#' can specify the genes of interest, or let the function use the top variable
#' genes by default. The function supports scRNA-seq, spatial proteomics, 
#' and spatial transcriptomics.
#'
#' @param data A Seurat object containing `celltype` and `sample` label
#' @param type The type of dataset, either "scrna", "spatial_t", or "spatial_p".
#' @param genes Default to NULL, in which case the top variable genes will be
#'  used. If provided by user, need to be in the format of a list containing the
#'  genes of interest, eg, genes <- c(GZMA", "GZMK", "CCR7", "RPL38" )
#' @param num_top_gene Number of top variable genes to use when genes is not provided.
#' Defaults to 1500.
#' @param ncores Number of cores for parallel processing.
#'
#' @return a dataframe of samples x features
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





#' Generate overall aggregated gene proportion expression
#'
#' @description 
#' This function computes the proportion of gene expression across samples. The user
#' can specify the genes of interest, or let the function use the top variable
#' genes by default. The function supports scRNA-seq, spatial proteomics, 
#' and spatial transcriptomics.
#'
#' @param data A Seurat object containing `celltype` and `sample` label
#' @param type The type of dataset, either "scrna", "spatial_t", or "spatial_p".
#' @param genes Default to NULL, in which case the top variable genes will be
#'  used. If provided by user, need to be in the format of a list containing the
#'  genes of interest, eg, genes <- c(GZMA", "GZMK", "CCR7", "RPL38" )
#' @param num_top_gene Number of top variable genes to use when genes is not provided.
#' Defaults to 1500.
#' @param ncores Number of cores for parallel processing.
#' 
#' @return a dataframe of samples x features
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





#' Generate overall aggregated gene correlation
#'
#' @description 
#' This function computes the correlation of gene expression across samples. The user
#' can specify the genes of interest, or let the function use the top variable
#' genes by default. The function supports scRNA-seq, spatial proteomics, 
#' and spatial transcriptomics.
#'
#' @param data A Seurat object containing `celltype` and `sample` label
#' @param type The type of dataset, either "scrna", "spatial_t", or "spatial_p".
#' @param genes Default to NULL, in which case the top variable genes will be
#'  used. If provided by user, need to be in the format of a list containing the
#'  genes of interest, eg, genes <- c(GZMA", "GZMK", "CCR7", "RPL38" )
#' @param num_top_gene Number of top variable genes to use when genes is not provided.
#' Defaults to 5.
#' @param ncores Number of cores for parallel processing.
#'
#' @return a dataframe of samples x features
#'
#' @examples
#'
#'  data <- readRDS(
#'    system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#'  )
#'  feature_gene_cor <- run_gene_cor(
#'    data, type = "scrna", num_top_gene = 5, ncores = 1
#'  )
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




#' Generate L stats
#' 
#' @description
#' This function calculates L-statistics to measure spatial autocorrelation.
#' The function supports  spatial proteomics and spatial transcriptomics.
#' 
#' @param data A Seurat object containing `celltype` and `sample` label
#' @param type The type of dataset, either "scrna", "spatial_t", or "spatial_p".
#' @param ncores Number of cores for parallel processing.
#'
#' @return a dataframe of samples x features
#'
#' @examples
#'
#' data <- readRDS(
#' system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#' )
#' x <- sample(1:100, ncol(data) , replace = TRUE)
#' y <- sample(1:100, ncol(data) , replace = TRUE)
#' data <- makeSeurat(data, spatialCoords = list(x,y))
#' data$sample <- sample( c("patient1", "patient2", "patient3"), ncol(data) , replace= TRUE )
#' 
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





#' Generate cell type interaction
#'
#' @description 
#' This function calculates the pairwise distance between cell types
#' for a sample by using the coordinates and cell types of the cells. 
#' The function supports spatial proteomics and spatial transcriptomics.
#' 
#' @param data A Seurat object containing `celltype` and `sample` label
#' @param type The type of dataset, either "scrna", "spatial_t", or "spatial_p".
#' @param ncores Number of cores for parallel processing.
#'
#' @return a dataframe of samples x features
#'
#' @examples
#'
#' 
#'data <- readRDS(
#' system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#' )
#' x <- sample(1:100, ncol(data) , replace = TRUE)
#' y <- sample(1:100, ncol(data) , replace = TRUE)
#' data <- makeSeurat(data, spatialCoords = list(x,y))
#' data$sample <- sample( c("patient1", "patient2", "patient3"), ncol(data) , replace= TRUE )
#' 
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








#' Generate Moran's I
#' 
#' @description 
#' This function calculates Moran's I to measure spatial autocorrelation,
#' which an indicattion of how strongly the feature(ie, genes/proteins) 
#' expression values in a sample cluster or disperse. The function supports 
#' spatial proteomics and spatial transcriptomics.
#'  
#' @param data A Seurat object containing `celltype` and `sample` label
#' @param type The type of dataset, either "scrna", "spatial_t", or "spatial_p".
#' @param ncores Number of cores for parallel processing.
#'
#' @return a dataframe of samples x features
#'
#' @examples
#'
#'data <- readRDS(
#' system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#' )
#' x <- sample(1:100, ncol(data) , replace = TRUE)
#' y <- sample(1:100, ncol(data) , replace = TRUE)
#' data <- makeSeurat(data, spatialCoords = list(x,y))
#' data$sample <- sample( c("patient1", "patient2", "patient3"), ncol(data) , replace= TRUE )
#' 
#' feature_celltype_interaction <- run_celltype_interaction(
#'     data,
#'     type = "spatial_p", ncores = 1
#' )
#'
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









#' Generate nearest neighbour correlation
#' 
#' @description 
#' This function calculates the nearest neighbour correlation for each feature 
#' (eg, proteins) in each sample. This is calculated by taking the correlation 
#' between each cell and its nearest neighbours cell for a particular feature.  
#' This function supports spatial proteomics, and spatial transcriptomics.
#'
#' @param data A Seurat object containing `celltype` and `sample` label
#' @param type The type of dataset, either "scrna", "spatial_t", or "spatial_p".
#' @param ncores Number of cores for parallel processing.
#' @param num_top_gene Number of top variable genes to use when genes is 
#' not provided. Defaults to 1500.
#' 
#' @return a dataframe of samples x features
#'
#' @examples
#'
#' data <- readRDS(
#' system.file("extdata", "example_scrnaseq.rds", package = "scFeatures")
#' )
#' x <- sample(1:100, ncol(data) , replace = TRUE)
#' y <- sample(1:100, ncol(data) , replace = TRUE)
#' data <- makeSeurat(data, spatialCoords = list(x,y))
#' data$sample <- sample( c("patient1", "patient2", "patient3"), ncol(data) , replace= TRUE )
#' 
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
run_nn_correlation <- function(data, type = "spatial_p", num_top_gene = NULL, ncores = 1) {
    check_data(data, type)

    if (type == "scrna") {
        cli::cli_warn(c(
            "`nn_correlation` currently does not support 'scrna'"
        ))
        return(NULL)
    } else if (type %in% c("spatial_p", "spatial_t")) {
        X <- helper_nncorr_protein(data, num_top_gene, ncores = ncores)
    } else {
        cli::cli_abort(c(
            "Parameter {.var type} must be 'scrna', 'spatial_p' or 'spatial_t'",
            "x" = "'{type}' is not a valid input data type"
        ))
    }

    X <- as.data.frame(X)
    return(X)
}
