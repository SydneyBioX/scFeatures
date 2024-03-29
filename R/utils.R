

#' Check if required metadata is in the object
#'
#' @description 
#' This function checks that the object contains all the 
#' necessary metadata for functions in scFeatures. 
#' For scRNA-seq, it must contain `sample` and `celltype`.
#' For spatial proteomics, it must contain `sample`, `celltype` and 
#' and spatial proteomics data, it must contain `sample` and `x_cord`
#' and `y_cord`. 
#' For spatial transcriptomics data, it must contain `sample`, `x_cord`
#' `y_cord` and a `predictions` assay.
#' 
#' @param alldata A list object. 
#' @param type Type of dataset, either "scrna" (for scRNA-seq), "spatial_t" 
#' (for spatial transcriptomics) or "spatial_p" (for spatial proteomics). 
#'
#' @return NULL
#' @importFrom methods is
#' 
#' @noRd
check_data <- function(alldata, type = "scrna") {
  
    if ( !"data" %in% names(alldata)) {
        stop("Please make sure you provide the data")
    }

    if (!"sample" %in% names(alldata)) {
        stop(
            "For scRNA-seq and spatial proteomics, ",
            "please make sure the data contains celltype and sample label. ",
            "For spatial proteomics, ensure the data contain sample information."
        )
    }

    if (type %in% c("spatial_t", "spatial_p")) {
        if (!"y_cord" %in% names(alldata) || !"x_cord" %in% names(alldata)) {
            stop("Please ensure the data contain x_cord and y_cord.")
        }
    }

    if (type == "spatial_t") {
        if (!"predictions" %in% names(alldata)) {
            stop(
                "Please make sure the data contain a predictions assay.\n",
                "See vignette's section on spatial transcriptomics for details."
            )
        }
    }
}





#' Enable parallel processing
#'
#' @description 
#' This function takes the number of cores to use for parallel processing and
#' generates a `BiocParallel` object that can be used to control the 
#' parallelization of functions. 
#' It automatically determines whether to use the `SnowParam` or 
#' `MulticoreParam` of the `BiocParallel` package based on the operating system.
#'
#' @param cores The number of cores to use for parallel processing.
#' @return A `BiocParallel` object that can be used to control the 
#' parallelization of functions.
#' 
#' @noRd
generateBPParam <- function(cores = 1) {
    seed <- .Random.seed[1]

    if (cores == 1) {
        BPparam <- BiocParallel::SerialParam(RNGseed = seed)
    } else { # Parallel processing is desired.
        # set BPparam RNGseed if the user ran set.seed(someNumber) themselves.
        if (Sys.info()["sysname"] == "Windows") { # Only SnowParam suits Windows.
            BPparam <- BiocParallel::SnowParam(
                min(cores, BiocParallel::snowWorkers("SOCK")),
                RNGseed = seed
            )
        } else if (Sys.info()["sysname"] %in% c("MacOS", "Linux")) {
            BPparam <- BiocParallel::MulticoreParam(
                min(cores, BiocParallel::multicoreWorkers()),
                RNGseed = seed
            ) # Multicore is faster than SNOW, but it doesn't work on Windows.
        } else { # Something weird.
            BPparam <- BiocParallel::bpparam() # BiocParallel will figure it out.
        }
    }

    return(BPparam)
}





#' Create pseudo-bulk for each cell type of each sample
#'
#' @description 
#' This function takes a list object as input and creates a pseudo-bulk 
#' for each cell type of each sample in the object. This is computed by
#' taking the row means of the expression for each cell type of each sample 
#' If a cell type does not exist in a sample, the expression values for that 
#' cell type will be 0 for all genes.
#'
#' @param alldata A list object. 
#' @param ncores Number of cores for parallel computation.
#' @return A matrix containing the pseudo-bulks for each cell type
#'  of each sample. 
#' 
#' @noRd
bulk_sample_celltype <- function(alldata, ncores = 1) {
  
    BPparam <- generateBPParam(ncores)

    # x <- unique( alldata$sample)[1]
    bulk <- BiocParallel::bplapply(unique(alldata$sample), function(x) {
        # for this patient
        this_sample <- alldata$data[, alldata$sample == x, drop=FALSE]
        this_sample_celltype <- alldata$celltype[alldata$sample == x, drop=FALSE]

        # loop through each cell type
        # y <-unique(alldata$celltype)[1]
        this_sample_bulk <- lapply(unique(alldata$celltype), function(y) {
            index <- which(this_sample_celltype == y)

            # if cell type does not exist in patient, expression is 0 for all genes
            if (length(index) == 0) {
                temp <- rep(0, nrow(alldata$data))
                # if there is only 1 cell, do not need to take mean
            } else if (length(index) == 1) {
                temp <- this_sample[, index]
                # if multiple cells, average across all cells
            } else {
                temp <- DelayedMatrixStats::rowMeans2(
                    DelayedArray(this_sample[, index])
                )
            }

            temp <- as.matrix(temp)
            rownames(temp) <- rownames(alldata$data)
            temp
        })

        this_sample_bulk <- as.data.frame(do.call(cbind, this_sample_bulk))

        colnames(this_sample_bulk) <- paste0(x, "--", unique(alldata$celltype))
        this_sample_bulk
    }, BPPARAM = BPparam)

    bulk <- as.data.frame(do.call(cbind, bulk))

    return(bulk)
}



#' Create pseudo-bulk for each sample in a list object
#'
#' @description  
#' This function takes a list object as input and creates 
#' a pseudo-bulk for each sample in the object. This is calculated by
#' taking row means of the expression for each sample. 
#'
#' @param alldata A list object. 
#' @param ncores Number of cores for parallel computation. 
#' @return A matrix containing the pseudo-bulks for each sample 
#' in the input data. 
#' 
#' @noRd
bulk_sample <- function(alldata, ncores = 1) {
    BPparam <- generateBPParam(ncores)

    
    # x <- unique(alldata$sample)[1]
    bulk <- BiocParallel::bplapply( unique(alldata$sample) , function(x) {
        index <- which( alldata$sample == x)
        temp <- DelayedMatrixStats::rowMeans2(
            DelayedArray(alldata$data[, index, drop=FALSE])
        )
        temp <- as.matrix(temp)
    }, BPPARAM = BPparam)

    bulk <- as.data.frame(do.call(cbind, bulk))
    rownames(bulk) <- rownames(alldata$data)
    colnames(bulk) <- unique(alldata$sample)
  
   
    return(bulk)
}



#' Estimate a relative number of cells per spot for
#' spatial transcriptomics data
#'
#' This function takes a list object containing spatial transcriptomics matrix as input 
#' and estimates the relative number of cells per spot in the data. 
#' The number of cells is estimated as the library size scaled to 
#' the range from 1 to 100.
#' This value stored in the `number_cells` attribute.
#' 
#' @param alldata A list object containing spatial transcriptomics
#'
#' @return a vector with the relative number of cells in each spot. 
#' 
#' @importFrom MatrixGenerics colSums2
#'
#' @examples
#'
#' utils::data("example_scrnaseq" , package = "scFeatures")
#' data <- example_scrnaseq@assays$RNA@data
#' data <- list(data = data)
#' number_of_cells <- get_num_cell_per_spot(data)
#'
#'
#' @export
get_num_cell_per_spot <- function(alldata) {
 
  readcount <- log2(colSums(alldata$data))
  
  linMap <- function(x, from, to) {
    (x - min(x)) / max(x - min(x)) * (to - from) + from
  }
  
  numberofcells <- linMap(readcount, 1, 100)
 
  return(numberofcells)
}


#' Rearrange string
#' 
#' @description 
#' This function takes a string as input and rearranges words in the string 
#' so that they are sorted alphabetically. 
#'
#' @param str A character string containing words separated by underscores.
#' 
#' @return A character string with the words sorted alphabetically and
#' separated by underscores.
#' 
#' @noRd
rearrange_string <- function(str) {
    unlist(lapply(strsplit(str, "_"), function(x) paste(sort(x), collapse = "_")))
}



#' Compute number of cells in each cell type in each spot 
#' for spatial transcriptomics data
#' 
#' @description 
#' This function computes the number of cells in each cell type by 
#' multiplying the cell type probability in each spot with the 
#' relative number of cells in each spot. The relative number of cells
#' are estimated using library size of each spot. See get_num_cell_per_spot().
#' 
#' @param alldata A list object containing spatial transcriptomics. 
#' 
#' @return A matrix with the number of cells per cell type at each spot.
#' 
#' @noRd
get_num_cell_per_celltype <- function(alldata) {
  
  
    number_cells <-  get_num_cell_per_spot(alldata)
    alldata$number_cells <- number_cells
    
    prob <- as.matrix(alldata$predictions)
    zero_celltype <- names(which(rowSums(prob) == 0))
    prob <- prob[!rownames(prob) %in% zero_celltype, ]

    MultVec <- alldata$number_cells
    num_cell_per_spot <- mapply(FUN = "*", as.data.frame(prob), MultVec)

    num_cell_per_spot <- round(num_cell_per_spot, 0)
    mode(num_cell_per_spot) <- "integer"

    rownames(num_cell_per_spot) <- rownames(prob)

    return(num_cell_per_spot)
}



#' Compute L statistic for a point pattern
#'
#' @description
#' This function computes the L statistic for a given point pattern. 
#'
#' @param ppp_obj a point pattern object
#' @param from define the window of the point pattern to compute
#' @param to define the window of the point pattern to compute
#' @param L_dist a numeric value specifying the maximum distance 
#' at which the L statistic will be computed.
#'
#' @return
#' A numeric value of the L statistic.
#'
#' @noRd
L_stats <- function(ppp_obj = NULL, from = NULL, to = NULL, L_dist = NULL) {
    L <- spatstat.explore::Lcross(ppp_obj,
        from = from, to = to,
        verbose = FALSE,
        correction = "best"
    )

    L_theo <- L$theo[L$r <= L_dist]
    L_iso <- L$iso[L$r <= L_dist]
    L_res <- mean(L_iso - L_theo)

    return(L_res)
}



# TODO: check if return value is html file.

#' Create an association study report in HTML format
#'
#' @description This function takes the feature matrix generated by
#' scFeatures as input and creates an HTML report containing the
#' results of the association study. The report is saved to the specified 
#' output folder.
#' 
#' @param scfeatures_result a named list storing the scFeatures feature output.
#' Note that the names of the list should be one or multiple of the following:
#' `proportion_raw`, `proportion_logit`, `proportion_ratio`,
#' `gene_mean_celltype`, `gene_prop_celltype`, `gene_cor_celltype`,
#'  `pathway_gsva`, `pathway_mean`, `pathway_prop`, `CCI`,
#' `gene_mean_aggregated`, `gene_cor_aggregated`, and `gene_prop_aggregated`.
#' 
#' @param output_folder the path to the folder where the HTML report 
#' will be saved
#'
#' @return an HTML file, saved to the directory defined in the `output_folder`
#' argument
#' 
#' @examples
#' \dontrun{
#' output_folder <- tempdir()
#' utils::data("scfeatures_result" , package = "scFeatures")
#' run_association_study_report(scfeatures_result, output_folder )
#' }
#' @import rmarkdown
#' 
#' @export
run_association_study_report <- function(scfeatures_result, output_folder) {
    # check name
    correct_name <- any(names(scfeatures_result) %in% c(
        "proportion_raw", "proportion_logit",
        "proportion_ratio", "gene_mean_celltype",
        "gene_prop_celltype", "gene_cor_celltype",
        "pathway_gsva", "pathway_mean",
        "pathway_prop", "CCI",
        "gene_mean_aggregated", "gene_cor_aggregated", "gene_prop_aggregated"
    ))
    if (correct_name == FALSE) {
        warning(
            "Please check you have named the feature types in correct naming format."
        )
    }

    # need to retrieve the output report structure from the package
    output_report <- system.file("extdata", "output_report.Rmd",
        package = "scFeatures"
    )

    file.copy(from = output_report, to = output_folder, overwrite = FALSE)

    # generate the html output
    rmarkdown::render(
        input = paste0(output_folder, "/", "output_report.Rmd"),
        output_format = "html_document",
        output_file = "output_report.html",
        output_dir = output_folder
    )
}
