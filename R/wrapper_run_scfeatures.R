


#' Wrapper function to run all feature types in scFeatures
#' 
#' @description The scFeatures function generates a variety of features from a Seurat object 
#' containing single cell RNA-sequencing data. By default, all feature types will be generated 
#' and returned in a single list containing multiple data frames.  
#' 
#' @param data input data, a Seurat object containing "sample" and "celltype" column.
#'  "x_cord" and "y_cord" is also required for constructing the features in the spatial metrics category.
#'
#' @param feature_types vector containing the name of the feature types to generate,
#' options are "proportion_raw", "proportion_logit" , "proportion_ratio",
#' "gene_mean_celltype", "gene_prop_celltype", "gene_cor_celltype",
#' "pathway_gsva" , "pathway_mean", "pathway_prop", 
#' "CCI", 
#' "gene_mean_aggregated", "gene_prop_aggregated", 'gene_cor_aggregated',
#' "L_stats" , "celltype_interaction" , "morans_I", "nn_correlation".
#'  If no value is provided, all the above feature types will be generated.
#'
#' @param type input data type, either "scrna" (stands for single-cell RNA-sequencing data),
#' "spatial_p" (stands for spatial proteomics data), or "spatial_t" (stands for single cell spatial data )
#' @param ncores number of cores , default to 1
#' @param species either "Homo sapiens" or "Mus musculus". Defaults to "Homo sapiens" if no value provided
#' @param celltype_genes the genes of interest for celltype specific gene expression feature category
#' If no value is provided, the top variable genes will be used
#' @param aggregated_genes the genes of interest for overall aggregated gene expression feature category
#' If no value is provided, the top variable genes will be used
#' @param geneset the geneset of interest for celltype specific pathway feature category
#' If no value is provided, the 50 hallmark pathways will be used
#' @param sample the sample identifier if using a SingleCellExperiment
#' @param celltype the celltype identifier if using a SingleCellExperiment
#' @param assay the assay identifier if using a SingleCellExperiment
#' @param spatialCoords the spatialCoords identifiers if using a SingleCellExperiment
#'
#' @return a list of dataframes containing the generated feature matrix in the form of sample x features
#'
#' @examples
#'
#' \dontrun{
#'  data <- readRDS(system.file("extdata", "example_scrnaseq.rds", package = "scFeatures"))
#'  scfeatures_result <- scFeatures(data, type = "scrna")
#' }
#' 
#' @export
scFeatures <- function(data, feature_types = NULL, type = "scrna", ncores = 1,
    species = "Homo sapiens", celltype_genes = NULL, aggregated_genes = NULL, geneset = NULL,
    sample = "sample", celltype = "celltype", assay = "logcounts", spatialCoords = NULL) {
    data <- makeSeurat(data, sample, celltype, assay, spatialCoords)

    return_list <- list()

    # if null, generate everything
    if (is.null(feature_types)) {
        feature_types <- c(
            "proportion_raw", "proportion_logit", "proportion_ratio",
            "gene_mean_celltype", "gene_prop_celltype", "gene_cor_celltype",
            "pathway_gsva", "pathway_mean", "pathway_prop",
            "CCI", 
            "gene_mean_aggregated", "gene_prop_aggregated", "gene_cor_aggregated",
            "L_stats", "celltype_interaction", "morans_I", "nn_correlation"
        )
    }


    for (thisfeature in feature_types) {
        try({
            if (thisfeature == "proportion_raw") {
                message("generating proportion raw features")
                return_list[["proportion_raw"]] <- run_proportion_raw(data = data, type = type, ncores = ncores)
            }

            if (thisfeature == "proportion_logit") {
                message("generating proportion logit features")
                return_list[["proportion_logit"]] <- run_proportion_logit(data = data, type = type, ncores = ncores)
            }

            if (thisfeature == "proportion_ratio") {
                message("generating proportion ratio features")
                return_list[["proportion_ratio"]] <- run_proportion_ratio(data = data, type = type, ncores = ncores)
            }


            data_remove_mito <- remove_mito(data)

            if (thisfeature == "gene_mean_celltype") {
                message("generating gene mean celltype features")
                return_list[["gene_mean_celltype"]] <- run_gene_mean_celltype(
                    data = data_remove_mito, type = type, ncores = ncores,
                    genes = celltype_genes
                )
            }


            if (thisfeature == "gene_prop_celltype") {
                message("generating gene prop celltype features")
                return_list[["gene_prop_celltype"]] <- run_gene_prop_celltype(data_remove_mito,
                    type = type, ncores = ncores,
                    genes = celltype_genes
                )
            }

            if (thisfeature == "gene_cor_celltype") {
                message("generating gene cor celltype features")
                return_list[["gene_cor_celltype"]] <- run_gene_cor_celltype(data_remove_mito,
                    type = type, ncores = ncores,
                    genes = celltype_genes
                )
            }

            if (thisfeature == "pathway_gsva") {
                message("generating pathway GSVA features")
                return_list[["pathway_gsva"]] <- run_pathway_gsva(data,
                    type = type, ncores = ncores,
                    species = species, geneset = geneset
                )
            }

            if (thisfeature == "pathway_mean") {
                message("generating pathway mean features")
                return_list[["pathway_mean"]] <- run_pathway_mean(data,
                    type = type, ncores = ncores,
                    species = species, geneset = geneset
                )
            }

            if (thisfeature == "pathway_prop") {
                message("generating pathway prop features")
                return_list[["pathway_prop"]] <- run_pathway_prop(data,
                    type = type, ncores = ncores,
                    species = species, geneset = geneset
                )
            }

            
            if (thisfeature == "CCI"){
                print("generating CCI features")
                 return_list[["CCI"]] <- run_CCI(data , 
                      type = type, ncores = ncores 
                )
            }
            
            
            if (thisfeature == "gene_mean_aggregated") {
                message("generating gene mean aggregated features")
                return_list[["gene_mean_bulk"]] <- run_gene_mean(data,
                    type = type, ncores = ncores,
                    genes = aggregated_genes
                )
            }

            if (thisfeature == "gene_prop_aggregated") {
                message("generating gene prop aggregated features")
                return_list[["gene_prop_bulk"]] <- run_gene_prop(data,
                    type = type, ncores = ncores,
                    genes = aggregated_genes
                )
            }

            if (thisfeature == "gene_cor_aggregated") {
                message("generating gene cor aggregated features")
                return_list[["gene_cor_bulk"]] <- run_gene_cor(data,
                    type = type, ncores = ncores,
                    genes = aggregated_genes
                )
            }

            if (thisfeature == "L_stats") {
                message("generating L function features")
                return_list[["L_stats"]] <- run_L_function(data, type = type, ncores = ncores)
            }

            if (thisfeature == "celltype_interaction") {
                message("generating cell type interaction features")
                return_list[["celltype_interaction"]] <- run_celltype_interaction(data, type = type, ncores = ncores)
            }

            if (thisfeature == "morans_I") {
                message("generating Moran's I features")
                return_list[["morans_I"]] <- run_Morans_I(data, type = type, ncores = ncores)
            }

            if (thisfeature == "nn_correlation") {
                message("generating nearest neighbour correlation features")
                return_list[["nn_correlation"]] <- run_nn_correlation(data, type = type, ncores = ncores)
            }
        })
    }

    lapply(return_list, as.data.frame)
    # return (return_list)
}



#' Format data into Seurat object structured for scFeatures functions
#' 
#' @description 
#' This function is used to convert a SingleCellExperiment, SpatialExperiment or a
#' Seurat object into Seurat object containing all required fields and structured 
#' for scFeatures functions. 
#'
#' @param data input data, either a SingleCellExperiment or SpatialExperiment
#'  object.
#'  The object needs to contain a column named "sample" and a column named
#'  "celltype".
#'  Alternatively, user can provide the name of the column containing sample
#'  and celltype into the `sample` and `celltype` argument.
#'  When passing as SingleCellExperiment or SpatialExperiment, by default we
#'  use the assay stored in "logcount".
#'  Alternatively, user can specify the assay to use in the `assay` argument.
#'  If users want to construct features from the spatial category, by default
#'  we need columns called "x_cord" and "y_cord".
#'  Alternatively, please specify the relevant column in the `spatialCoords`
#'  argument.
#'  For spot-based spatial transcriptomics, we also requires a matrix
#'  containing cell type prediction probability of each spot, in the format of
#'  celltype x spot
#'
#' @param sample a vector providing sample identifier for each cell. If not
#' provided, we assume the data contain a metadata column "sample" for running
#' scFeatures. 
#' @param celltype a vector providing celltype identifier. If not
#' provided, we assume the data contain a metadata column "celltype" for 
#' running scFeatures. 
#' @param assay the assay identifier if using a SingleCellExperiment or
#'  SpatialExperiment object.
#' @param spatialCoords the spatialCoords identifiers provided in a list of 
#' two vectors, if users want to construct features from the spatial category.
#' If not provided, we assume the data contain the metadata columns "x_cord" and 
#' "y_cord" for constructing spatial features. 
#' @param spotProbability a matrix in the format of celltype x spot, where each
#' entry is the prediction probability of that cell type for each spot. This is
#' needed by spatial transcriptomics data. 
#'
#' @return A `Seurat` dataset containing required metadata for running scFeatures. 
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom SpatialExperiment spatialCoords
#'
#' @examples
#'
#' \dontrun{
#'  data <- readRDS(system.file("extdata", "example_scrnaseq.rds", package = "scFeatures"))
#'  coordinate <- list(x = rep(1, ncol(data)), y = rep(1, ncol(data)))
#'  data <- makeSeurat(data, spatialCoords = coordinate)
#' }
#' @export
makeSeurat <- function(data,
    sample,
    celltype,
    assay,
    spatialCoords,
    spotProbability = NULL) {
    if (is(data, "Seurat")) {
        if (!is.null(celltype)){
            data$celltype <- celltype
        }
        if (!is.null(sample)){
            data$sample <- sample
        }

        if (!is.null(spatialCoords)) {
            data$x_cord <-  spatialCoords[1]
            data$y_cord <-  spatialCoords[2]
        }
        if (!is.null(spotProbability)) {
            data[["predictions"]] <- data[["RNA"]]
            data@assays$predictions@data <- spotProbability
        }
        return(data)
    }

    if (is(data, "SingleCellExperiment")) {
        df <- data

        if(!is.null(celltype)){
            df$celltype <- celltype
        }else{
            df$celltype <- SummarizedExperiment::colData(df)[, "celltype"]
        }

        if(!is.null(sample)){
            df$sample <- sample
        }else{
            df$sample <- SummarizedExperiment::colData(df)[, "sample"]
        }
        if (!is.null(spatialCoords)) {
            df$x_cord <- spatialCoords[1]
            df$y_cord <- spatialCoords[2]
        }else{
            df$x_cord <- SummarizedExperiment::colData(df)[, "x_cord"]
            df$y_cord <- SummarizedExperiment::colData(df)[, "y_cord"]
        }

        
        if (!is.null(assay)){
            data <- Seurat::as.Seurat(df, data = assay)
        }else{
            data <- Seurat::as.Seurat(df)
        }
        
        data@assays$RNA <- data@assays$originalexp

        if (!is.null(spotProbability)) {
            data@assays$predictions <- spotProbability
        }
        return(data)
    }


    if (is(data, "SpatialExperiment")) {
        df <- data

        if(!is.null(celltype)){
            df$celltype <- celltype
        }else{
            df$celltype <- SummarizedExperiment::colData(df)[, "celltype"]
        }

        if(!is.null(sample)){
            df$sample <- sample
        }else{
            df$sample <- SummarizedExperiment::colData(df)[, "sample"]
        }

        if(!is.null(spatialCoords)){
            df$x_cord <- spatialCoords[1]
            df$y_cord <- spatialCoords[2]
        }else{
            df$x_cord <- SpatialExperiment::spatialCoords(df)[, 1]
            df$y_cord <- SpatialExperiment::spatialCoords(df)[, 2]
        }
        
        if (!is.null(assay)){
            data <- Seurat::as.Seurat(df, data = assay)
        }else{  
            data <- Seurat::as.Seurat(df)
        }
        
        data@assays$RNA <- data@assays$originalexp

        if (!is.null(spotProbability)) {
            data@assays$predictions <- spotProbability
        }

        return(data)
    }

    data
}
