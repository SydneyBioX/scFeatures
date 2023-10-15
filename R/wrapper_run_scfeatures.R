


#' Wrapper function to run all feature types in scFeatures
#' 
#' @description The scFeatures function generates a variety of features from a Seurat object 
#' containing single cell RNA-sequencing data. By default, all feature types will be generated 
#' and returned in a single list containing multiple data frames.  
#' 
#' @param data input data, a matrix of genes by cells 
#' @param sample a vector of sample information
#' @param celltype a vector of cell type information
#' @param spatialCoords a list of two vectors containing the x and y coordinates of each cell   
#' @param spotProbability a matrix of spot probability, each row represents a celltype and each column represents a spot  
#' 
#' @param feature_types vector containing the name of the feature types to generate,
#' options are "proportion_raw", "proportion_logit" , "proportion_ratio",
#' "gene_mean_celltype", "gene_prop_celltype", "gene_cor_celltype",
#' "pathway_gsva" , "pathway_mean", "pathway_prop", 
#' "CCI", 
#' "gene_mean_aggregated", "gene_prop_aggregated", 'gene_cor_aggregated',
#' "L_stats" , "celltype_interaction" , "morans_I", "nn_correlation".
#'  If no value is provided, all the above feature types will be generated.
#' @param type input data type, either "scrna" (stands for single-cell RNA-sequencing data),
#' "spatial_p" (stands for spatial proteomics data), or "spatial_t" (stands for single cell spatial data )
#' @param ncores number of cores , default to 1
#' 
#' @param species either "Homo sapiens" or "Mus musculus". Defaults to "Homo sapiens" if no value provided
#' @param celltype_genes the genes of interest for celltype specific gene expression feature category
#' If no value is provided, the top variable genes will be used
#' @param aggregated_genes the genes of interest for overall aggregated gene expression feature category
#' If no value is provided, the top variable genes will be used
#' @param geneset the geneset of interest for celltype specific pathway feature category
#' If no value is provided, the 50 hallmark pathways will be used

#'
#' @return a list of dataframes containing the generated feature matrix in the form of sample x features
#'
#' @examples
#' data("example_scrnaseq" , package = "scFeatures") 
#' data <- example_scrnaseq
#' celltype <- data$celltype
#' sample <- data$sample
#' data <- data@assays$RNA@data
#' scfeatures_result <- scFeatures(data, celltype = celltype, sample = sample, type = "scrna", feature_types = "proportion_raw")
#' 
#' @export
scFeatures <- function(data = NULL, sample = NULL ,  celltype = NULL, 
                       spatialCoords = NULL,  spotProbability = NULL , 
                       feature_types = NULL, type = "scrna", ncores = 1,
    species = "Homo sapiens", celltype_genes = NULL, aggregated_genes = NULL, geneset = NULL ) {
    
    alldata <- formatData( data, sample, celltype,  spatialCoords , spotProbability)
  
    
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
                return_list[["proportion_raw"]] <- run_proportion_raw(data = alldata, type = type, ncores = ncores)
            }

            if (thisfeature == "proportion_logit") {
                message("generating proportion logit features")
                return_list[["proportion_logit"]] <- run_proportion_logit(data = alldata, type = type, ncores = ncores)
            }

            if (thisfeature == "proportion_ratio") {
                message("generating proportion ratio features")
                return_list[["proportion_ratio"]] <- run_proportion_ratio(data = alldata, type = type, ncores = ncores)
            }


            data_remove_mito <- remove_mito_ribo(alldata)

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
                return_list[["pathway_gsva"]] <- run_pathway_gsva(alldata,
                    type = type, ncores = ncores,
                    species = species, geneset = geneset
                )
            }

            if (thisfeature == "pathway_mean") {
                message("generating pathway mean features")
                return_list[["pathway_mean"]] <- run_pathway_mean(alldata,
                    type = type, ncores = ncores,
                    species = species, geneset = geneset
                )
            }

            if (thisfeature == "pathway_prop") {
                message("generating pathway prop features")
                return_list[["pathway_prop"]] <- run_pathway_prop(alldata,
                    type = type, ncores = ncores,
                    species = species, geneset = geneset
                )
            }

            
            if (thisfeature == "CCI"){
                message("generating CCI features")
                 return_list[["CCI"]] <- run_CCI(alldata , 
                      type = type, ncores = ncores 
                )
            }
            
            
            if (thisfeature == "gene_mean_aggregated") {
                message("generating gene mean aggregated features")
                return_list[["gene_mean_bulk"]] <- run_gene_mean(alldata,
                    type = type, ncores = ncores,
                    genes = aggregated_genes
                )
            }

            if (thisfeature == "gene_prop_aggregated") {
                message("generating gene prop aggregated features")
                return_list[["gene_prop_bulk"]] <- run_gene_prop(alldata,
                    type = type, ncores = ncores,
                    genes = aggregated_genes
                )
            }

            if (thisfeature == "gene_cor_aggregated") {
                message("generating gene cor aggregated features")
                return_list[["gene_cor_bulk"]] <- run_gene_cor(alldata,
                    type = type, ncores = ncores,
                    genes = aggregated_genes
                )
            }

            if (thisfeature == "L_stats") {
                message("generating L function features")
                return_list[["L_stats"]] <- run_L_function(alldata, type = type, ncores = ncores)
            }

            if (thisfeature == "celltype_interaction") {
                message("generating cell type interaction features")
                return_list[["celltype_interaction"]] <- run_celltype_interaction(alldata, type = type, ncores = ncores)
            }

            if (thisfeature == "morans_I") {
                message("generating Moran's I features")
                return_list[["morans_I"]] <- run_Morans_I(alldata, type = type, ncores = ncores)
            }

            if (thisfeature == "nn_correlation") {
                message("generating nearest neighbour correlation features")
                return_list[["nn_correlation"]] <- run_nn_correlation(alldata, type = type, ncores = ncores)
            }
        })
    }

    lapply(return_list, as.data.frame)
    # return (return_list)
}


formatData <- function(data = NULL,
                       sample = NULL ,
                       celltype = NULL,
                       spatialCoords = NULL,
                       spotProbability = NULL) {
  
  alldata <- list()
  
  if (is.null(data)){
    print("Please provide data as gene x cell matrix")
  }else{
    alldata$data <- data
  }
  if (is.null(sample)){
    print("Please provide the sample name of each cell")
  }else{
    alldata$sample <- sample
  }
  
  
  if (is.null(celltype)){
    
  }else{
    alldata$celltype <- celltype
  }
  
 
  if (is.null(spatialCoords)) {
    
  }else{
    alldata$x_cord <- spatialCoords[[1]]
    alldata$y_cord <- spatialCoords[[2]]
  }
 
  if (is.null( spotProbability)){
    
  }else{
    alldata$predictions <- spotProbability
  }
  
 return(alldata) 
}
  

 


