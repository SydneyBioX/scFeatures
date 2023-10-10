


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
#' data("example_scrnaseq" , package = "scFeatures") 
#' data <- example_scrnaseq
#' scfeatures_result <- scFeatures(data, type = "scrna", feature_types = "proportion_raw")
#' 
#' @export
scFeatures <- function(data = NULL, sample = NULL ,  celltype = NULL, 
                       spatialCoords = NULL,  spotProbability = NULL , 
                       feature_types = NULL, type = "scrna", ncores = 1,
    species = "Homo sapiens", celltype_genes = NULL, aggregated_genes = NULL, geneset = NULL ) {
  
    # data <- makeSeurat(data, sample, celltype, assay, spatialCoords)
    
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
#'  data("example_scrnaseq" , package = "scFeatures")
#'  data <- example_scrnaseq
#'  coordinate <- list(x = rep(1, ncol(data)), y = rep(1, ncol(data)))
#'  data <- makeSeurat(data, spatialCoords = coordinate)
#' 
#' @export
makeSeurat <- function(data,
                       sample = NULL ,
                       celltype = NULL,
                       select_assay = NULL,
                       spatialCoords = NULL,
                       spotProbability = NULL) {
  
  
  if (!is.null(spotProbability)){
    spotProbability  <- t(spotProbability )
    predictions <- data.frame( 
      as.matrix(spotProbability ),
      row.names = colnames(data),  stringsAsFactors = FALSE )
    
    predictions <- Seurat::CreateAssayObject( data = t(x = as.matrix(x = predictions)), 
                                              check.matrix = FALSE )
  }
  
  
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
      data[["predictions"]] <-  predictions
    }
    return(data)
  }
  
  if (is(data, "SingleCellExperiment")) {
    df <- data
    
    if(!is.null(celltype)){
      df$celltype <- celltype
    } 
    
    if(!is.null(sample)){
      df$sample <- sample
    }
    
    if (!is.null(spatialCoords)) {
      df$x_cord <- spatialCoords[[1]]
      df$y_cord <- spatialCoords[[2]]
    } 
    
    
    
    if (!is.null(select_assay)){
      data_seurat <-  CreateSeuratObject( counts = assay(df,  select_assay))
    }else{
      data_seurat <- CreateSeuratObject( counts = assay(df,  "logcounts"))
    }
    
    
    data_seurat@meta.data <-  data.frame( colData(df) )

    if (!is.null(spotProbability)) {
      data_seurat[["predictions"]] <-  predictions
    }
    return(data_seurat)
  }

 
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
  

 


