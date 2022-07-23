


#' wrapper function to run all feature types in scFeatures
#'
#' @param data input data, a Seurat object containing `celltype` and `sample` label
#' @param feature_types vector containing the name of the feature types to generate, 
#' options are "proportion_raw", "proportion_logit" , "proportion_ratio", 
#' "gene_mean_celltype", "gene_prop_celltype", "gene_cor_celltype",  
#' "pathway_gsva" , "pathway_mean", "pathway_prop" , "CCI", 
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
#' data <- readRDS(system.file("extdata", "example_scrnaseq.rds", package = "scFeatures"))
#' scfeatures_result <- scFeatures(data, type = "scrna" )
#' 
#' @import parallel
#' @export
scFeatures <- function( data , feature_types = NULL , type =  "scrna", ncores  = 1 ,
                        species  = "Homo sapiens", celltype_genes = NULL, aggregated_genes = NULL , geneset = NULL,
                      sample = "sample", celltype = "celltype", celltype = "logcounts", spatialCoords = c("x_cord", "y_cord")){
  
  data <- makeSeurat(data, sample, celltype, assay, spatialCoords)
  
  return_list <- list()
  
  # if null, generate everything 
  if ( is.null(feature_types)){
    feature_types  <-   c("proportion_raw", "proportion_logit" , "proportion_ratio", 
                    "gene_mean_celltype", "gene_prop_celltype", "gene_cor_celltype", 
                    "pathway_gsva" , "pathway_mean", "pathway_prop" , 
                    "CCI", 
                    "gene_mean_aggregated", "gene_prop_aggregated", 'gene_cor_aggregated', 
                    "L_stats" , "celltype_interaction" , "morans_I", "nn_correlation")
  }
 
  
  for (thisfeature in feature_types ){
    
    try({
      
       if (thisfeature == "proportion_raw"){
         print("generating proportion raw features")
         return_list[["feature_proportion_raw"]] <- run_proportion_raw( data = data , type = type, ncores = ncores )
       } 
      
       if (thisfeature == "proportion_logit"){
         print("generating proportion logit features")
         return_list[["feature_proportion_logit"]] <-  run_proportion_logit( data = data , type = type, ncores = ncores )
       }
      
      if (thisfeature ==  "proportion_ratio" ){ 
        print("generating proportion ratio features")
        return_list[["feature_proportion_ratio"]] <- run_proportion_ratio( data = data , type = type, ncores = ncores )
      }
        
      
      data_remove_mito <- remove_mito(data)
      
      if (thisfeature == "gene_mean_celltype"){
        print("generating gene mean celltype features")
        return_list[["feature_gene_mean_celltype"]] <- run_gene_mean_celltype( data = data_remove_mito , type = type, ncores = ncores, 
                                                                               genes = celltype_genes ) 
      }
      
      
      if (thisfeature == "gene_prop_celltype"){
        print("generating gene prop celltype features")
        return_list[["feature_gene_prop_celltype"]] <- run_gene_prop_celltype(data_remove_mito , type = type, ncores = ncores,
                                                                              genes = celltype_genes )
      }
       
      if (thisfeature ==  "gene_cor_celltype"){
        print("generating gene cor celltype features")
        return_list[["feature_gene_cor_celltype"]] <- run_gene_cor_celltype(data_remove_mito , type = type, ncores = ncores,
                                                                             genes = celltype_genes )
      }
      
      if (thisfeature == "pathway_gsva"  ){
        print("generating pathway GSVA features")
        return_list[["feature_pathway_gsva"]] <- run_pathway_gsva(data, type = type, ncores = ncores ,
                                                                  species = species  , geneset = geneset ) 
      }
        
      if (thisfeature ==  "pathway_mean" ){
        print("generating pathway mean features")
        return_list[["feature_pathway_mean"]] <- run_pathway_mean(data ,  type = type, ncores = ncores , 
                                                                  species = species , geneset = geneset )
      }
      
      if (thisfeature ==  "pathway_prop" ){
        print("generating pathway prop features")
        return_list[["feature_pathway_prop"]] <- run_pathway_prop(data,  type = type, ncores = ncores ,
                                                                  species = species  , geneset = geneset )
      }
        
     
      if (thisfeature == "CCI"){
        print("generating CCI features")
        return_list[["feature_CCI"]] <- run_CCI(data , type = type, ncores = ncores ,   species = species   )
      }
      
      if (thisfeature ==  "gene_mean_aggregated" ){
        print("generating gene mean aggregated features")
        return_list[["feature_gene_mean_bulk"]] <- run_gene_mean(data,  type = type, ncores = ncores,
                                                                 genes = aggregated_genes )
      }

      if (thisfeature == "gene_prop_aggregated" ){
        print("generating gene prop aggregated features")
        return_list[["feature_gene_prop_bulk"]] <- run_gene_prop(data, type = type, ncores = ncores,
                                                                 genes= aggregated_genes )
      }
      
      if (thisfeature == "gene_cor_aggregated" ){
        print("generating gene cor aggregated features")
        return_list[["feature_gene_cor_bulk"]] <- run_gene_cor(data , type = type, ncores = ncores,
                                                               genes = aggregated_genes )
      }
      
      if (thisfeature ==  "L_stats" ){
        print("generating L function features")
        return_list[["feature_L_stats"]] <- run_L_function( data, type = type, ncores = ncores  )
      }
      
      if (thisfeature ==  "celltype_interaction"){
        print("generating cell type interaction features")
        return_list[["feature_celltype_interaction"]] <- run_celltype_interaction( data,  type = type, ncores = ncores  )
      }
      
      if (thisfeature == "morans_I"){
        print("generating Moran's I features")
        return_list[["feature_morans_I"]] <- run_Morans_I( data, type = type, ncores = ncores  )
      }
    
      if (thisfeature ==  "nn_correlation"){
        print("generating nearest neighbour correlation features")
        feature_nn_correlation <- run_nn_correlation( data,  type = type, ncores = ncores   )
      }
    
   
    })
  }
  
  
  return (return_list)
}



#' @importFrom SummarizedExperiment colData
#' @importFrom SpatialExperiment spatialCoords
makeSeurat <- function(data, sample, celltype, assay, spatialCoords){

  if(is(data, "SingleCellExperiment")){
    df <- data
    df$celltype <- SummarizedExperiment::colData(df)[celltype]
    df$sample <- SummarizedExperiment::colData(df)[sample]
    df$x_cord <- SummarizedExperiment::colData(df)[spatialCoords[1]] 
    df$y_cord <- SummarizedExperiment::colData(df)[spatialCoords[2]] 
    data <- Seurat::as.Seurat(cells, data = assay)
    data@assays$RNA <- data@assays$originalexp
    return(data)
    }
  
  if(is(data, "SpatialExperiment")){
  df <- data
  df$celltype <- SummarizedExperiment::colData(df)[celltype]
  df$sample <- SummarizedExperiment::colData(df)[sample]
  df$x_cord <- SpatialExperiment::spatialCoords(df)[,1] 
  df$y_cord <- SpatialExperiment::spatialCoords(df)[,2] 
  data <- Seurat::as.Seurat(cells, data = assay)
  data@assays$RNA <- data@assays$originalexp
  return(data)
  }
    
  data
}
    

