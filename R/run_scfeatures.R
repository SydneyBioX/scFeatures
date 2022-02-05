# this file contains the 17 functions that generate the 17 feature classes 



#' generate cell type proportion raw  
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
#' feature_proportion_raw <- run_proportion_raw(data, type = "scnrna", ncores = 1)
#' 
#' @import gtools
#' @import tidyr
#' @import parallel
#' 
#' @export
run_proportion_raw <- function( data, type = "scrna" , ncores = 1 ){
 
  check_data(data, type)
  
  if ( type %in% c( "scrna" , "spatial_p"))  {
        X <- helper_proportion_raw( data, logit = F )
  }
  
  
  if ( type == "spatial_t" ) {
       X <- helper_proportion_raw_st( data, logit = F , ncores )
  }
 
  return (X)
  
}




#' generate cell type proportion raw and cell type proportion logit transformed
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
#' feature_proportion_logit <- run_proportion_logit(data, type = "scnrna", ncores = 1)
#' 
#' @import gtools
#' @import tidyr
#' @import parallel
#' 
#' @export
run_proportion_logit <- function( data, type = "scrna" ,  ncores = 1 ){
 
  check_data(data, type)
  
  if ( type %in% c( "scrna" , "spatial_p"))  {
     X <- helper_proportion_raw(data , logit = T)
  }
  
  
  if ( type == "spatial_t" ) {
    X <- helper_proportion_raw_st(data , logit = T , ncores )
  }
  
  
  
  return (X)
  
}



#' generate cell type proportion ratio 
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
#' feature_proportion_ratio <- run_proportion_ratio(data, type = "scnrna", ncores = 1)
#' 
#' @import gtools
#' @import tidyr
#' @import parallel
#' 
#' @export
run_proportion_ratio <- function( data, type = "scrna" , ncores = 1 ){
 
  check_data(data, type)
  
  if ( type %in% c( "scrna" , "spatial_p") )  {
     X <- helper_proportion_ratio(data , ncores   ) 
  }
  
  
  if ( type == "spatial_t" ) {
    X <- helper_proportion_ratio_st(  data , ncores )
  }
  
  return (X)
  
}










#' generate cell type specific gene mean expression 
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
#' data_remove_mito <- remove_mito(data) #optional step, if mito and ribo genes are not of interest
#' feature_gene_mean_celltype <- run_gene_mean_celltype(data_remove_mito, type = "scnrna", num_top_gene = 100, ncores = 1) 
#' 
#' @import proxyC
#' @import DelayedMatrixStats
#' @import parallel
#' @import DelayedArray
#' 
#' @export
run_gene_mean_celltype  <- function( data, type = "scrna" , num_top_gene = NULL , ncores = 1 ){
 
  check_data(data, type)
  
  if ( type %in% c( "scrna" , "spatial_p") )  {
    X <- helper_gene_mean_celltype(data ,  num_top_gene , ncores)
  }
  
  if ( type == "spatial_t" ) {
    X <- helper_gene_mean_celltype_st(data)
  }

  return (X)
  
}





#' generate cell type specific gene proportion expression
#'
#' @param data input data, a Seurat object containing `celltype` and `sample` label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param num_top_gene determines how many genes to include. By default, picks around 100 most variable genes per cell type.   
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#'
#' @examples
#'
#' data <- readRDS(system.file("extdata", "example_scrnaseq.rds", package = "scFeatures"))
#' data_remove_mito <- remove_mito(data) #optional step, if mito and ribo genes are not of interest
#' feature_gene_prop_celltype <- feature_gene_prop_celltype(data_remove_mito, type = "scnrna", num_top_gene = 100,  ncores = 1) 
#' 
#' @import proxyC
#' @import DelayedMatrixStats
#' @import parallel
#' @import DelayedArray
#' 
#' @export
run_gene_prop_celltype  <- function( data, type = "scrna" , num_top_gene = NULL, ncores = 1 ){
  
  check_data(data, type)
  
  if ( type %in% c( "scrna" , "spatial_p") )  {
    X <- helper_gene_prop_celltype(data , num_top_gene , ncores)
  }
  
  
  if ( type == "spatial_t" ) {
    print("This feature class currently does not support spatial transcriptomics")
    return(NULL)
  }
  
  
  return (X)
  
}





#' generate cell type specific gene correlation 
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
#' data_remove_mito <- remove_mito(data) #optional step, if mito and ribo genes are not of interest
#' feature_gene_cor_celltype <- feature_gene_cor_celltype(data_remove_mito, type = "scnrna", num_top_gene = 100, ncores = 1) 
#' 
#' @import proxyC
#' @import DelayedMatrixStats
#' @import parallel
#' @import DelayedArray
#' 
#' @export
run_gene_cor_celltype  <- function( data, type = "scrna" , num_top_gene = NULL, ncores = 1 ){
  
  check_data(data, type)
  
  if ( type %in% c( "scrna" , "spatial_p") )  {
    X <- helper_gene_cor_celltype(data ,  num_top_gene , ncores)
  }
  
  if ( type == "spatial_t" ) {
    print("This feature class currently does not support spatial transcriptomics")
    return(NULL)
  }
  
  return (X)
  
}






#' generate pathway score using gene set enrichement analysis 
#'
#' @param data input data, a Seurat object containing `celltype` and `sample` label
#' @param geneset By default (when the `geneset` argument is not specified),  we use the 50 hallmark gene set from msigdb. 
#' The users can also provide their geneset of interest in a list format, with each list entry containing a vector of the names of genes in a gene set.
#' eg, geneset <- list("pathway_a" = c("CAPNS1", "TLCD1"), "pathway_b" = c("PEX6","DPRXP4" )) 
#' @param species whether the species is "Homo sapiens" or "Mus musculus". Default is "Homo sapiens". 
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param subsample whether to subsample, either T or F. For larger datasets (eg, over 30,000 cells), the subsample function can be used to increase speed. 
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#'
#' @examples
#'
#' data <- readRDS(system.file("extdata", "example_scrnaseq.rds", package = "scFeatures"))
#' feature_pathway_gsva <- run_pathway_gsva(data, geneset = NULL, species = "Homo sapiens",  type = "scrna" , subsample = F,  ncores = 1 )
#' 
#' @import msigdbr
#' @import ensembldb
#' @import GSVA
#' @import EnsDb.Hsapiens.v79
#' @import EnsDb.Mmusculus.v79
#' @import ensembldb
#' @import dplyr
#' @import reshape2
#' @import DelayedArray
#' @import parallel
#' 
#' @export
run_pathway_gsva <- function( data, geneset = NULL , 
                              species = "Homo sapiens" ,
                              type = "scrna" ,  subsample = T , ncores = 1  ){
  
  check_data(data, type)
  
  if ( is.null(geneset) ){
    geneset <- get_geneset(species = species)
  }
  
  if (subsample & ncol(data) > 90000 ){
     data <- data[ , sample(1:ncol(data), 90000) ]
  }
  
  if ( type == "scrna" )  {
    # if the user does not provide geneset, need to get the geneset from msigdb
    X <- helper_pathway_gsva(data, geneset = geneset , ncores = ncores )
  }
  
  if ( type ==  "spatial_p") {
    print("This feature class currently does not support spatial proteomics")
    return(NULL)
  }
  
  if ( type == "spatial_t" ) {
    print("This feature class currently does not support spatial transcriptomics")
    return(NULL)
  }
  
  return (X)
  
}






#' generate pathway score using expression level
#' 
#' @param data input data, a Seurat object containing `celltype` and `sample` label
#' @param geneset By default (when the `geneset` argument is not specified),  we use the 50 hallmark gene set from msigdb. 
#' The users can also provide their geneset of interest in a list format, with each list entry containing a vector of the names of genes in a gene set.
#' eg, geneset <- list("pathway_a" = c("CAPNS1", "TLCD1"), "pathway_b" = c("PEX6","DPRXP4" )) 
#' @param species whether the species is "Homo sapiens" or "Mus musculus". Default is "Homo sapiens". 
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#'
#' @examples
#'
#' data <- readRDS(system.file("extdata", "example_scrnaseq.rds", package = "scFeatures"))
#' feature_pathway_mean <- run_pathway_mean(data, geneset = NULL, species = "Homo sapiens",  type = "scrna" ,  ncores = 1 )
#' 
#' @import msigdbr
#' @import ensembldb
#' @import GSVA
#' @import EnsDb.Hsapiens.v79
#' @import EnsDb.Mmusculus.v79
#' @import ensembldb
#' @import dplyr
#' @import reshape2
#' @import DelayedArray
#' @import parallel
#' 
#' @export
run_pathway_mean <- function( data, geneset = NULL , 
                              species = "Homo sapiens" , type = "scrna" , ncores = 1  ){
  
  check_data(data, type)
  
  if ( is.null(geneset) ){
    geneset <- get_geneset(species = species)
  }
  
  if ( type == "scrna" )  {
    X <- helper_pathway_mean(data  , geneset = geneset , ncores = ncores )
  }
  
  if ( type == "spatial_p" )  {
    print("This feature class currently does not support spatial proteomics")
    return(NULL)
  }
  
  if ( type == "spatial_t" ) {
    X <- helper_pathway_mean_st( data  , geneset = geneset , ncores = ncores )
  }

  return (X)
  
}





#' generate pathway score using proportion of expression 
#' 
#' @param data input data, a Seurat object containing `celltype` and `sample` label
#' @param geneset By default (when the `geneset` argument is not specified),  we use the 50 hallmark gene set from msigdb. 
#' The users can also provide their geneset of interest in a list format, with each list entry containing a vector of the names of genes in a gene set.
#' eg, geneset <- list("pathway_a" = c("CAPNS1", "TLCD1"), "pathway_b" = c("PEX6","DPRXP4" ))  
#' @param species whether the species is "Homo sapiens" or "Mus musculus". Default is "Homo sapiens". 
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#'
#' @examples
#'
#' data <- readRDS(system.file("extdata", "example_scrnaseq.rds", package = "scFeatures"))
#' feature_pathway_prop <- run_pathway_prop(data, geneset = NULL, species = "Homo sapiens",  type = "scrna" ,  ncores = 1 )
#' 
#' @import msigdbr
#' @import ensembldb
#' @import GSVA
#' @import EnsDb.Hsapiens.v79
#' @import EnsDb.Mmusculus.v79
#' @import ensembldb
#' @import dplyr
#' @import reshape2
#' @import DelayedArray
#' @import parallel
#' 
#' @export
run_pathway_prop <- function( data, geneset = NULL , 
                              species = "Homo sapiens" ,  type = "scrna" , ncores = 1 ){
 
  check_data(data, type)
  
  if ( is.null(geneset) ){
    geneset <- get_geneset(species = species)
  }
  
  if ( type == "scrna" )  {
    X <- helper_pathway_prop(data  , geneset = geneset ,  ncores = ncores )
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








#' generate cell cell communication score
#'
#' @param data input data, a Seurat object containing `celltype` and `sample` label
#' @param species whether the species is "Homo sapiens" or "Mus musculus". Default is "Homo sapiens". 
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @examples
#'
#' data <- readRDS(system.file("extdata", "example_scrnaseq.rds", package = "scFeatures"))
#' feature_CCI <- run_CCI(data, species = "Homo sapiens",  type = "scrna" ,  ncores = 1 )
#' 
#' @import dplyr
#' @import reshape2
#' @import DelayedArray
#' @import parallel
#' @import CellChat
#' @import plyr
#' 
#' @export
run_CCI <- function( data,  species = "Homo sapiens" , type = "scrna" , ncores = 1  ){
  
  check_data(data, type)
  
  if ( type == "scrna" )  {
    X <- helper_CCI(data, species = species, ncores =  ncores )
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








#' generate overall aggregated mean expression 
#'
#' @param data input data, a Seurat object containing `celltype` and `sample` label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param num_top_gene determines how many genes to include. By default, picks 1500 most variable genes.
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @examples
#'
#' data <- readRDS(system.file("extdata", "example_scrnaseq.rds", package = "scFeatures"))
#' feature_gene_mean <- run_gene_mean(data, type = "scrna" , num_top_gene= 1500, ncores = 1 )
#' 
#' @import proxyC
#' @import DelayedMatrixStats
#' @import parallel
#' 
#' @export
run_gene_mean  <- function( data, type = "scrna" , num_top_gene= NULL, ncores = 1 ){
 
  check_data(data, type)
  
  if ( type %in% c( "scrna" , "spatial_p" ,  "spatial_t") )  {
    X <- helper_gene_mean(data,  num_top_gene , ncores )
  }
 
  return (X)
  
}





#' generate overall aggregated gene proportion expression
#'
#' @param data input data, a Seurat object containing `celltype` and `sample` label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param num_top_gene determines how many genes to include. By default, picks 1500 most variable genes.
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @examples
#'
#' data <- readRDS(system.file("extdata", "example_scrnaseq.rds", package = "scFeatures"))
#' feature_gene_prop <- run_gene_prop(data, type = "scrna" , num_top_gene= 1500, ncores = 1 )
#' 
#' @import proxyC
#' @import DelayedMatrixStats
#' @import parallel
#' @import DelayedArray
#' 
#' @export
run_gene_prop <- function( data, type = "scrna" , num_top_gene  = NULL, ncores = 1 ){
  
  check_data(data, type)
  
  if ( type %in% c( "scrna" , "spatial_p" , "spatial_t") )  {
    X <- helper_gene_prop(data ,  num_top_gene , ncores)
  }

  return (X)
  
}





#' generate overall aggregated gene correlation 
#'
#' @param data input data, a Seurat object containing `celltype` and `sample` label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param num_top_gene determines how many genes to include. By default, picks 1500 most variable genes.
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @examples
#'
#' data <- readRDS(system.file("extdata", "example_scrnaseq.rds", package = "scFeatures"))
#' feature_gene_cor <- run_gene_cor(data, type = "scrna" , num_top_gene= 1500, ncores = 1 )
#' 
#' @import proxyC
#' @import DelayedMatrixStats
#' @import parallel
#' @import DelayedArray
#' 
#' @export
run_gene_cor  <- function( data, type = "scrna" , num_top_gene = NULL, ncores = 1 ){
  
  check_data(data, type)
  
  if ( type %in% c( "scrna" , "spatial_p"  , "spatial_t" ) )  {
    X <- helper_gene_cor(data  , num_top_gene , ncores)
  }
 
  return (X)
  
}



 
#' generate L stats
#'
#' @param data input data, a Seurat object containing `celltype` and `sample` label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @examples
#'
#' data <- readRDS(system.file("extdata", "example_spatial_proteomics.rds", package = "scFeatures"))
#' feature_L_function <- run_L_function(data, type = "spatial_p", ncores = 1 )
#' 
#' @import spatstat.geom
#' @import parallel
#' 
#' @export
run_L_function <- function( data, type = "spatial_p" , ncores = 1 ){
  
  check_data(data, type)
  
  if ( type == "scrna" )  {
    print("This feature class currently does not support scRNA-seq")
    return(NULL)
  }
  
  if ( type == "spatial_p" )  {
    X <- helper_L_stat_sp(data)
  }
  
  if ( type == "spatial_t" ) {
    X <- helper_L_stat_st(data)
  }
  
  return (X)
  
}





#' generate cell type interaction 
#'
#' @param data input data, a Seurat object containing `celltype` and `sample` label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @examples
#'
#' data <- readRDS(system.file("extdata", "example_spatial_proteomics.rds", package = "scFeatures"))
#' feature_celltype_interaction <- run_celltype_interaction(data, type = "spatial_p", ncores = 1 )
#' 
#' @import spatstat.geom
#' @import parallel
#' 
#' @export
run_celltype_interaction <- function( data, type = "spatial_p" , ncores = 1  ){

  check_data(data, type)
  
  if ( type == "scrna" )  {
    print("This feature class currently does not support scRNA-seq")
    return(NULL)
  }
  
  if ( type == "spatial_p" )  {
    X <- helper_celltype_interaction_sp(data)
  }
  
  if ( type == "spatial_t" ) {
    X <- helper_celltype_interaction_st(data)
  }
  
  return (X)
  
}








#' generate Moran's I
#'
#' @param data input data, a Seurat object containing `celltype` and `sample` label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @examples
#'
#' data <- readRDS(system.file("extdata", "example_spatial_proteomics.rds", package = "scFeatures"))
#' feature_Morans_I <- run_Morans_I(data, type = "spatial_p", ncores = 1 )
#' 
#' @import spatstat.geom
#' @import parallel
#' @import ape
#' 
#' @export
run_Morans_I <- function( data, type = "spatial_p" , ncores = 1  ){
  
  check_data(data, type)
  
  if ( type == "scrna" )  {
    print("This feature class currently does not support scRNA-seq")
    return(NULL)
  }
  
  if ( type %in% c( "spatial_p" , "spatial_t" ))   {
    X <- helper_moran( data , num_top_gene =  NULL,   ncores = 1)
  }
  
  return (X)
  
}









#' generate nearest neighbour correlation
#'
#' @param data input data, a Seurat object containing `celltype` and `sample` label
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @examples
#'
#' data <- readRDS(system.file("extdata", "example_spatial_proteomics.rds", package = "scFeatures"))
#' feature_nn_correlation <- run_nn_correlation(data, type = "spatial_p", ncores = 1 )
#' 
#' @import spatstat.geom
#' @import spatstat.core
#' @import parallel
#' 
#' @export
run_nn_correlation <- function( data, type = "spatial_p", ncores = 1){

  check_data(data, type)
  
  if ( type == "scrna" )  {
    print("This feature class currently does not support scRNA-seq")
    return(NULL)
  }
  
  if ( type %in% c( "spatial_p" , "spatial_t" ))   {
    X <- helper_nncorr_protein(data, num_top_gene  =  NULL , ncores = ncores  )
  }

  return (X)
  
}








