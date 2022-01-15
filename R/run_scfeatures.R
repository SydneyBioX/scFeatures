# this file contains the 17 functions that generate the 17 feature classes 


# scrna, spatial_p , spatial_t 





#' generate cell type proportion raw  
#'
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @import gtools
#' @import tidyr
#' @import parallel
#' 
#' @export
run_proportion_raw <- function( data, type = "scrna" , ncores = 8 ){
   
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
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @import gtools
#' @import tidyr
#' @import parallel
#' 
#' @export
run_proportion_logit <- function( data, type = "scrna" ,  ncores = 8  ){
  
  
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
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @import gtools
#' @import tidyr
#' @import parallel
#' 
#' @export
run_proportion_ratio <- function( data, type = "scrna" , ncores = 8  ){
  
  
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
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @import proxyC
#' @import DelayedMatrixStats
#' @import parallel
#' @import DelayedArray
#' 
#' @export
run_gene_mean_celltype  <- function( data, type = "scrna" , num_top_gene = NULL , ncores = 8 ){
  
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
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @import proxyC
#' @import DelayedMatrixStats
#' @import parallel
#' @import DelayedArray
#' 
#' @export
run_gene_prop_celltype  <- function( data, type = "scrna" , num_top_gene = NULL, ncores = 8 ){
  
  
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
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @import proxyC
#' @import DelayedMatrixStats
#' @import parallel
#' @import DelayedArray
#' 
#' @export
run_gene_cor_celltype  <- function( data, type = "scrna" , num_top_gene = NULL, ncores = 8 ){
  
  
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
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
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
                              type = "scrna" ,  subsample = T , ncores = 8  ){
  
  
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
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
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
                              species = "Homo sapiens" , type = "scrna" , ncores = 8  ){
  
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
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
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
                              species = "Homo sapiens" ,  type = "scrna" , ncores = 8 ){
  
  
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
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' 
#' @import dplyr
#' @import reshape2
#' @import DelayedArray
#' @import parallel
#' @import CellChat
#' @import plyr
#' 
#' @export
run_CCI <- function( data,  species = "Homo sapiens" , type = "scrna" , ncores = 8  ){
  
  
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








#' generate bulk gene mean expression 
#'
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @import proxyC
#' @import DelayedMatrixStats
#' @import parallel
#' 
#' @export
run_gene_mean  <- function( data, type = "scrna" , num_top_gene= NULL, ncores = 8 ){
  
  if ( type %in% c( "scrna" , "spatial_p" ,  "spatial_t") )  {
    X <- helper_gene_mean(data,  num_top_gene , ncores )
  }
 
  return (X)
  
}





#' generate bulk gene proportion expression
#'
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @import proxyC
#' @import DelayedMatrixStats
#' @import parallel
#' @import DelayedArray
#' 
#' @export
run_gene_prop <- function( data, type = "scrna" , num_top_gene  = NULL, ncores = 8 ){
  
  
  if ( type %in% c( "scrna" , "spatial_p" , "spatial_t") )  {
    X <- helper_gene_prop(data ,  num_top_gene , ncores)
  }
  
 
  return (X)
  
}





#' generate bulk gene correlation 
#'
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @import proxyC
#' @import DelayedMatrixStats
#' @import parallel
#' @import DelayedArray
#' 
#' @export
run_gene_cor  <- function( data, type = "scrna" , num_top_gene = NULL, ncores = 8 ){
  
  
  if ( type %in% c( "scrna" , "spatial_p"  , "spatial_t" ) )  {
    X <- helper_gene_cor(data  , num_top_gene , ncores)
  }
 
  return (X)
  
}



 
#' generate L stats
#'
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @import spatstat.geom
#' @import parallel
#' 
#' @export
run_L_function <- function( data, type = "spatial_p" , ncores = 8 ){
  
  
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
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @import spatstat.geom
#' @import parallel
#' 
#' @export
run_celltype_interaction <- function( data, type = "spatial_p" , ncores = 8  ){
  
  
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
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @import spatstat.geom
#' @import parallel
#' @import ape
#' 
#' @export
run_Morans_I <- function( data, type = "spatial_p" , ncores = 8  ){
  
  
  if ( type == "scrna" )  {
    print("This feature class currently does not support scRNA-seq")
    return(NULL)
  }
  
  
  if ( type %in% c( "spatial_p" , "spatial_t" ))   {
    X <- helper_moran( data , num_top_gene =  NULL,   ncores = 8)
  }
  
  return (X)
  
}









#' generate Moran's I
#'
#' @param data input data
#' @param type input data type, either scrna, spatial_p, or spatial_t
#' @param ncores number of cores 
#' 
#' @return a matrix of samples x features 
#' 
#' @import spatstat.geom
#' @import spatstat.core
#' @import parallel
#' 
#' @export
run_nn_correlation <- function( data, type = "spatial_p", ncores = 8  ){
  
  
  if ( type == "scrna" )  {
    print("This feature class currently does not support scRNA-seq")
    return(NULL)
  }
  
  
  
  if ( type %in% c( "spatial_p" , "spatial_t" ))   {
    X <- helper_nncorr_protein(data, num_top_gene  =  NULL , ncores = ncores  )
  }
  
  
  
  return (X)
  
}








