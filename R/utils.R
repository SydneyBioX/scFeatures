

# check that all metadata are in the data
check_data <- function(data, type = "scrna"){
  
  if (class(data) != "Seurat"){
    print("please make sure the data is in a Seurat object")
  }
 
  if (  !"sample" %in% names(data@meta.data)){
      print("For scRNA-seq and spatial proteomics, please make sure the data contains celltype and sample label.
            For spatial proteomics, please make sure the data contains sample information.")
      stop()
  }
  
  if (type %in% c( "spatial_t" , "spatial_p")){
    if ( !"y_cord"  %in% names(data@meta.data) ||  !"x_cord" %in% names(data@meta.data) ) {
     print("please make sure the data contains x_cord and y_cord")
     stop()
   }
  }
  
  if (type == "spatial_t"){
    if ( ! "predictions" %in% names(data@assays) ) {
      print("please make sure the data contains a predictions assay")
      print("see vignette's section on spatial transcriptomics for more explanation")
      stop()
    }
  }

}






generateBPParam <- function(cores = 1){
  
  set.seed(1)
  seed <- .Random.seed[1]
  
  if(cores == 1){
    BPparam <- BiocParallel::SerialParam(RNGseed = seed)
  } else { # Parallel processing is desired.
    # Also set the BPparam RNGseed if the user ran set.seed(someNumber) themselves.
    if(Sys.info()["sysname"] == "Windows") {# Only SnowParam suits Windows.
      BPparam <- BiocParallel::SnowParam(min(cores, BiocParallel::snowWorkers("SOCK")), RNGseed = seed)
    } else if (Sys.info()["sysname"] %in% c("MacOS", "Linux")) {
      BPparam <- BiocParallel::MulticoreParam(min(cores, BiocParallel::multicoreWorkers()), RNGseed = seed) # Multicore is faster than SNOW, but it doesn't work on Windows.
    } else { # Something weird.
      BPparam <- BiocParallel::bpparam() # BiocParallel will figure it out.
    }
  }
  
  return(  BPparam  )
  
}
  
  
  
 

# create pseudo-bulk for each cell type of each sample 
bulk_sample_celltype <- function(data , ncores = 1  ){
  
  BPparam <- generateBPParam(ncores)
  
  # x <- unique( data$sample)[1]
  bulk <- BiocParallel::bplapply( unique( data$sample),  function(x) {
      
      # for this patient
      this_sample <- data[ , data$sample == x]
      
      # loop through each cell type 
      this_sample_bulk <- lapply( unique(data$celltype), function(y) {

        index <-  which(this_sample$celltype == y )
        
        # if this cell type does not exist in this patient, expression is 0 for all genes 
        if(length(index) == 0){ 
          temp <-  rep(0, nrow(data)) 
        # if there is only 1 cell, do not need to take mean
        }else if (length(index) == 1){ 
          temp <-    this_sample@assays$RNA@data[, index] 
        # if multiple cells, average across all cells 
        }else{
          temp <-  DelayedMatrixStats::rowMeans2( DelayedArray( this_sample@assays$RNA@data[, index] ) )
        }
        
        temp <- as.matrix(temp)
        rownames(temp ) <- rownames(data)
        temp
      })
      
      this_sample_bulk <- as.data.frame(do.call(cbind, this_sample_bulk))
      
      colnames( this_sample_bulk) <- paste0(x , "--" , unique(data$celltype) )
      this_sample_bulk
  
  }, BPPARAM =  BPparam) 
  
  bulk <- as.data.frame(do.call(cbind,  bulk))
  bulk <- CreateSeuratObject( bulk )
 
  sample <- unlist( lapply( strsplit(  colnames(bulk), "--") , `[`, 1))
  celltype <- unlist( lapply( strsplit(  colnames(bulk), "--") , `[`, 2))
  #condition <- unlist( lapply( strsplit(  colnames(bulk), "--") , `[`, 3))
  bulk$sample <- sample
  bulk$celltype <- celltype
  #bulk$condition <- condition 
  
  return( bulk )
  
}



# create pseudo-bulk for each sample 
bulk_sample  <- function(data,  ncores = 1 ){
  
  BPparam <- generateBPParam(ncores)
  
  bulk <-  BiocParallel::bplapply(   unique(data$sample), function(x ) {
    index <-  which(data$sample == x )
    temp <-  DelayedMatrixStats::rowMeans2(DelayedArray( data@assays$RNA@data[, index]))
    temp <- as.matrix(temp)
  }, BPPARAM =  BPparam) 
  
  bulk <- as.data.frame(do.call(cbind,  bulk ))
  rownames(bulk) <- rownames(data)
  colnames(  bulk ) <-   unique(data$sample)
  bulk <- CreateSeuratObject(bulk)
  
  bulk$sample <- unique(data$sample)
 # bulk$condition <- data[ , match( bulk$sample, data$sample)]$condition 

  return( bulk )
}



#' estimate a relative number of cells per spot 
#'
#' @param data input data
#' 
#' @return data with the relative number of cells per spot stored in the meta.data
#' 
#' 
#' @export
get_num_cell_per_spot <- function(data){
  
  readcount <- log2( colSums( data ) )
  
  linMap <- function(x, from, to)
    (x - min(x)) / max(x - min(x)) * (to - from) + from
  
  numberofcells <- linMap(readcount, 1, 100)
  
  data$number_cells <- numberofcells
  
  return ( data )
}


rearrange_string <- function(str) {
  unlist(lapply(strsplit(str, "_"), function(x) paste(sort(x), collapse = "_")))
}



# get number of cells in each cell type in each spot 
# calculated by cell type probability in each spot times the relative number of cells in each spot 
# relative number of cells are estimated using library size of each spot 
get_num_cell_per_celltype <- function(data){
  
  prob <- data@assays$predictions
  prob <- as.matrix(prob@data)
  prob <- prob[ !rownames(prob ) == "max", ]
  zero_celltype <- names( which ( rowSums ( prob )  == 0 ))
  prob <-  prob[ !rownames(prob) %in% zero_celltype , ]
  
  MultVec <-  data$number_cells
  num_cell_per_spot  <- mapply(FUN = "*", as.data.frame( prob ), MultVec)
  
  
  num_cell_per_spot <- round(num_cell_per_spot, 0)
  mode( num_cell_per_spot ) <- "integer"
  
  rownames(num_cell_per_spot ) <- rownames(  prob)
  
  return( num_cell_per_spot )
  
}




 


L_stats <- function(ppp_obj = NULL, from = NULL, to = NULL, L_dist = NULL) {
  L <-  spatstat.core::Lcross(ppp_obj,  from = from, to = to,
                              verbose = FALSE,
                              correction = "best")
  
  
  L_theo <- L$theo[L$r <= L_dist]
  L_iso <- L$iso[L$r <= L_dist]
  L_res <- mean(L_iso - L_theo)
  
  return(L_res)
}





#' perform pre-processing
#'
#' @param data input data
#' @param normalise whether to normalise the data. Note if the data has already been normalised (eg, log2CPM), there is no need to normalise again. 
#' 
#' @return a matrix of samples x features 
#' 
#' @import Seurat
#' 
#' @export
process_data <- function(data, normalise = T){
  
  if (!is.null(data@meta.data$celltype) ){
    
    data$celltype <- gsub("\\+", "plus", data$celltype)
    data$celltype <- gsub("\\-", "minus", data$celltype)
    
    data$celltype <- as.character(data$celltype)
    
    # remove "small" patient that has less than 10 cells across all cell types 
    celltype_per_patient <- table(data$celltype, data$sample)
    a  <- celltype_per_patient <= 10
    a <-  colSums(a)
    small_patient <- names( which (a == length(unique(data$celltype))))
    
    if (length(small_patient) > 0 ){
      data <- data[, -c(which(data$sample %in% small_patient))]
    }
    
  }
 
  if (!is.null(data@meta.data$sample) ){
    data$sample <- as.character(data$sample)
  }
  
  if (!is.null(data@meta.data$condition)){
    data$condition <- as.character(data$condition)
 
  }
 
  
  if(normalise == T){
    data <- Seurat::NormalizeData(data)
  }

  return(data)
}



#' automatically generate the association study report in an html format
#'
#' @param scfeatures_list a named list storing the feature output from scfeatures
#' @param output_folder directory for saving the html report file
#' 
#' @return the html file will be saved in the directory defined in the output_folder
#' 
#' @import rmarkdown
#' 
#' @export
run_association_study_report <- function( scfeatures_result, output_folder ){
  # check name
  correct_name <- any( names(scfeatures_result) %in%  c(
     "proportion_raw"  , "proportion_logit"   , 
      "proportion_ratio",  "gene_mean_celltype"  , 
     "gene_prop_celltype" , "gene_cor_celltype",   
     "pathway_gsva" ,"pathway_mean" ,       
      "pathway_prop" ,  "CCI"   ,
     "gene_mean_aggregated", "gene_cor_aggregated"  ,"gene_prop_aggregated"))
  if (correct_name == FALSE){
    print("Please check you have named the feature types in correct naming format.")
  }
  
  # need to retrieve the output report structure from the package
  output_report <- system.file("extdata",   "output_report.Rmd",   
                               package = "scFeatures")
  
  file.copy(from = output_report, to = output_folder, overwrite = FALSE)
  
  # generate the html output 
  rmarkdown::render(input = "output_report.Rmd", 
                    output_format = "html_document", 
                    output_file = "output_report.html")
}


