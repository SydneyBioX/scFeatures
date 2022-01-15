

# create pseudo-bulk for each cell type of each sample 
bulk_sample_celltype <- function(data , ncores = 8  ){
  
  bulk <- mclapply( unique( data$sample),  function(x) {
      
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
        temp
      })
      
      this_sample_bulk <- as.data.frame(do.call(cbind, this_sample_bulk))
      
      colnames( this_sample_bulk) <- paste0(x , "--" , unique(data$celltype) )
      this_sample_bulk
  
  }, mc.cores = ncores)
  
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
bulk_sample  <- function(data,  ncores = 8 ){
  
  bulk <- mclapply( unique(data$sample), function(x ) {
    index <-  which(data$sample == x )
    temp <-  DelayedMatrixStats::rowMeans2(DelayedArray( data@assays$RNA@data[, index]))
    temp <- as.matrix(temp)
  }, mc.cores= ncores)
  
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
#' @import gtools
#' @import tidyr
#' @import parallel
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
  num_cell_per_spot  <- mapply(FUN = "*", as.data.frame(      prob ), MultVec)
  
  
  num_cell_per_spot <- round(num_cell_per_spot, 0)
  mode( num_cell_per_spot ) <- "integer"
  
  rownames(num_cell_per_spot ) <- rownames(  prob)
  
  return( num_cell_per_spot )
  
}




 


L_stats <- function(ppp_obj = NULL, from = NULL, to = NULL, L_dist = NULL) {
  L <-  spatstat.core::Lcross(ppp_obj, 
                              from = from,
                              to = to,
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



