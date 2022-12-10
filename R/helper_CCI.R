 
# helper function to run CCI 
helper_CCI <- function( data , ncores = 1  ){
  
  
  BPparam <- generateBPParam(ncores)
  
  # x <- unique(data$sample)

  individual_cci <- BiocParallel::bplapply(  unique(data$sample), function(x){
    
                        thissample <- data[, data$sample == x]
                        data_dataframe <- thissample@assays$RNA@data
                        celltype <- as.factor(thissample$celltype)
                        celltype_numeric <- as.numeric(  celltype)
                        
                        signal <- SingleCellSignalR::cell_signaling(data = data_dataframe, genes = rownames(data_dataframe), 
                                                                   cluster =   celltype_numeric,
                                                                   c.names = levels(celltype), write = FALSE)
                        
                        # concat interaction from each cell type
                        all_interaction <- NULL
                        for ( i in 1:length(signal)){
                            this_celltype <- signal[[i]]
                            this_celltype$feature <- paste0( colnames( this_celltype )[1] , "->" , colnames( this_celltype )[2],
                                                             "--", 
                                                             this_celltype[, 1]  , "->", this_celltype[, 2])
                            this_celltype <-   this_celltype[, c("LRscore", "feature")]
                            all_interaction <- rbind( all_interaction,    this_celltype )
                        }
                        all_interaction
                     }, BPPARAM = BPparam) 
   
   
   # gather the cell - cell interaction probability into sample x interaction probability matrix 
   X <- NULL
   for (i in c(1:length(individual_cci))){
     temp <-   individual_cci[[i]]
     temp <- temp[ !duplicated(temp$feature) ,]
   
     if (is.null(X)){
       X <- temp
     }else{
       X <- merge(X, temp , by="feature", all = T)
     }
   }
   rownames(X) <- X$feature
   X <- X[, -1]
   colnames(X) <-  unique(data$sample)
   X[is.na(X)] <- 0
   
   X <- t(X)
   
   return(X)
   
}

 