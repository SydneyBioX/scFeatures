 
# helper function to run CCI 
helper_CCI <- function( alldata , ncores = 1  ){
 
  BPparam <-  generateBPParam(ncores)
 
  data(LRdb, package = "SingleCellSignalR")
  
  # x <- unique(alldata$sample)[1]
 
  capture.output( suppressMessages( individual_cci <- BiocParallel::bplapply(  unique(alldata$sample), function(x){
    
                        data_dataframe  <- alldata$data[, alldata$sample == x]
                        
                        celltype <- as.factor( alldata$celltype[ alldata$sample == x])
                        celltype_numeric <- as.numeric(  celltype)
                        
                        signal <-  SingleCellSignalR::cell_signaling(data = data_dataframe,
                                                                     genes = rownames(data_dataframe), 
                                                                    cluster =   celltype_numeric,
                                                                    c.names = levels(celltype), write = FALSE)
                        
                        if ( length(signal) == 0 ){
                           all_interaction <- data.frame(LRscore = 0, feature = "placeholder" )
                        }else{
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
                        }
                       
                        all_interaction
               }, BPPARAM = BPparam) ) )
   
   
   # gather the cell - cell interaction probability into sample x interaction probability matrix 
   X <- NULL
   for (i in c(1:length( individual_cci))){
     temp <-   individual_cci[[i]]
     temp <- temp[ !duplicated(temp$feature) ,]
   
     if (is.null(X)){
       X <- temp
     }else{
       X <- merge(X, temp , by="feature", all = TRUE)
     }
   }
   rownames(X) <- X$feature
   X <- X[!rownames(X) == "placeholder", ]
   X <- X[, -1]
   colnames(X) <-  unique(alldata$sample)
   X[is.na(X)] <- 0
   
   X <- t(X)
 
   return(X)
   
}

 