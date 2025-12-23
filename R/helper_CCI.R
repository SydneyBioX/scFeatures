 
# helper function to run CCI 
helper_CCI <- function( alldata , species, ncores = 1  ){
 
  BPparam <-  generateBPParam(ncores)
 
  
  if (species == "human"){
    CellChatDB <- CellChatDB.human
  }else{
    CellChatDB <- CellChatDB.mouse
  }
  
  
  # x <- unique(alldata$sample)[1]
 
  capture.output( suppressMessages( individual_cci <- BiocParallel::bplapply(  unique(alldata$sample), function(x){
    
              this_sample_data <- alldata$data[, alldata$sample == x]
              colnames(  this_sample_data) <- make.unique( colnames(   this_sample_data) ) 
              
              
              meta = data.frame(labels = alldata$celltype[ alldata$sample == x]   )   
              rownames(meta) <-  colnames(   this_sample_data) 
          
              cellchat  <- createCellChat(object =  this_sample_data , meta = meta, 
                                          group.by = "labels")
          
              
              cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
              
              
              
              groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
           
              cellchat@DB <- CellChatDB # set the used database in the object
              
              cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
              # do parallel
              cellchat <- identifyOverExpressedGenes(cellchat)
              cellchat <- identifyOverExpressedInteractions(cellchat)
              cellchat <- projectData(cellchat, PPI.human)
              cellchat <- computeCommunProb(cellchat)
              cellchat <- computeCommunProbPathway(cellchat)
              cellchat <- aggregateNet(cellchat)
              
              cellchat_score <-   netVisual_bubble(   cellchat,   return.data = TRUE )   
              cellchat_score  <- cellchat_score $communication
              
     
              if ( nrow(  cellchat_score ) == 0 ){
                cellchat_score <- data.frame(LRscore = 0, feature = "placeholder" )
              }else{
                cellchat_score$feature <- paste0(  cellchat_score$source   , "->" , 
                                                   cellchat_score$target,
                                                   "--", 
                                                  cellchat_score$ligand  , "->", 
                                                  cellchat_score$receptor)

                }
             
              
              
              cellchat_score 
              
              
        }, BPPARAM = BPparam) ) )
  
  
  
  
   
   
   # gather the cell - cell interaction probability into sample x interaction probability matrix 
   X <- NULL
   for (i in c(1:length( individual_cci))){
     temp <-   individual_cci[[i]]
     
     temp  <-    temp [, c("feature", "prob" )]
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

 
