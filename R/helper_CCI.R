 
# helper function to run CCI 
helper_CCI <- function( alldata , species, ncores = 1  ){
 
  BPparam <-  generateBPParam(ncores)
 
  data("CellChatDB.human", package="CellChat")
  data("CellChatDB.mouse", package="CellChat")
  data("PPI.human", package="CellChat")
  data("PPI.mouse",package="CellChat" )
  
  if (species == "Homo sapiens"){
    CellChatDB <- CellChatDB.human
    PPI <- PPI.human
  }else{
    CellChatDB <- CellChatDB.mouse
    PPI <- PPI.mouse
  }
  
  
  # x <- unique(alldata$sample)[2]
 
 
  capture.output( suppressMessages( individual_cci <- BiocParallel::bplapply(  unique(alldata$sample), function(x){
           
 
    
    err <- try({
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
              cellchat <- smoothData(cellchat, adj = PPI)
              
            
              cellchat <- computeCommunProb(cellchat)
              cellchat <- computeCommunProbPathway(cellchat)
              cellchat <- aggregateNet(cellchat)
              
              
          
              cellchat_score <-   netVisual_bubble(   cellchat,   return.data = TRUE )   
              cellchat_score  <- cellchat_score $communication
              
            
              
                cellchat_score$feature <- paste0(  cellchat_score$source   , "->" , 
                                                   cellchat_score$target,
                                                   "--", 
                                                  cellchat_score$ligand  , "->", 
                                                  cellchat_score$receptor)

               
           
              cellchat_score 
              
             
    })
    
    if (class(err) == "try-error"){
 
        cellchat_score <- data.frame(source = "placeholder", 
                                     target = "placeholder",
                                     ligand = "placeholder",
                                     receptor = "placeholder",
                                     prob =  0,
                                     pval = 0,
                                     interaction_name = "placeholder",
                                     interaction_name_2 = "placeholder",
                                     pathway_name = "placeholder",
                                     annotation =  "placeholder",
                                     evidence =  "placeholder",
                                     source.target =  "placeholder",
                                     prob.original =  0,
                                     feature = "placeholder" )
        
        
        cellchat_score 
    
    }else{
      cellchat_score 
      }
  
       }, BPPARAM = BPparam) ))
  
  
  
  
   
   
   # gather the cell - cell interaction probability into sample x interaction probability matrix 
   X <- NULL
   for (i in c(1:length( individual_cci))){
     temp <-   individual_cci[[i]]
     
     temp  <-    temp [, c("feature", "prob" )]
     temp <- temp[ !is.na(temp$prob), ]
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

 
