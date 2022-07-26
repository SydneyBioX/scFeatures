


individual_celltype_interaction_sp <- function(this_sample   ){
  
 
  cell_points <- spatstat.geom::ppp(x = this_sample$x_cord , 
                     y =  this_sample$y_cord , 
                     check = FALSE,
                     yrange = c( 0 , max(this_sample$y_cord  ) ), 
                     xrange = c(0,  max(this_sample$x_cord) ),
                     marks = as.factor( this_sample$celltype))
  
  tab <- table(this_sample$celltype)
  
  cellTypes_toTest <- names(tab[which(tab > 10)])
  cellTypes_pair <- expand.grid(cellTypes_toTest, cellTypes_toTest,
                                stringsAsFactors = FALSE)
 
  
  # Calcualte the pairwise distance
  d <- spatstat.geom::pairdist(cell_points, squared = FALSE)
  diag(d) <- Inf
  
  nn_list <- apply(d, 1, function(x) which(x < 50))
  
  nn_list_cellTypes <- lapply(seq_along(nn_list), function(idx) {
    if (length(nn_list[[idx]]) > 0) {
      paste( this_sample$celltype [idx],
             this_sample$celltype[nn_list[[idx]]], sep = "_")
    }
  })
  
  nn_list_cellTypes <- unlist(nn_list_cellTypes)
  nn_list_cellTypes <- rearrange_string(nn_list_cellTypes)
  nn_list_cellTypes <- table(nn_list_cellTypes)
 
  return(  nn_list_cellTypes  )
}





helper_celltype_interaction_sp <- function(data,  ncores = 1){
  
  BPparam <- generateBPParam(ncores)
  
  s <- unique(data$sample)[1]
 
  nn_list_cellTypes <- BiocParallel::bplapply( unique(data$sample) , function(s) {
 
    this_sample <- data[ , data$sample == s]
    nn_list_cellTypes <- individual_celltype_interaction_sp(  this_sample )
    
  }, BPPARAM =  BPparam)
  
  
  temp <- NULL
  
  for( i in c(1: length(  nn_list_cellTypes ))){
    
  
    err <-  try({
      
      a <- nn_list_cellTypes[[i]] 
      a <-  data.frame( a ) 
      a$Freq <-  a$Freq / sum (a$Freq)
      
      if( is.null(temp)){
        temp <- a
      }else{
        temp <-  suppressWarnings( merge(temp, a, by="nn_list_cellTypes", all = T) ) 
        } 
    }, silent = T )
    
    
    if ( class(err) == "try-error") {
      a <- data.frame( rep(0, nrow(temp)))
      a$nn_list_cellTypes <- temp$nn_list_cellTypes
      temp <- suppressWarnings( merge( temp , a ,   by="nn_list_cellTypes", all = T ) ) 
      
    }
    
    colnames(temp ) <- make.names(colnames(temp), unique=T)
    
  }
  
  
  rownames(temp) <- temp$nn_list_cellTypes
  temp <- temp[, -1]
  
  colnames(  temp ) <- unique(data$sample)
  
  temp  <- t( temp )
  
  temp[is.na(temp)] <- 0 
  nn_list_cellTypes <- temp
  
  
  return (  nn_list_cellTypes  )
  
  
}






individual_celltype_interaction_st <- function(  thisprob ){
  
 
  x <- 1
  temp <- lapply( 1:ncol( thisprob) , function(x) {
    
    thisspot <- thisprob[ , x]
    thisspot <- thisspot %*% t(thisspot )
    rownames(thisspot ) <- colnames(thisspot)
    thisspot <- reshape2::melt(     thisspot)
    a <- data.frame(thisspot$value)
    rownames(a) <- paste0(thisspot$Var1 , "-with-", thisspot$Var2)
    a
  })
  
  temp <- do.call(cbind, temp)
  nn_list_cellTypes <- rowSums(temp)
  
  
  return (   nn_list_cellTypes )
}




helper_celltype_interaction_st <- function(data,  ncores = 1){
  
  BPparam <- generateBPParam(ncores)
  
  s <- unique(data$sample)[1]
  
  prob <- data@assays$predictions
  prob <- as.matrix(prob@data)
  prob <- prob[ !rownames(prob ) == "max", ]
  zero_celltype <- names( which ( rowSums ( prob )  == 0 ))
  prob <-  prob[ !rownames(prob) %in% zero_celltype , ]
  
  
  nn_list_cellTypes <- BiocParallel::bplapply( unique(data$sample) , function(s) {
    index <- which(data$sample == s)
    thisprob <- prob[, index]
    nn_list_cellTypes <- individual_celltype_interaction_st(    thisprob  )
    
  }, BPPARAM =  BPparam)

  
  nn_list_cellTypes <- do.call(cbind,   nn_list_cellTypes)
  colnames(  nn_list_cellTypes  ) <- unique(data$sample)
  
  nn_list_cellTypes <- t(nn_list_cellTypes)
  
  return (  nn_list_cellTypes  )
  
  
}








 


individual_L_stat_st  <- function( thissample  , this_num_cell_per_spot ){
  
  
  # expand each spot into its number of cells 
  x <- c()
  y <- c()
  
  celltype <- c()
  
  i <-1 
  
  
  gap_x <- (  max (thissample$x_cord ) -   min ( thissample$x_cord) ) /length(thissample$x_cord) / 2
  gap_y <- (  max (thissample$y_cord ) -   min ( thissample$y_cord) ) /length(thissample$y_cord) / 2
  
  
  for ( i in 1:ncol(thissample)){
    
    thisspot <- thissample[ , i]
    thisspot_num_cell <- this_num_cell_per_spot[, i ]
    
    total_num_cell <- sum(thisspot_num_cell)
    
    this_gap_x <-    gap_x / total_num_cell
    this_gap_y <-    gap_y / total_num_cell
    
    x <- c( x ,   seq( thissample$x_cord[i],   by =   this_gap_x , length.out =     total_num_cell  )   )
    y <- c( y,   seq( thissample$y_cord[i],   by = this_gap_y, length.out =     total_num_cell  ) )
    
    celltype <- c( celltype, rep(rownames(this_num_cell_per_spot ),  thisspot_num_cell ))
    
  }
  
  
  cell_points_threecelltype <- spatstat.geom::ppp(x = x , 
                                   y =  y , 
                                   check = FALSE,
                                   yrange = c( min (as.numeric( thissample$y_cord )), 
                                               max (as.numeric( thissample$y_cord ))  ) , 
                                   xrange = c( min (as.numeric( thissample$x_cord )), 
                                               max (as.numeric( thissample$x_cord ))  ) , 
                                   marks =   as.factor(   celltype ))
  
  
  cellTypes_toTest <- names ( which(table(celltype) > 10))
  cellTypes_pair <- expand.grid(cellTypes_toTest, cellTypes_toTest,
                                stringsAsFactors = FALSE)
  
  
  
  L_patient  <- lapply ( 1:nrow(cellTypes_pair), function(i) {
    L_stats(    cell_points_threecelltype  , 
                from = cellTypes_pair[i, 1],
                to = cellTypes_pair[i, 2],
                L_dist = 4)
  })
  
  L_patient <- do.call(c, L_patient)
  names(L_patient) <- paste(cellTypes_pair[, 1], cellTypes_pair[, 2], sep = "_")
  
  return(L_patient)
  
  
}




helper_L_stat_st <- function(data, ncores = 1){
 
  BPparam <- generateBPParam(ncores)
  
  num_cell_per_spot  <- get_num_cell_per_celltype(data)
  
  s <- unique(data$sample)[1]
  
  
  L_stats  <- BiocParallel::bplapply( unique(data$sample), function(s) {
    
    index <- which( data$sample == s)
    
    thissample <- data[ ,    index ]
    this_num_cell_per_spot <- num_cell_per_spot[ , index]
    
    L_patient <- individual_L_stat_st(thissample  , this_num_cell_per_spot )
    
    
  }, BPPARAM =  BPparam)
  
  
  temp <- NULL
  
  for( i in c(1:length(L_stats))){
  
    err <-  try({
      
      a <-   L_stats[[i]] 
      a <- data.frame( a)
      a$rowname <- rownames(a)
      
      if( is.null(temp)){
        temp <- a
      }else{
        temp <-  suppressWarnings(  merge(temp, a, by="rowname", all = T))
        }
    }, silent = T)
    
    
    if ( class(err) == "try-error") {
      a <- data.frame( rep(0, nrow(temp)))
      a$rowname <- temp$rowname
      temp <-  suppressWarnings( merge( temp , a ,   by="rowname", all = T ) ) 
      
    }
    
    colnames(temp ) <- make.names(colnames(temp), unique=T)
  }
  
  
  rownames(temp) <- temp$rowname
  temp <- temp[, -1]
  
  colnames(  temp ) <- unique(data$sample)
  
  temp  <- t( temp )
  
  temp[is.na(temp)] <- 0 
  L_patient <- temp
  
  return (  L_patient )
  
}




individual_L_stat_sp <- function( this_sample){
 
 
  cell_points <- spatstat.geom::ppp(x = this_sample$x_cord , 
                     y =  this_sample$y_cord , 
                     check = FALSE,
                     yrange = c( 0 , max(this_sample$y_cord  ) ), 
                     xrange = c(0,   max(this_sample$x_cord) ),
                     marks = as.factor( this_sample$celltype))
  
  tab <- table(this_sample$celltype)
  
  cellTypes_toTest <- names(tab[which(tab > 10)])
  cellTypes_pair <- expand.grid(cellTypes_toTest, cellTypes_toTest,
                                stringsAsFactors = FALSE)
  
  
  
  L_patient <- list()
  for (i in 1:nrow(cellTypes_pair)) {
 
    L_patient[[i]] <- L_stats(cell_points, 
                              from = cellTypes_pair[i, 1],
                              to = cellTypes_pair[i, 2],
                              L_dist = 50)
  }
  
  L_patient <- do.call(c, L_patient)
  names(L_patient) <- paste(cellTypes_pair[, 1], cellTypes_pair[, 2], sep = "_")
  
  return (L_patient )
  
  
}



helper_L_stat_sp <- function(data, ncores = 1){
  
  
  BPparam <- generateBPParam(ncores)
  
  s <- unique(data$sample)[1]
  
  
  L_stats_result  <- BiocParallel::bplapply( unique(data$sample), function(s) {
    
    index <- which( data$sample == s)
    
    thissample <- data[ ,    index ]
   
    L_patient <- individual_L_stat_sp(thissample   )
    
    
  } , BPPARAM =  BPparam)
  
  
  temp <- NULL
  
  for( i in c(1:length(L_stats_result))){
    
    err <-  try({
      
      a <-   L_stats_result[[i]] 
      a <- data.frame( a)
      a$rowname <- rownames(a)
      
      if( is.null(temp)){
        temp <- a
      }else{
        temp <-  suppressWarnings( merge(temp, a, by="rowname", all = T) ) 
        } 
    }, silent = T)
    
    
    if ( class(err) == "try-error") {
      a <- data.frame( rep(0, nrow(temp)))
      a$rowname <- temp$rowname
      temp <- suppressWarnings( merge( temp , a ,   by="rowname", all = T ) ) 
      
    }
    
    colnames(temp ) <- make.names(colnames(temp), unique=T)
  }
  
  
  rownames(temp) <- temp$rowname
  temp <- temp[, -1]
  
  colnames(  temp ) <- unique(data$sample)
  
  temp  <- t( temp )
  
  temp[is.na(temp)] <- 0 
  L_patient <- temp
  
  return (  L_patient )
  
}










individual_nncorr_protein <- function(  thissample ){
  
  exprsMat <- thissample@assays$RNA@data
  
  cell_points_cts <- spatstat.geom::ppp(x = as.numeric( thissample$x_cord) , y = as.numeric( thissample$y_cord) ,
                         check = FALSE,
                         xrange = c( min (as.numeric(thissample$x_cord ) ), 
                                     max (as.numeric(thissample$x_cord)  )) , 
                         yrange = c( min (as.numeric(thissample$y_cord)) , 
                                     max ( as.numeric(thissample$y_cord )  )) , 
                         marks = t(as.matrix( exprsMat)))
  err <- try({
    nncorr_protein <- spatstat.core::nncorr(cell_points_cts)["correlation", ]
  }, silent = T)
  if (class(err) == "try-error"){
    nncorr_protein <- NA
  }
  
  return(  nncorr_protein )
  
}





helper_nncorr_protein <- function(data, num_top_gene  =  NULL , ncores = 1 ){
  
  BPparam <- generateBPParam(ncores)
  
  if ( is.null( num_top_gene ) ){
    num_top_gene = min(nrow(data), 1500 ) 
  }
  
  top_gene <- find_var_gene(data,  num_top_gene  =  num_top_gene , 
                            ncores = ncores , celltype = F)
  
  data@assays$RNA@data <-   data@assays$RNA@data[ rownames( data@assays$RNA@data) %in% top_gene, ]
 
  
  s <- unique(data$sample)[1]
  
  nncorr_protein  <- BiocParallel::bplapply( unique(data$sample) , function(s) {
    thissample <- data[ , data$sample == s]
    L_patient <- individual_nncorr_protein(  thissample  )
    
  } , BPPARAM =  BPparam)
  
  
  temp <- NULL
  
  for( i in c(1:length( nncorr_protein))){

    err <-  try({
      
      a <-  nncorr_protein[[i]] 
      a <- data.frame( a)
      a$rowname <- rownames(a)
      
      if( is.null(temp)){
        temp <- a
      }else{
        temp <-  suppressWarnings( merge(temp, a, by="rowname", all = T) )
        } 
    }, silent = T)
    
    
    if ( class(err) == "try-error") {
      a <- data.frame( rep(0, nrow(temp)))
      a$rowname <- temp$rowname
      temp <- suppressWarnings( merge( temp , a ,   by="rowname", all = T ) ) 
      
    }
    
    colnames(temp ) <- make.names(colnames(temp), unique=T)
    
  }
  
  
  rownames(temp) <- temp$rowname
  temp <- temp[, -1]
  
  colnames(  temp ) <-   unique(data$sample)
  
  temp  <- t( temp )
  
  temp[is.na(temp)] <- 0 
  nncorr_protein <- temp
  
  
  return (  nncorr_protein  )
   
  
}





individual_moran_cor <- function(    thissample ){
  
  exprsMat <- thissample@assays$RNA@data
  
  cell_points_cts <- spatstat.geom::ppp(x = as.numeric( thissample$x_cord) , y = as.numeric( thissample$y_cord ),
                         check = FALSE,
                         xrange = c( min (as.numeric(thissample$x_cord )), 
                                     max (as.numeric(thissample$x_cord)  )) , 
                         yrange = c( min (as.numeric(thissample$y_cord)), 
                                     max ( as.numeric( thissample$y_cord )  )) , 
                         marks = t(as.matrix( exprsMat)))
  
  d <- spatstat.geom::pairdist(cell_points_cts, squared = FALSE)
  diag(d) <- Inf
  
  
  w <- 1/d
  
  
  moran_cor <-  lapply( 1:nrow(exprsMat), function(x){
    err <- try( val <- ape::Moran.I(exprsMat[x, ], w)$observed, silent = TRUE)
    if (is( err, "try-error")) {
      NA
    }else{
      val
    }
  } )
  
  names(moran_cor) <- rownames(exprsMat)
  moran_cor <- unlist(moran_cor)
  
  return(    moran_cor )
  
}





helper_moran <- function(data,  num_top_gene =  NULL, ncores = 1){
  
  BPparam <- generateBPParam(ncores)
  
  if ( is.null( num_top_gene ) ){
    num_top_gene = min(nrow(data), 1500 ) 
  }
  
  top_gene <- find_var_gene(data,  num_top_gene  =  num_top_gene , 
                            ncores = ncores , celltype = F)
  
  data@assays$RNA@data <-  data@assays$RNA@data[ rownames( data@assays$RNA@data) %in% top_gene, ]
  
  
  s <- unique(data$sample)[1]
  
  moran_cor  <- BiocParallel::bplapply( unique(data$sample) , function(s) {
    thissample <- data[ , data$sample == s]
    moran_cor  <- individual_moran_cor (   thissample)
    
  } , BPPARAM =  BPparam)
   
  temp <- NULL
  
  for( i in c(1:length( moran_cor))){
   
    err <-  try({
      
      a <- moran_cor[[i]] 
      a <- data.frame( a)
      a$rowname <- rownames(a)
      
      if( is.null(temp)){
        temp <- a
      }else{
        temp <-  suppressWarnings( merge(temp, a, by="rowname", all = T) ) 
        } 
    }, silent = T)
    
    
    if ( class(err) == "try-error") {
      a <- data.frame( rep(0, nrow(temp)))
      a$rowname <- temp$rowname
      temp <- suppressWarnings( merge( temp , a ,   by="rowname", all = T ) ) 
    }
    
    colnames(temp ) <- make.names(colnames(temp), unique=T)
  }
  
  
  rownames(temp) <- temp$rowname
  temp <- temp[, -1]
  
  colnames(  temp ) <- unique(data$sample)
  
  temp  <- t( temp )
  
  temp[is.na(temp)] <- 0 
  moran_cor <- temp
  
  
  return(  moran_cor  )
}












