





# get hallmark geneset if the user does not provide one 
get_geneset <- function(species = "Homo sapiens"){
  
  m_df <- msigdbr::msigdbr(species = species , category = "H")
  m_t2g <- m_df %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
  geneset <- split(x = m_t2g$entrez_gene, f = m_t2g$gs_name)
  
  # convert the entrez id to gene symbol 
  geneset  <-  lapply( geneset, function(x){
    
    if (species == "Homo sapiens"){
      geneID  <- ensembldb::select(EnsDb.Hsapiens.v79 , 
                                        keys= as.character(x) , 
                                        keytype = "ENTREZID",  columns = c("SYMBOL", "GENEID"))
      
    }else{
      geneID <- ensembldb::select(EnsDb.Mmusculus.v79 , 
                                        keys= as.character(x) , 
                                        keytype = "ENTREZID",  columns = c("SYMBOL", "GENEID"))
    }
    this_geneset <- unique ( geneID$SYMBOL)
    this_geneset
  })
  
  return(   geneset )
}




# helper function to run the pathway gsva 
helper_pathway_gsva <- function(data,  method = "ssgsea"  , geneset, ncores = 1 ){
  
  if (method == "ssgsea"){
    
      # if the dataset has greater than 30000 cells 
      #  then it is actually too large to be computed in one go in gsva
      #  split into multiple set 
      if (ncol(data) > 30000){
        
        index <- seq(1, ncol(data), by = 30000)
        index <- c( index, ncol(data))
        topMatrixGSVA <- NULL
        
        for (i in c(2:length(index))){
          start = index[i-1]
          finish = index[i ] - 1
          
          print(paste0("calculating ", start, " to ", finish , " cells" ))
          thesecell  <- as.matrix( data@assays$RNA@data[, start: finish])
          temp_topMatrixGSVA <- GSVA::gsva(thesecell  ,  geneset,  method="ssgsea", 
                                     min.sz=10, max.sz=999999, abs.ranking=F, verbose=T,
                                     parallel.sz=ncores)
          topMatrixGSVA  <- cbind( topMatrixGSVA , temp_topMatrixGSVA )
        }
        
      }else{ # if less than 30000, can directly be used as input into the GSVA function 
        topMatrixGSVA <- GSVA::gsva( as.matrix( data@assays$RNA@data) ,  geneset,  method="ssgsea", 
                              min.sz=10, max.sz=999999, abs.ranking=F, verbose=T,
                              parallel.sz=ncores)
      }
      
       
      X <- format_pathway(data, topMatrixGSVA, ncores )
      
  }
  
  if (method == "aucell"){
    
    cells_rankings <- AUCell::AUCell_buildRankings(data@assays$RNA@data, nCores = ncores, plotStats = F)
    
    
    # geneSets <- GSEABase::GeneSet(genes, setName="geneSet1") # alternative
    cells_AUC <- AUCell::AUCell_calcAUC(geneSets =  geneset,  rankings = cells_rankings, 
                                aucMaxRank=nrow(cells_rankings)*0.05,     nCores = ncores)
    
    cells_AUC <- AUCell::getAUC(cells_AUC)
    
    X <- format_pathway(data, cells_AUC, ncores)
    
  }
  
  
  return (  X )
  
}


 


# helper function to run the pathway mean 
helper_pathway_mean <- function(data, geneset , ncores = 1 ){
  
  BPparam <- generateBPParam(ncores)
  
  geneset_score_all <- BiocParallel::bplapply( geneset, function(x){
    
    exprsMat_geneset <- data[rownames(data) %in% x ,  ]
    
    geneset  <-   DelayedMatrixStats::colMeans2(DelayedArray::DelayedArray( exprsMat_geneset@assays$RNA@data )) - 
                  DelayedMatrixStats::colMeans2(DelayedArray::DelayedArray( data@assays$RNA@data ))
    geneset <- as.matrix(geneset)
    
    geneset_score <- data.frame(  geneset =   geneset  ) 
    geneset_score
  } , BPPARAM =  BPparam)
  
  geneset_score_all  <- as.data.frame(do.call(cbind,   geneset_score_all ))
  colnames( geneset_score_all ) <- names(geneset)
  geneset_score_all <- t(  geneset_score_all)
  colnames( geneset_score_all ) <- colnames(data)
  X <- format_pathway(data,   geneset_score_all, ncores)
  
  
  return ( X )
}






# helper function to calculate proportion that a gene is expressed in each cell type in a sample 
individual_geneset_proportion_celltype  <- function(data, this_geneset){
  
  # first find the average expression of the genes across cells 
  expression_level <- DelayedMatrixStats::colMeans2( DelayedArray::DelayedArray( 
    data [ rownames(data) %in% this_geneset, ]@assays$RNA@data))
  
  # find the third quantile as the threshold 
  third_quantile <-  as.numeric( quantile(  expression_level, 0.75))
  # if the expression is higher than the third quantile,
  # this gene is expressed 
  expression_level <- as.numeric(expression_level >= third_quantile)
  data$expression_high <-   expression_level
  
  # loop through all samples 
  geneset_prop <- NULL
  for (this_sample_name in unique( data$sample)){
    
    this_sample <- data@meta.data[data$sample == this_sample_name, ]
    
    # for each cell type, find the proportion of cells that this gene is expressed 
    this_geneset_prop <- table( this_sample$celltype , this_sample$expression_high)
    this_geneset_prop <- this_geneset_prop / rowSums(  this_geneset_prop)
    this_geneset_prop <- as.data.frame( this_geneset_prop)
    this_geneset_prop <-   this_geneset_prop[   this_geneset_prop$Var2 == 1, ]
    
    if (nrow( this_geneset_prop) == 0 ){
      this_geneset_prop <- data.frame(celltype =  unique(this_sample$celltype),
                                      proportion_expressed = 0 ,
                                      sample = this_sample_name)
    }else{
      this_geneset_prop <- this_geneset_prop[, c(1,3)]
      this_geneset_prop$sample <- this_sample_name
      colnames(  this_geneset_prop) <- c("celltype" , "proportion_expressed", "sample")
    }
    
    # check if any cell type in this patient is missing
    missing_celltype <-  setdiff(unique(data$celltype),  unique(this_geneset_prop$celltype))
    
    if (length(missing_celltype)> 0){
      for (this_missing_celltype in missing_celltype){
        this_geneset_prop <- rbind(this_geneset_prop, 
                                   data.frame(celltype = this_missing_celltype, proportion_expressed =  0,     sample = this_sample_name ) )
      }
    }
    this_geneset_prop <-   this_geneset_prop[ match(unique(data$celltype), this_geneset_prop$celltype) ,  ]
    
    geneset_prop  <- rbind(  geneset_prop , this_geneset_prop)
  }
  
  
  return(geneset_prop)
  
}



helper_pathway_prop <- function(data, geneset ,   ncores = 1 ){
  
  BPparam <- generateBPParam(ncores)
  
  geneset_prop_df <-  BiocParallel::bplapply( c(1:length(geneset)), function(i){
    this_geneset <- geneset[[i]]
    geneset_prop <- individual_geneset_proportion_celltype(data, this_geneset )
    geneset_prop$condition <-  data$condition[ match(geneset_prop$sample, data$sample) ]
    geneset_prop$geneset <- names(geneset)[i]
    geneset_prop
    
  }, BPPARAM =  BPparam)
  
  
  geneset_prop_df <- do.call(rbind,  geneset_prop_df)
  geneset_prop_df$celltype <- paste0( geneset_prop_df$geneset , "--", geneset_prop_df$celltype )
  geneset_prop_df <- geneset_prop_df[,  c("sample", "celltype",  "proportion_expressed")]
  geneset_prop_df <-   geneset_prop_df %>%  pivot_wider(names_from = celltype, values_from = proportion_expressed)
  
  geneset_prop_df <- as.data.frame( geneset_prop_df)
  rownames(  geneset_prop_df ) <- geneset_prop_df$sample
  geneset_prop_df  <-  geneset_prop_df [, -1]
  geneset_prop_df  <-   geneset_prop_df [ unique(data$sample), ]
  
  return (  geneset_prop_df)
}






# The output from GSVA and geneset mean is in the form of pathway score x cell 
# Need to convert to patient x pathway features 
format_pathway <- function(data, topMatrixGSVA = geneset_score_all , ncores ){
  
  # first need to aggregate the pathway score of each cell type
  topMatrixGSVA <- CreateSeuratObject(topMatrixGSVA)
  topMatrixGSVA$celltype <- data$celltype
  topMatrixGSVA$sample <- data$sample
  topMatrixGSVA$condition <- data$condition 
  bulk_data <- bulk_sample_celltype(topMatrixGSVA, ncores)
  
  # then convert to patient x (pathway a -- cell type a , pathway a -- cell type b)
  bulk_data <- reshape2:: melt(t ( as.data.frame( bulk_data@assays$RNA@data)))
  celltype <- unlist( lapply( strsplit(as.character( bulk_data$Var1 ) , "--") , `[`, 2) )
  bulk_data$Var2  <- paste0( bulk_data$Var2 , "--", celltype)
  bulk_data$Var1  <- unlist( lapply( strsplit(as.character( bulk_data$Var1 ) , "--") , `[`, 1) )
  bulk_data <- bulk_data %>% pivot_wider(names_from = Var2, values_from = value)
  bulk_data[is.na(bulk_data)] <- 0
  
  bulk_data <- as.data.frame(bulk_data)
  rownames( bulk_data) <- bulk_data$Var1
  bulk_data <-   bulk_data[, -1 ]
  
  bulk_data <- bulk_data[ unique(data$sample), ]
  
  return (  bulk_data  )
  
}







helper_pathway_mean_st  <- function( data , geneset, ncores = 1 ){
  
  
  BPparam <- generateBPParam(ncores)
  
  prob <- as.matrix(data@assays$predictions@data)
  prob <- prob[ !rownames(prob ) == "max", ]
  zero_celltype <- names( which ( rowSums( prob )  == 0 ))
  prob <-  prob[ !rownames(prob) %in% zero_celltype , ]
 
  rownames(prob) <- make.names(rownames(prob))
  
  celltype <- sort(  rownames(prob))
  
  
  x = 3
  
  allgeneset <- BiocParallel::bplapply( 1:length(geneset) ,  function(x){
    
    thisgeneset <- geneset[[x]]
    
    s <- unique(data$sample)[1]
    
    temp <-  lapply( unique(data$sample) , function(s){
      
      index <- which(data$sample == s)
      gene <- intersect( rownames(data@assays$RNA@data) ,  thisgeneset)
      
      gene_count <- as.matrix( data@assays$RNA@data[   gene, index])
      gene_name <- rownames(gene_count)
      
      thisprob <- prob[ , index]
      
      
      result_coef <- NULL
      
      i <- 1
      for ( i in c(1:nrow(gene_count))) {
        
        thisgene <- data.frame( count = gene_count[i , ]) 
        thisgene <- cbind(thisgene, t( thisprob))
        
        for ( thiscelltype in  celltype  ) {
          
          
          model <- lm( glue::glue('count ~ {thiscelltype}') , thisgene)
          a <-  summary(model)
          
          if ( nrow(a$coefficients ) == 1){
            value <- data.frame(0 )
          }else{
            value  <-  data.frame(  a$coefficients [2 , 2] )
          }
          
          rownames(   value) <- paste0( thiscelltype,  "-" , gene_name[i])
          colnames(   value ) <- "val"
          
          if (is.null( result_coef )){
            result_coef <-  data.frame( value) 
          }else{
            result_coef  <- rbind(  result_coef  , value)
          }
        }
      }
      
      result_coef
      
    }  )
    
    
    temp   <- do.call(cbind, temp  )
    colnames( temp  ) <-  unique(data$sample) 
    
    ct <-  unlist( lapply ( strsplit(rownames( temp ), "-"),`[`, 1) )
    
    final <- NULL
    i <- unique(  celltype ) [1]
    for ( i in  unique(  celltype ) ){
      thiscelltype <- which(   ct == i )
      thiscelltype <-  temp [thiscelltype , ]
      this_val <-  colSums( thiscelltype) 
      final  <-rbind(final,     this_val )
    }
    
    rownames(final) <- paste0( unique(  celltype )  ,  "-", names(geneset)[x] )
    final
    
    
  } , BPPARAM =  BPparam)
  
  allgeneset <- do.call(rbind, allgeneset)
  
  colnames(  allgeneset  ) <- unique(data$sample)
  
  allgeneset  <- t( allgeneset )
 
  return(  allgeneset)
}











