#' Remove mitochondrial and ribosomal genes, and other highly correlated genes
#'
#' @description
#' This function removes mitochondria and ribosomal genes and genes
#' highly correlated with these genes, as mitochondria and ribosomal
#' genes are typically not  interesting to look at.
#'
#' @param data A Seurat object containing expression data
#'
#' @return The Seurat object with the mitochrondrial and ribosomal genes and other highly
#' correlated genes removed
#'
#' @importFrom proxyC simil
#'
#' @examples
#'
#' data("example_scrnaseq" , package = "scFeatures")
#' data <- example_scrnaseq[, 1:20]
#' 
#' data <- remove_mito_ribo(data)
#'
#' @export
remove_mito_ribo <- function(alldata) {
  
  if ( ncol( alldata$data ) > 20000) {
    temp <- alldata$data[, sample(1:ncol(alldata$data), 20000)]
  } else {
    temp <- alldata$data
  }

  nms <- rownames(temp)

  bad_genes <- unique(
    c(
      grep("^MT-", nms, value = TRUE),
      grep("^MTMR", nms, value = TRUE),
      grep("^MTND", nms, value = TRUE),
      grep("RPL|RPS", nms, value = TRUE),
      "NEAT1", "TMSB4X", "TMSB10"
    ),
    c(
      grep("^mt-", nms, value = TRUE),
      grep("^Mtmr", nms, value = TRUE),
      grep("^Mtnd", nms, value = TRUE),
      grep("Rpl|Rps", nms, value = TRUE),
      "Neat1", "Tmsb4x", "Tmsb10"
    )
  )

  check <- sum(rownames(alldata$data) %in% bad_genes)



  if (check > 0) {
    gene_corMat <- proxyC::simil(
      temp[rownames(alldata$data) %in% bad_genes, ],
      temp[!rownames(alldata$data) %in% bad_genes, ],
      method = "correlation"
    )

    gene_corMat_max <- apply(gene_corMat, 2, max, na.rm = TRUE)
    exclude_genes <- c(bad_genes, names(gene_corMat_max)[gene_corMat_max > 0.7])

    alldata$data <- alldata$data[-c(which(rownames(alldata$data) %in% exclude_genes)), ]
  }


  return( alldata )
}




#' Identify the highly variable genes (HVGs) in the input data. By default,
#' the function calculates the HVG across the cells within each cell type,
#' as well across all cells). This is done for each sample separately, then
#' taking the union of the HVGs across all samples. The ouput is a dataframe
#' with two columns: marker and celltype. When celltype is set to FALSE, the
#' function only calculates the HVG across all cells and returns a vector of HVGs.
#' @noRd
find_var_gene <- function(alldata, num_top_gene = 1500, ncores = 1, celltype = TRUE) {
  BPparam <- generateBPParam(ncores)

  if (celltype == TRUE) {
    # here calculates the HVG across all cells across all cell types
    
    # thissample <- unique(alldata$sample)[1]
    hvg_across_all_cells <- BiocParallel::bplapply(unique(alldata$sample), function(thissample) {
      this <- alldata$data[, alldata$sample == thissample]
      gene_var <- DelayedMatrixStats::rowVars(DelayedArray::DelayedArray(this))
      top_gene <- order(gene_var, decreasing = TRUE)[1:num_top_gene]
      thisgene <- rownames(alldata$data)[top_gene]
    }, BPPARAM = BPparam)

    hvg_across_all_cells <- unique(unlist(hvg_across_all_cells))


    # below calculates the HVG within each cell type
    # thiscelltype <- unique( alldata$celltype)[1]
    gene <- BiocParallel::bplapply(unique(alldata$celltype), function(thiscelltype) {
      this_data <- alldata$data[, alldata$celltype == thiscelltype]
      this_data_sample <-  alldata$sample[ alldata$celltype == thiscelltype]
      thisgene <- c()

      #  thissample <-  unique( this_data_sample)[1]
      for (thissample in unique(this_data_sample  )) {
        this <- this_data[, this_data_sample == thissample , drop=FALSE ]
        # to calculate hvg, we want at least 10 cells
        if (ncol(this) > 10) {
          gene_var <- DelayedMatrixStats::rowVars(DelayedArray::DelayedArray(this))
          top_gene <- order(gene_var, decreasing = TRUE)[1:num_top_gene]
          temp <- rownames(alldata$data)[top_gene]
          thisgene <- c(thisgene, temp)
        }
      }

      thisgene <- unique(thisgene)

      # add the HVG  across all cells to this HVG within each cell type
      thisgene <- unique(c(thisgene, hvg_across_all_cells))

      thisgene
    }, BPPARAM = BPparam)

    names(gene) <- unique(alldata$celltype)

    all_marker <- NULL

    for (i in c(1:length(gene))) {
      try(
        {
          this <- gene[[i]]
          all_marker <- rbind(all_marker, data.frame(
            marker = this,
            celltype = names(gene)[[i]]
          ))
        },
        silent = TRUE
      )
    }

    gene <- all_marker
  } else {
    
    # thissample <- unique(alldata$sample)[1]
    gene <- BiocParallel::bplapply( unique(alldata$sample), function(thissample) {
      this <- alldata$data[, alldata$sample == thissample]
      gene_var <- DelayedMatrixStats::rowVars(DelayedArray::DelayedArray(this))
      top_gene <- order(gene_var, decreasing = TRUE)[1:num_top_gene]
      thisgene <- rownames(data)[top_gene]
    }, BPPARAM = BPparam)

    gene <- unique(unlist(gene))
  }

  return(gene)
}









#' This function is used to calculate gene expression levels in a Seurat object
#' containing expression values. By default, the function first finds the variable
#' genes per cell type using the find_var_gene function, then calculates the gene
#' expression levels for these genes in their respective cell type.
#' The output is a returns a matrix of samples by features.
#' @noRd
helper_gene_mean_celltype <- function( alldata, genes = NULL, num_top_gene = NULL, ncores = 1) {
  BPparam <- generateBPParam(ncores)


  if (is.null(num_top_gene)) {
    num_top_gene <- min(nrow( alldata$data), 100)
  }


  if (is.null(genes)) {
    all_marker <- find_var_gene(alldata ,
      num_top_gene = num_top_gene,
      ncores = ncores, celltype = TRUE
    )
  } else {
    all_marker <- genes
  }

  # j <- unique(all_marker$celltype)[1]
  final_matrix <- BiocParallel::bplapply(unique(all_marker$celltype), function(j) {
    #   i <- unique( data$sample) [1]

    gene <- all_marker[ all_marker$celltype == j, ]$marker

    # i <- unique(alldata$sample)[1]
    X_this_celltype <- lapply( unique(alldata$sample), function(i) {
      index <- intersect(
        which(alldata$celltype == j),
        which(alldata$sample == i)
      )

      if (length(index) > 0) {
        this_patient_data <- alldata$data[gene, index] 

        if (length(index) == 1) {
          this_patient_prop <- as.numeric(this_patient_data)
        } else {
          this_patient_prop <- DelayedMatrixStats::rowMeans2(DelayedArray::DelayedArray(this_patient_data))
        }
        data.frame(value = as.vector(this_patient_prop))
      } else {
        data.frame(value = as.vector(rep(0, length(gene))))
      }
    })


    X_this_celltype <- do.call(cbind, X_this_celltype)
    colnames(X_this_celltype) <- unique(alldata$sample)
    rownames(X_this_celltype) <- paste0(j, "--", gene)


    X_this_celltype
  }, BPPARAM = BPparam)

  final_matrix <- do.call(rbind, final_matrix)
  colnames(final_matrix) <- unique(alldata$sample)

  final_matrix <- t(final_matrix)
  # a <- CreateSeuratObject(  final_matrix)
  # a$sample <- unique(data$sample)
  # a$condition <- data[ , match( a$sample, data$sample)]$condition


  return(final_matrix)
}







#' This function is used to calculate the proportion of gene expression levels in a Seurat object
#' containing expression values. By default, the function first finds the variable
#' genes per cell type using the find_var_gene function, then calculates the gene
#' expression levels for these genes in their respective cell type.
#' The output is a returns a matrix of samples by features.
#' @noRd
helper_gene_prop_celltype <- function(alldata, genes = NULL, num_top_gene = NULL, ncores = 1) {
  BPparam <- generateBPParam(ncores)

  if (is.null(num_top_gene)) {
    num_top_gene <- min(nrow(alldata$data), 100)
  }

  if (is.null(genes)) {
    all_marker <- find_var_gene(alldata,
      num_top_gene = num_top_gene,
      ncores = ncores, celltype = TRUE
    )
  } else {
    all_marker <- genes
  }

  # j <- unique(data$celltype)[3]
  final_matrix <- BiocParallel::bplapply(unique(all_marker$celltype), function(j) {
    #   i <- unique( data$sample) [2]

    gene <- all_marker[all_marker$celltype == j, ]$marker

    X_this_celltype <- lapply(unique( alldata$sample), function(i) {
      index <- intersect(
        which(alldata$celltype == j),
        which(alldata$sample == i)
      )

      if (length(index) > 0) {
        this_patient_data <- alldata$data[gene, index] 
        this_patient_data <- +(this_patient_data > 1)
        if (length(index) == 1) {
          this_patient_prop <- as.numeric(this_patient_data)
        } else {
          this_patient_prop <- DelayedMatrixStats::rowMeans2(DelayedArray::DelayedArray(this_patient_data))
        }
        data.frame(value = as.vector(this_patient_prop))
      } else {
        data.frame(value = as.vector(rep(0, length(gene))))
      }
    })


    X_this_celltype <- do.call(cbind, X_this_celltype)
    colnames(X_this_celltype) <- unique(alldata$sample)
    rownames(X_this_celltype) <- paste0(j, "--", gene)


    X_this_celltype
  }, BPPARAM = BPparam)

  final_matrix <- do.call(rbind, final_matrix)

  colnames(final_matrix) <- unique(alldata$sample)

  final_matrix <- t(final_matrix)
  #
  # final_matrix <- CreateSeuratObject(  final_matrix  )
  #
  # final_matrix$sample <- unique(data$sample)
  # final_matrix$condition <- data[ , match(   final_matrix$sample, data$sample)]$condition
  #
  return(final_matrix)
}







#' This function is used to calculate the gene expression correlation in a Seurat object
#' containing expression values. By default, the function first finds the variable
#' genes per cell type using the find_var_gene function, then calculates the gene
#' expression levels for these genes in their respective cell type.
#' The output is a returns a matrix of samples by features.
#' @noRd
helper_gene_cor_celltype <- function(alldata, genes = NULL, num_top_gene = NULL, ncores = 1) {
  BPparam <- generateBPParam(ncores)

  if (is.null(num_top_gene)) {
    num_top_gene <- min(nrow(alldata$data), 5)
  }

  if (is.null(genes)) {
    all_marker <- find_var_gene(alldata,
      num_top_gene = num_top_gene,
      ncores = ncores, celltype = TRUE
    )
  } else {
    all_marker <- genes
  }

  # for each cell type, get the top 100 most variable correlation pair

  #  thiscelltype <- unique( alldata$celltype) [1]

  cor_thiscelltype <- BiocParallel::bplapply(unique(all_marker$celltype), function(thiscelltype) {
    thisdata <- alldata$data[, alldata$celltype == thiscelltype]
    thisdata_sample <- alldata$sample[alldata$celltype == thiscelltype]
    gene <- all_marker[all_marker$celltype == thiscelltype, ]$marker
    thisdata <- thisdata[gene, ]


    # x <- unique(alldata$sample)[1]

    cor_data <- lapply(unique(alldata$sample), function(x) {
      index <- which(thisdata_sample == x)

      if (length(index) > 1) {
        thispatient <- thisdata[, index]
        # cor_data <-  cor( as.matrix( t(  as.matrix( thispatient@assays$RNA@data) )))

        cor_data <- proxyC::simil(thispatient, method = "correlation")
        cor_data <- as.matrix(cor_data)

        cor_data[is.na(cor_data)] <- 0
        diag(cor_data) <- NA
        cor_data[lower.tri(cor_data)] <- NA
        cor_data <- reshape2::melt(cor_data)
        cor_data <- cor_data[!is.na(cor_data$value), ]
        temp <- data.frame(cor_data$value)
        rownames(temp) <- paste0(
          thiscelltype, "--",
          cor_data$Var1, "-with-", cor_data$Var2
        )
        temp
      } else {
        temp <- data.frame(rep(0, choose(length(gene), 2)))
        temp
      }
    })

    cor_data <- do.call(cbind, cor_data)
    colnames(cor_data) <- unique(alldata$sample)


    cor_data
  }, BPPARAM = BPparam)



  cor_thiscelltype <- do.call(rbind, cor_thiscelltype)
  colnames(cor_thiscelltype) <- unique(alldata$sample)

  cor_thiscelltype <- t(cor_thiscelltype)
  # a <- CreateSeuratObject( cor_thiscelltype)
  # a$sample <- unique(data$sample)
  # a$condition <- data[ , match( a$sample, data$sample)]$condition
  #
  return(cor_thiscelltype)
}












#' This function is used to calculate the expression of genes in each cell type
#' in a Seurat object containing expression values for spatial transcriptomics.
#' It performs a linear regression on the gene expression and cell type composition
#' at each spot to obtain regression coefficients and p-values. The regression
#' coefficients is considered as the cell type's contribution to the expression the gene.
#' Similar to the bulk deconvolution concept.
#'
#' @importFrom glue glue
#' @importFrom stats lm
#'
#' @noRd
helper_gene_mean_celltype_st <- function( alldata, genes = NULL, num_top_gene = NULL, ncores = 1) {
  BPparam <- generateBPParam(ncores)

  if (is.null(num_top_gene)) {
    num_top_gene <- min(nrow(alldata$data), 1500)
  }

  prob <- as.matrix(alldata$predictions)
  zero_celltype <- names(which(rowSums(prob) == 0))
  prob <- prob[!rownames(prob) %in% zero_celltype, ]

 
  if (is.null(genes)) {
    top_gene <- find_var_gene(alldata,
      num_top_gene = num_top_gene,
      ncores = ncores, celltype = FALSE
    )
  } else {
    top_gene <- genes
  }



  alldata$data <- alldata$data[rownames( alldata$data) %in% top_gene, ]

  rownames(prob) <- make.names(rownames(prob))

  celltype <- sort(rownames(prob))


  # s <- unique(alldata$sample)[1]

  temp <- BiocParallel::bplapply(unique(alldata$sample), function(s) {
    index <- which(alldata$sample == s)
    thispatient <- alldata$data[, index]

    gene_count <- as.matrix(alldata$data[, index])
    gene_name <- rownames(gene_count)

    thisprob <- prob[, index]

    i <- 1

    result_coef <- NULL
    result_p <- NULL

    for (i in c(1:nrow(gene_count))) {
      thisgene <- data.frame(count = gene_count[i, ])
      thisgene <- cbind(thisgene, t(thisprob))

      for (thiscelltype in celltype) {
        model <- lm(glue::glue("count ~ {thiscelltype}"), thisgene)
        a <- summary(model)

        if (nrow(a$coefficients) == 1) {
          temp <- data.frame(0)
          value <- data.frame(0)
        } else {
          #  model <-  data.frame(  coef(model) )
          temp <- data.frame(a$coefficients[2, 4])
          value <- data.frame(a$coefficients[2, 2])
        }

        rownames(temp) <- paste0(thiscelltype, "-", gene_name[i])
        colnames(temp) <- "val"

        rownames(value) <- paste0(thiscelltype, "-", gene_name[i])
        colnames(value) <- "val"


        if (is.null(result_p)) {
          result_p <- data.frame(temp)
        } else {
          result_p <- rbind(result_p, temp)
        }


        if (is.null(result_coef)) {
          result_coef <- data.frame(value)
        } else {
          result_coef <- rbind(result_coef, value)
        }
      }
    }

    list(result_p, result_coef)
  }, BPPARAM = BPparam)


  result_coef <- lapply(temp, `[[`, 2)
  result_coef <- do.call(cbind, result_coef)
  colnames(result_coef) <- unique(alldata$sample)

  result_coef <- t(result_coef)
  # final <- CreateSeuratObject(   result_coef)
  # final <- CreateSeuratObject(final)
  #
  # final$sample <- unique(data$sample)
  # final$condition <- data[ , match(  final$sample, data$sample)]$condition
  #
  return(result_coef)
}
