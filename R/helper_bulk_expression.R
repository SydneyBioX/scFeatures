
#' This function calculates the mean expression of genes in each sample in 
#' a list object. If no list of genes is provided, it will use find_var_gene 
#' function to select the top variable genes. 
#' 
#' @noRd
helper_gene_mean <- function(alldata, genes = NULL, num_top_gene = NULL, ncores = 1) {
  
    if (is.null(num_top_gene)) {
        num_top_gene <- min( nrow( alldata$data), 1500)
    }

    if (is.null(genes)) {
        gene <- find_var_gene(alldata ,
            num_top_gene = num_top_gene,
            ncores = ncores, celltype = FALSE
        )
    } else {
        gene <- genes
    }


    alldata$data <- alldata$data[gene, ]

    X <- bulk_sample(alldata, ncores)
    X <- t(X)

    return(X)
}

  
#' This function calculates the proportion of expression of genes in 
#' each sample in a list object. If no list of genes is provided, it will 
#' use find_var_gene function to select the top variable genes. 
#' @noRd
helper_gene_prop <- function( alldata, genes = NULL, num_top_gene = NULL, ncores = 1) {
  
    BPparam <- generateBPParam(ncores)

    if (is.null(num_top_gene)) {
        num_top_gene <- min(nrow(alldata$data), 1500)
    }


    if (is.null(genes)) {
        gene <- find_var_gene( alldata ,
            num_top_gene = num_top_gene,
            ncores = ncores, celltype = FALSE
        )
    } else {
        gene <- genes
    }


    alldata$data <- alldata$data[ gene , ]

    # thispatient  <- unique( alldata$sample )[1]
    gene_prop <- BiocParallel::bplapply(  unique(alldata$sample), function(thispatient) {
            this_patient_data <- alldata$data[, alldata$sample == thispatient, drop= FALSE]
            this_patient_data <- this_patient_data 
            this_patient_data <- +(this_patient_data > 1)
            this_patient_prop <- DelayedMatrixStats::rowMeans2(
                DelayedArray::DelayedArray(this_patient_data)
            )
        },
        BPPARAM = BPparam
    )

    gene_prop <- do.call(cbind, gene_prop)

    colnames(gene_prop) <- unique(alldata$sample)

    rownames(gene_prop) <- rownames(alldata$data)

    gene_prop <- t(gene_prop)

    return(gene_prop)
}







#' This function calculates the correlation of gene expression in 
#' each sample in a list object. If no list of genes is provided, it will 
#' use find_var_gene function to select the top variable genes. 
#' @noRd
helper_gene_cor <- function(alldata,  genes = NULL,   num_top_gene = NULL,   ncores = 1) {
  
    BPparam <- generateBPParam(ncores)

    if (is.null(num_top_gene)) {
        num_top_gene <- min(nrow(alldata$data), 50)
    }


    if (is.null(genes)) {
        gene <- find_var_gene(alldata,
            num_top_gene = num_top_gene,
            ncores = ncores, celltype = FALSE
        )
    } else {
        gene <- genes
    }


    alldata$data <- alldata$data[gene, ]


    #  x <- unique(alldata$sample) [1]
    cor_data <- BiocParallel::bplapply(unique(alldata$sample), function(x) {
      
        thisdata <- alldata$data[, alldata$sample == x, drop=FALSE]

        cor_data <- proxyC::simil( thisdata ,   method = "correlation"    )
        cor_data <- as.matrix(cor_data)

        cor_data[is.na(cor_data)] <- 0
        diag(cor_data) <- NA
        cor_data[lower.tri(cor_data)] <- NA
        cor_data <- reshape2::melt(cor_data)
        cor_data <- cor_data[!is.na(cor_data$value), ]
        temp <- data.frame(cor_data$value)
        rownames(temp) <- paste0(cor_data$Var1, "-with-", cor_data$Var2)
        temp
    }, BPPARAM = BPparam)

    cor_data <- do.call(cbind, cor_data)
    colnames(cor_data) <- unique(alldata$sample)


    cor_data <- t(cor_data)

    return(cor_data)
}
