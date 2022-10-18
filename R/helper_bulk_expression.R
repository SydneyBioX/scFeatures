

#  gene mean bulk, default to 1500 variable genes per sample
helper_gene_mean <- function(data, 
                             genes = NULL,
                             num_top_gene = NULL,
                             ncores = 1) {
    if (is.null(num_top_gene)) {
        num_top_gene <- min(nrow(data), 1500)
    }

    if (is.null(genes)) {
        gene <- find_var_gene(data,
            num_top_gene = num_top_gene,
            ncores = ncores, celltype = FALSE
        )
    } else {
        gene <- genes
    }


    data <- data[gene, ]

    X <- bulk_sample(data, ncores)
    X <- as.matrix(X@assays$RNA@data)

    X <- t(X)

    return(X)
}







# gene prop bulk
# use the proportion expression as the prediction feature
# for each variable genes, calcalate the proportion that it is expressed in each
# patient
helper_gene_prop <- function(data,
                             genes = NULL,
                             num_top_gene = NULL,
                             ncores = 1) {
    BPparam <- generateBPParam(ncores)

    if (is.null(num_top_gene)) {
        num_top_gene <- min(nrow(data), 1500)
    }


    if (is.null(genes)) {
        gene <- find_var_gene(data,
            num_top_gene = num_top_gene,
            ncores = ncores, celltype = FALSE
        )
    } else {
        gene <- genes
    }


    data <- data[gene, ]

    # thispatient  <- unique( data$sample )[1]
    gene_prop <- BiocParallel::bplapply(
        unique(data$sample), function(thispatient) {
        this_patient_data <- data[, data$sample == thispatient]
        this_patient_data <- this_patient_data@assays$RNA@data
        this_patient_data <- +(this_patient_data > 1)
        this_patient_prop <- DelayedMatrixStats::rowMeans2(
            DelayedArray::DelayedArray(this_patient_data)
        )
    }, BPPARAM = BPparam)

    gene_prop <- do.call(cbind, gene_prop)

    colnames(gene_prop) <- unique(data$sample)

    rownames(gene_prop) <- rownames(data)

    gene_prop <- t(gene_prop)

    return(gene_prop)
}







# gene correlation bulk
helper_gene_cor <- function(data,
                            genes = NULL,
                            num_top_gene = NULL,
                            ncores = 1) {

    BPparam <- generateBPParam(ncores)

    if (is.null(num_top_gene)) {
        num_top_gene <- min(nrow(data), 50)
    }


    if (is.null(genes)) {
        gene <- find_var_gene(data,
            num_top_gene = num_top_gene,
            ncores = ncores, celltype = FALSE
        )
    } else {
        gene <- genes
    }


    data <- data[gene, ]


    #  x <- unique(data$sample) [5]

    cor_data <- BiocParallel::bplapply(unique(data$sample), function(x) {
        thisdata <- data[, data$sample == x]

        cor_data <- proxyC::simil(
            thisdata@assays$RNA@data,
            method = "correlation"
        )
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
    colnames(cor_data) <- unique(data$sample)


    cor_data <- t(cor_data)

    return(cor_data)
}
