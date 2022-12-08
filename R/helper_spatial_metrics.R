

#' This function calculates the pairwise distance between cell types
#' for a sample by using the coordinates and cell types of the cells. 
#' It keeps the cell types with more than 10 cells and calculates
#' the pairwise distance between these cell types. The function returns 
#' a list of the pairwise distances between the cell types.
#' This function is designed for spatial proteomics. 
individual_celltype_interaction_sp <- function(this_sample) {
    cell_points <- spatstat.geom::ppp(
        x = this_sample$x_cord,
        y = this_sample$y_cord,
        check = FALSE,
        yrange = c(0, max(this_sample$y_cord)),
        xrange = c(0, max(this_sample$x_cord)),
        marks = as.factor(this_sample$celltype)
    )

    tab <- table(this_sample$celltype)

    cellTypes_toTest <- names(tab[which(tab > 10)])
    cellTypes_pair <- expand.grid(cellTypes_toTest, cellTypes_toTest,
        stringsAsFactors = FALSE
    )


    # Calcualte the pairwise distance
    d <- spatstat.geom::pairdist(cell_points, squared = FALSE)
    diag(d) <- Inf

    nn_list <- apply(d, 1, function(x) which(x < 50))

    nn_list_cellTypes <- lapply(seq_along(nn_list), function(idx) {
        if (length(nn_list[[idx]]) > 0) {
            paste(this_sample$celltype[idx],
                this_sample$celltype[nn_list[[idx]]],
                sep = "_"
            )
        }
    })

    nn_list_cellTypes <- unlist(nn_list_cellTypes)
    nn_list_cellTypes <- rearrange_string(nn_list_cellTypes)
    nn_list_cellTypes <- table(nn_list_cellTypes)

    return(nn_list_cellTypes)
}




#' Calculates the pairwise distance between cell types for spatial proteomics. 
#' It applies the individual_celltype_interaction_sp function to each sample,
#' then merges the results from each sample. 
helper_celltype_interaction_sp <- function(data, ncores = 1) {
    BPparam <- generateBPParam(ncores)

    s <- unique(data$sample)[1]

    nn_list_cellTypes <- BiocParallel::bplapply(unique(data$sample), function(s) {
        this_sample <- data[, data$sample == s]
        nn_list_cellTypes <- individual_celltype_interaction_sp(this_sample)
    }, BPPARAM = BPparam)


    temp <- NULL

    for (i in seq_along(nn_list_cellTypes)) {
        err <- try(
            {
                a <- nn_list_cellTypes[[i]]
                a <- data.frame(a)
                a$Freq <- a$Freq / sum(a$Freq)

                if (is.null(temp)) {
                    temp <- a
                } else {
                    temp <- suppressWarnings(
                        merge(temp, a, by = "nn_list_cellTypes", all = TRUE)
                    )
                }
            },
            silent = TRUE
        )


        if (is(err, "try-error")) {
            a <- data.frame(rep(0, nrow(temp)))
            a$nn_list_cellTypes <- temp$nn_list_cellTypes
            temp <- suppressWarnings(
                merge(temp, a, by = "nn_list_cellTypes", all = TRUE)
            )
        }

        colnames(temp) <- make.names(colnames(temp), unique = TRUE)
    }


    rownames(temp) <- temp$nn_list_cellTypes
    temp <- temp[, -1]

    colnames(temp) <- unique(data$sample)

    temp <- t(temp)

    temp[is.na(temp)] <- 0
    nn_list_cellTypes <- temp


    return(nn_list_cellTypes)
}





#' This function calculates the cell-type interactions for spatial 
#' transcriptomic. It takes probabilities of cell-type assignments 
#' for each spot and perform matrix multiplication to get the probability of 
#' cell-type co-occurrences for each spot. The output is a  vector containing 
#' the co-occurrence probabilities for each cell-type pair.
individual_celltype_interaction_st <- function(thisprob) {
    x <- 1
    temp <- lapply(seq_len(ncol(thisprob)), function(x) {
        thisspot <- thisprob[, x]
        thisspot <- thisspot %*% t(thisspot)
        rownames(thisspot) <- colnames(thisspot)
        thisspot <- reshape2::melt(thisspot)
        a <- data.frame(thisspot$value)
        rownames(a) <- paste0(thisspot$Var1, "-with-", thisspot$Var2)
        a
    })

    temp <- do.call(cbind, temp)
    nn_list_cellTypes <- rowSums(temp)


    return(nn_list_cellTypes)
}



#' Calculates the cell-type interactions for spatial transcriptomic.
#' It applies the individual_celltype_interaction_st function to each sample,
#' then merges the results from each sample. 
helper_celltype_interaction_st <- function(data, ncores = 1) {
    BPparam <- generateBPParam(ncores)

    s <- unique(data$sample)[1]

    prob <- data@assays$predictions
    prob <- as.matrix(prob@data)
    prob <- prob[!rownames(prob) == "max", ]
    zero_celltype <- names(which(rowSums(prob) == 0))
    prob <- prob[!rownames(prob) %in% zero_celltype, ]


    nn_list_cellTypes <- BiocParallel::bplapply(unique(data$sample), function(s) {
        index <- which(data$sample == s)
        thisprob <- prob[, index]
        nn_list_cellTypes <- individual_celltype_interaction_st(thisprob)
    }, BPPARAM = BPparam)


    nn_list_cellTypes <- do.call(cbind, nn_list_cellTypes)
    colnames(nn_list_cellTypes) <- unique(data$sample)

    nn_list_cellTypes <- t(nn_list_cellTypes)

    return(nn_list_cellTypes)
}










#' This function calculates the L-statistics for a given sample for spatial
#' transcriptomics data using the spatial coordinates, cell type information and 
#' the number of cells in each spatial location. The L-statistic measures the 
#' degree of spatial clustering between two objects in a particular area. The output 
#' vector represents the L-statistic for a pair of cell types.
individual_L_stat_st <- function(thissample, this_num_cell_per_spot) {
    # expand each spot into its number of cells
    x <- c()
    y <- c()

    celltype <- c()

    i <- 1


    gap_x <- (max(x_c <- thissample$x_cord) - min(x_c)) / length(x_c) / 2
    gap_y <- (max(y_c <- thissample$y_cord) - min(y_c)) / length(y_c) / 2


    for (i in seq_len(ncol(thissample))) {
        thisspot <- thissample[, i]
        thisspot_num_cell <- this_num_cell_per_spot[, i]

        total_num_cell <- sum(thisspot_num_cell)

        this_gap_x <- gap_x / total_num_cell
        this_gap_y <- gap_y / total_num_cell

        x <- c(
            x,
            seq(thissample$x_cord[i], by = this_gap_x, length.out = total_num_cell)
        )
        y <- c(
            y,
            seq(thissample$y_cord[i], by = this_gap_y, length.out = total_num_cell)
        )

        celltype <- c(
            celltype, rep(rownames(this_num_cell_per_spot), thisspot_num_cell)
        )
    }


    cell_points_threecelltype <- spatstat.geom::ppp(
        x = x,
        y = y,
        check = FALSE,
        yrange = c(
            min(as.numeric(thissample$y_cord)),
            max(as.numeric(thissample$y_cord))
        ),
        xrange = c(
            min(as.numeric(thissample$x_cord)),
            max(as.numeric(thissample$x_cord))
        ),
        marks = as.factor(celltype)
    )


    cellTypes_toTest <- names(which(table(celltype) > 10))
    cellTypes_pair <- expand.grid(cellTypes_toTest, cellTypes_toTest,
        stringsAsFactors = FALSE
    )



    L_patient <- lapply(seq_len(nrow(cellTypes_pair)), function(i) {
        L_stats(cell_points_threecelltype,
            from = cellTypes_pair[i, 1],
            to = cellTypes_pair[i, 2],
            L_dist = 4
        )
    })

    L_patient <- do.call(c, L_patient)
    names(L_patient) <- paste(cellTypes_pair[, 1], cellTypes_pair[, 2], sep = "_")

    return(L_patient)
}




helper_L_stat_st <- function(data, ncores = 1) {
    BPparam <- generateBPParam(ncores)

    num_cell_per_spot <- get_num_cell_per_celltype(data)

    s <- unique(data$sample)[1]


    L_stats <- BiocParallel::bplapply(unique(data$sample), function(s) {
        index <- which(data$sample == s)

        thissample <- data[, index]
        this_num_cell_per_spot <- num_cell_per_spot[, index]

        L_patient <- individual_L_stat_st(thissample, this_num_cell_per_spot)
    }, BPPARAM = BPparam)


    temp <- NULL

    for (i in seq_along(L_stats)) {
        err <- try(
            {
                a <- L_stats[[i]]
                a <- data.frame(a)
                a$rowname <- rownames(a)

                if (is.null(temp)) {
                    temp <- a
                } else {
                    temp <- suppressWarnings(merge(temp, a, by = "rowname", all = TRUE))
                }
            },
            silent = TRUE
        )


        if (is(err, "try-error")) {
            a <- data.frame(rep(0, nrow(temp)))
            a$rowname <- temp$rowname
            temp <- suppressWarnings(merge(temp, a, by = "rowname", all = TRUE))
        }

        colnames(temp) <- make.names(colnames(temp), unique = TRUE)
    }


    rownames(temp) <- temp$rowname
    temp <- temp[, -1]

    colnames(temp) <- unique(data$sample)

    temp <- t(temp)

    temp[is.na(temp)] <- 0
    L_patient <- temp

    return(L_patient)
}




individual_L_stat_sp <- function(this_sample) {
    cell_points <- spatstat.geom::ppp(
        x = this_sample$x_cord,
        y = this_sample$y_cord,
        check = FALSE,
        yrange = c(0, max(this_sample$y_cord)),
        xrange = c(0, max(this_sample$x_cord)),
        marks = as.factor(this_sample$celltype)
    )

    tab <- table(this_sample$celltype)

    cellTypes_toTest <- names(tab[which(tab > 10)])
    cellTypes_pair <- expand.grid(cellTypes_toTest, cellTypes_toTest,
        stringsAsFactors = FALSE
    )



    L_patient <- list()
    for (i in seq_len(nrow(cellTypes_pair))) {
        L_patient[[i]] <- L_stats(cell_points,
            from = cellTypes_pair[i, 1],
            to = cellTypes_pair[i, 2],
            L_dist = 50
        )
    }

    L_patient <- do.call(c, L_patient)
    names(L_patient) <- paste(cellTypes_pair[, 1], cellTypes_pair[, 2], sep = "_")

    return(L_patient)
}



helper_L_stat_sp <- function(data, ncores = 1) {
    BPparam <- generateBPParam(ncores)

    s <- unique(data$sample)[1]


    L_stats_result <- BiocParallel::bplapply(unique(data$sample), function(s) {
        index <- which(data$sample == s)

        thissample <- data[, index]

        L_patient <- individual_L_stat_sp(thissample)
    }, BPPARAM = BPparam)


    temp <- NULL

    for (i in seq_along(L_stats_result)) {
        err <- try(
            {
                a <- L_stats_result[[i]]
                a <- data.frame(a)
                a$rowname <- rownames(a)

                if (is.null(temp)) {
                    temp <- a
                } else {
                    temp <- suppressWarnings(merge(temp, a, by = "rowname", all = TRUE))
                }
            },
            silent = TRUE
        )


        if (is(err, "try-error")) {
            a <- data.frame(rep(0, nrow(temp)))
            a$rowname <- temp$rowname
            temp <- suppressWarnings(merge(temp, a, by = "rowname", all = TRUE))
        }

        colnames(temp) <- make.names(colnames(temp), unique = TRUE)
    }


    rownames(temp) <- temp$rowname
    temp <- temp[, -1]

    colnames(temp) <- unique(data$sample)

    temp <- t(temp)

    temp[is.na(temp)] <- 0
    L_patient <- temp

    return(L_patient)
}









#' Calculates the nearest neighbor correlation for a given sample 
#' using the spatial coordinates and expression values of cells. 
#' It returns a vector of correlation values for each gene. 
#' This function is used in as a helper function in helper_nncorr_protein 
#' function to generate the nearest neighbor correlation for all samples.
individual_nncorr_protein <- function(thissample) {
    exprsMat <- thissample@assays$RNA@data

    cell_points_cts <- spatstat.geom::ppp(
        x = as.numeric(thissample$x_cord), y = as.numeric(thissample$y_cord),
        check = FALSE,
        xrange = c(
            min(as.numeric(thissample$x_cord)),
            max(as.numeric(thissample$x_cord))
        ),
        yrange = c(
            min(as.numeric(thissample$y_cord)),
            max(as.numeric(thissample$y_cord))
        ),
        marks = t(as.matrix(exprsMat))
    )
    err <- try(
        {
            nncorr_protein <- spatstat.explore::nncorr(cell_points_cts)["correlation", ]
        },
        silent = TRUE
    )
    if (is(err, "try-error")) {
        nncorr_protein <- NA
    }

    return(nncorr_protein)
}



#' Calculates the nearest neighbor correlation all samples 
#' using the spatial coordinates and expression values of cells. 
#' It applies the individual_nncorr_protein function to each sample,
#' then merges the results from each sample. 
#' @importFrom methods is
helper_nncorr_protein <- function(data, num_top_gene = NULL, ncores = 1) {
    BPparam <- generateBPParam(ncores)

    if (is.null(num_top_gene)) {
        num_top_gene <- min(nrow(data), 1500)
    }

    top_gene <- find_var_gene(data,
        num_top_gene = num_top_gene,
        ncores = ncores, celltype = FALSE
    )

    data@assays$RNA@data <- data@assays$RNA@data[rownames(data@assays$RNA@data) %in% top_gene, ]


    s <- unique(data$sample)[1]

    nncorr_protein <- BiocParallel::bplapply(unique(data$sample), function(s) {
        thissample <- data[, data$sample == s]
        L_patient <- individual_nncorr_protein(thissample)
    }, BPPARAM = BPparam)


    temp <- NULL

    for (i in seq_along(nncorr_protein)) {
        err <- try(
            {
                a <- nncorr_protein[[i]]
                a <- data.frame(a)
                a$rowname <- rownames(a)

                if (is.null(temp)) {
                    temp <- a
                } else {
                    temp <- suppressWarnings(merge(temp, a, by = "rowname", all = TRUE))
                }
            },
            silent = TRUE
        )


        if (is(err, "try-error")) {
            a <- data.frame(rep(0, nrow(temp)))
            a$rowname <- temp$rowname
            temp <- suppressWarnings(merge(temp, a, by = "rowname", all = TRUE))
        }

        colnames(temp) <- make.names(colnames(temp), unique = TRUE)
    }


    rownames(temp) <- temp$rowname
    temp <- temp[, -1]

    colnames(temp) <- unique(data$sample)

    temp <- t(temp)

    temp[is.na(temp)] <- 0
    nncorr_protein <- temp


    return(nncorr_protein)
}




#' This function computes the Moran's I for each gene in a sample.
#' Moran's I is a measure of spatial autocorrelation, which indicates 
#' how strongly the gene expression values in a sample cluster or disperse. 
individual_moran_cor <- function(thissample) {
    exprsMat <- thissample@assays$RNA@data

    cell_points_cts <- spatstat.geom::ppp(
        x = as.numeric(thissample$x_cord), y = as.numeric(thissample$y_cord),
        check = FALSE,
        xrange = c(
            min(as.numeric(thissample$x_cord)),
            max(as.numeric(thissample$x_cord))
        ),
        yrange = c(
            min(as.numeric(thissample$y_cord)),
            max(as.numeric(thissample$y_cord))
        ),
        marks = t(as.matrix(exprsMat))
    )

    d <- spatstat.geom::pairdist(cell_points_cts, squared = FALSE)
    diag(d) <- Inf


    w <- 1 / d


    moran_cor <- lapply(seq_len(nrow(exprsMat)), function(x) {
        err <- try(val <- ape::Moran.I(exprsMat[x, ], w)$observed, silent = TRUE)
        if (is(err, "try-error")) {
            NA
        } else {
            val
        }
    })

    names(moran_cor) <- rownames(exprsMat)
    moran_cor <- unlist(moran_cor)

    return(moran_cor)
}




#' Calculates the Moran's I value for all samples. 
#' It applies the individual_moran_cor function to each sample,
#' then merges the results from each sample. 
helper_moran <- function(data, num_top_gene = NULL, ncores = 1) {
    BPparam <- generateBPParam(ncores)

    if (is.null(num_top_gene)) {
        num_top_gene <- min(nrow(data), 1500)
    }

    top_gene <- find_var_gene(data,
        num_top_gene = num_top_gene,
        ncores = ncores, celltype = FALSE
    )

    data@assays$RNA@data <- data@assays$RNA@data[
        rownames(data@assays$RNA@data) %in% top_gene,
    ]


    s <- unique(data$sample)[1]

    moran_cor <- BiocParallel::bplapply(unique(data$sample), function(s) {
        thissample <- data[, data$sample == s]
        moran_cor <- individual_moran_cor(thissample)
    }, BPPARAM = BPparam)

    temp <- NULL

    for (i in seq_along(moran_cor)) {
        err <- try(
            {
                a <- moran_cor[[i]]
                a <- data.frame(a)
                a$rowname <- rownames(a)

                if (is.null(temp)) {
                    temp <- a
                } else {
                    temp <- suppressWarnings(merge(temp, a, by = "rowname", all = TRUE))
                }
            },
            silent = TRUE
        )


        if (is(err, "try-error")) {
            a <- data.frame(rep(0, nrow(temp)))
            a$rowname <- temp$rowname
            temp <- suppressWarnings(merge(temp, a, by = "rowname", all = TRUE))
        }

        colnames(temp) <- make.names(colnames(temp), unique = TRUE)
    }


    rownames(temp) <- temp$rowname
    temp <- temp[, -1]

    colnames(temp) <- unique(data$sample)

    temp <- t(temp)

    temp[is.na(temp)] <- 0
    moran_cor <- temp


    return(moran_cor)
}
