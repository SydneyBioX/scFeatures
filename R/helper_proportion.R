

#' This function calculates the proportion of cells belonging to 
#' each cell type in each sample, applicable to scRNA-seq and 
#' spatial proteomics. It also takes an optional logit argument 
#' that specifies whether the proportions should be logit-transformed.
#' @noRd
helper_proportion_raw <- function( alldata, logit = TRUE) {
  
    df <- table(alldata$sample, alldata$celltype)
    df <- df / rowSums(df)

    # if logit transformation is needed,
    # need to do do the following to avoid infinite or NA
    if (logit) {
        try(df[df == 0] <- 0.001, silent = TRUE)
        try(df[df == 1] <- 0.999, silent = TRUE)
        df <- gtools::logit(df)
    }


    df <- as.data.frame(df)

    colnames(df) <- c("sample", "celltype", "proportion")


    df <- df |>
        tidyr::pivot_wider(names_from = "celltype", values_from = "proportion")
    df <- as.data.frame(df)
    rownames(df) <- df$sample
    df <- df[, -1]

    df <- df[unique(alldata$sample), ]

    return(df)
}



#' This function calculates the ratio of cell type proportion between 
#' two cell types in each sample, applicable to scRNA-seq and spatial
#' proteomics. The ratio is log2 transformed. 
#' @noRd
helper_proportion_ratio <- function( alldata, ncores = 1) {
  
    BPparam <- generateBPParam(ncores)

    allcelltype <- unique( alldata$celltype)

    #  x =  unique( alldata$sample)
    df <- BiocParallel::bplapply(unique( alldata$sample), function(x) {
        # loop through each sample
        this_sample_celltype <- alldata$celltype[alldata$sample == x ]

        # keep track of ratio for this sample
        temp_df <- NULL
        for (i in (seq_along(allcelltype))) {
            for (j in (seq_along(allcelltype))) {
                if (i < j) {
                    celltype1 <- allcelltype[i]
                    celltype2 <- allcelltype[j]

                    celltype_1 <- sum(this_sample_celltype %in% celltype1)
                    celltype_2 <- sum(this_sample_celltype %in% celltype2)


                    # if one of the cell type is missing,  add 1
                    # this is to avoid issues such as division by zero
                    if (celltype_1 == 0 || celltype_2 == 0) {
                        celltype_1 <- celltype_1 + 1
                        celltype_2 <- celltype_2 + 2
                    }

                    ratio <- celltype_1 / celltype_2

                    temp <- data.frame(
                        sample = x,
                        celltype = paste0(celltype1, "-vs-", celltype2),
                        ratio = ratio
                    )

                    temp_df <- rbind(temp_df, temp)
                }
            }
        }
        temp_df
    }, BPPARAM = BPparam)

    df <- as.data.frame(do.call(rbind, df))

    df$ratio <- log2(df$ratio + 1)


    df <- as.data.frame(df)

    colnames(df) <- c("sample", "celltype", "proportion")


    df <- df |>
        tidyr::pivot_wider(names_from = "celltype", values_from = "proportion")
    df <- as.data.frame(df)
    rownames(df) <- df$sample
    df <- df[, -1, drop = FALSE]

    df <- df[unique(alldata$sample), , drop = FALSE]


    return(df)
}














#' This function calculates the proportion of cells belonging to 
#' each cell type in each sample, applicable to spatial transcriptomics. 
#' It takes an optional logit argument that specifies whether the proportions 
#' should be logit-transformed.
#' @noRd
#' 
helper_proportion_raw_st <- function(alldata, logit = TRUE, ncores = 1) {
    BPparam <- generateBPParam(ncores)

 
    num_cell_spot <- get_num_cell_per_celltype(alldata)

    prop_table <- BiocParallel::bplapply(unique(alldata$sample), function(s) {
        index <- which(alldata$sample == s)
        celltype <- num_cell_spot[, index]
        celltype <- rowSums(celltype) / sum(celltype)
        data.frame(sample = s, celltype = names(celltype), proportion = celltype)
    }, BPPARAM = BPparam)

    prop_table <- do.call(rbind, prop_table)

    colnames(prop_table) <- c("sample", "celltype", "proportion")

    tab <- prop_table

    if (logit) {
        try(tab[tab$proportion == 0, ]$proportion <- 0.001, silent = TRUE)
        try(tab[tab$proportion == 1, ]$proportion <- 0.999, silent = TRUE)
        tab$proportion <- gtools::logit(tab$proportion) # this computes the logit transformation
    }

    tab <- as.data.frame(tab)

    colnames(tab) <- c("sample", "celltype", "proportion")


    tab <- tab |>
        tidyr::pivot_wider(names_from = "celltype", values_from = "proportion")
    tab <- as.data.frame(tab)
    rownames(tab) <- tab$sample
    tab <- tab[, -1]

    tab <- tab[unique(alldata$sample), ]

    return(tab)
}




#' This function calculates the ratio of cell type proportion between 
#' two cell types in each sample, applicable to spatial transcriptomics. 
#' The ratio is log2 transformed. 
#' @noRd
helper_proportion_ratio_st <- function(alldata, ncores = 1) {
    BPparam <- generateBPParam(ncores)

    num_cell_spot <- get_num_cell_per_celltype(alldata)
    allcelltype <- rownames(num_cell_spot)

    # x <- unique(alldata$sample)[1]
    df <- BiocParallel::bplapply(unique(alldata$sample), function(x) {
        # loop through each sample

        index <- which(alldata$sample == x)
        ct <- num_cell_spot[, index]
        ct <- rowSums(ct)


        # keep track of ratio for this sample
        temp_df <- NULL
        for (i in (seq_along(allcelltype))) {
            for (j in (seq_along(allcelltype))) {
                if (i < j) {
                    celltype1 <- allcelltype[i]
                    celltype2 <- allcelltype[j]

                    celltype_1 <- sum(ct[celltype1])
                    celltype_2 <- sum(ct[celltype2])


                    # if one of the cell type is missing,  add 1
                    # this is to avoid issues such as division by zero
                    if (celltype_1 == 0 || celltype_2 == 0) {
                        celltype_1 <- celltype_1 + 1
                        celltype_2 <- celltype_2 + 2
                    }

                    ratio <- celltype_1 / celltype_2

                    temp <- data.frame(
                        sample = x,
                        celltype = paste0(celltype1, "-vs-", celltype2),
                        ratio = ratio
                    )

                    temp_df <- rbind(temp_df, temp)
                }
            }
        }
        temp_df
    }, BPPARAM = BPparam)

    tab <- as.data.frame(do.call(rbind, df))

    tab$ratio <- log2(tab$ratio + 1)

    colnames(tab) <- c("sample", "celltype", "ratio")

    tab <- tab |> tidyr::pivot_wider(names_from = "celltype", values_from = "ratio")
    tab <- as.data.frame(tab)
    rownames(tab) <- tab$sample
    tab <- tab[, -1]

    tab <- tab[unique(alldata$sample), ]


    return(tab)
}
