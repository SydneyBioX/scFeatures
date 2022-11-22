helper_gene_mean_celltype_v2 <- function(data,
                                         genes = NULL,
                                         find_variable_genes = FALSE,
                                         num_top_gene = NULL,
                                         ncores = 1) {
    BPparam <- generateBPParam(ncores)

    lookup_table <- data@meta.data  |>
        dplyr::select(sample, celltype) |>
        tibble::rownames_to_column(var = "cell_id") |>
        dplyr::group_by(sample, celltype) |>
        dplyr::summarise(
            "n" = dplyr::n(),
            "ids" = paste(cell_id, collapse = ",")
        ) |>
        dplyr::ungroup()

    lookup_ids <- strsplit(lookup_table$ids, ",")

    matrix <- do.call("rbind", BiocParallel::bplapply(
        lookup_ids,
        function(ids) {
            # return object as numeric if only one id is present
            if (length(ids) == 1) {
                as.numeric(GetAssayData(data)[ ,ids])
            } else {
                MatrixGenerics::rowMeans2(GetAssayData(data)[ ,ids])
            }
        },
        BPPARAM = BPparam
    ))

    colnames(matrix) <- rownames(data)
    matrix <- cbind(
        sample = lookup_table$sample, celltype = lookup_table$celltype, matrix
    )
    matrix <- matrix |>
        tidyr::as_tibble() |>
        tidyr:: pivot_wider(
            id_cols = sample,
            names_from = celltype,
            values_from = -c('sample', 'celltype'),
            names_glue = "{celltype}--{.value}"
        ) |>
        as.matrix()

    rownames(matrix) <- matrix[,'sample']
    matrix <- matrix[,colnames(matrix) != 'sample']
    return(matrix)
} 

helper_gene_mean_celltype_v3 <- function(data,
                                         genes = NULL,
                                         find_variable_genes = FALSE,
                                         num_top_gene = NULL,
                                         ncores = 1) {
    BPparam <- generateBPParam(ncores)

    # an iterator function for bpiterate
    list_iter <- function(list) {
        list_iter <- iterators::iter(list)
        finished <- FALSE

        function() {
                if (finished) {
                    NULL
                } else {
                    yld <- try(iterators::nextElem(list_iter), silent = TRUE)
                    if (methods::is(yld, "try-error")) {
                        finished <<- TRUE
                        NULL
                    } else {
                        yld
                    }
                }
        }
    }

    lookup_table <- data@meta.data  |>
        dplyr::select(sample, celltype) |>
        tibble::rownames_to_column(var = "cell_id") |>
        dplyr::group_by(sample, celltype) |>
        dplyr::summarise(
            "n" = dplyr::n(),
            "ids" = paste(cell_id, collapse = ",")
        ) |>
        dplyr::ungroup()

    lookup_ids <- strsplit(lookup_table$ids, ",")

    ITER <- list_iter(lookup_ids)

    feature_matrix <- BiocParallel::bpiterate(
        ITER,
        function(ids) {
            # return object as numeric if only one id is present
            if (length(ids) == 1) {
                as.numeric(GetAssayData(data)[ ,ids])
            } else {
                MatrixGenerics::rowMeans2(GetAssayData(data)[ ,ids])
            }
        },
        BPPARAM = BPparam,
        REDUCE = rbind
    )

    colnames(feature_matrix) <- rownames(data)
    feature_matrix <- cbind(
        sample = lookup_table$sample, celltype = lookup_table$celltype, feature_matrix
    )
    feature_matrix <- feature_matrix |>
        tidyr::as_tibble() |>
        tidyr:: pivot_wider(
            id_cols = "sample",
            names_from = "celltype",
            values_from = -c('sample', 'celltype'),
            names_glue = "{celltype}--{.value}"
        ) |>
        as.matrix()

    rownames(feature_matrix) <- feature_matrix[,'sample']
    feature_matrix <- feature_matrix[,colnames(feature_matrix) != 'sample']
    return(feature_matrix)
}