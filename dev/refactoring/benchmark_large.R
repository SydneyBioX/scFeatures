library(microbenchmark)
library(SeuratObject)
library(HDF5Array)
library(refactor)
devtools::load_all()

# return normalised data with appropriate metadata
load_data <- function(filename, path = NULL, metadata = NULL) {

    if (is.null(path)) {
        path <- "/albona/nobackup2/biostat/datasets/singlecell/COVID19/withmetadata/"
    }

    all_h5_files <- list.files(
        paste0(path, "h5files"), full.names = TRUE
    )

    i <- grepl(filename, all_h5_files)
    if (!any(i)) {
       cli::cli_abort("Invalid filename: {filename} does not exist!")
    }

    print("Reading metadata")
    if (is.null(metadata)) {
       metadata <- paste0(path, "COVID_PBMC_all_cell_meta.rds")
    }
    meta <- readRDS(metadata)

    print("Reading .h5 file")
    thisfile_name <-  all_h5_files[i]
    thisfile_name  <- gsub(".h5", "", thisfile_name)

    thisfile <- HDF5Array(all_h5_files[i], "logcounts")
    thisfile_meta <- meta[match(colnames(thisfile), rownames(meta)), ]


    if (ncol(thisfile) > 800000 ){
        index <- sample(1:ncol(thisfile), 800000 )
        thisfile <- thisfile[ , index]
        thisfile_meta <- thisfile_meta[index, ]
    }
    thisfile <- CreateSeuratObject(as.matrix( thisfile ) )
    thisfile@meta.data <- thisfile_meta 

    thisfile$sample <- thisfile$meta_sample_id2
    thisfile$celltype <- thisfile$level2
    thisfile$condition <- thisfile$meta_severity

    data <- process_data(thisfile, normalise = TRUE) #  perform normalisation
    return(data)
}

# source the new function(s)
source("dev/opt_helper.R")

dataset <- "sinha"
data <- load_data(dataset)

cli::cli_inform("Running benchmarks.")

res <- microbenchmark(
    cores_1 = {
        helper_gene_mean_celltype(
            data, ncores = 1, genes = "all"
    )},
    cores_4 = {
        helper_gene_mean_celltype(
            data, ncores = 4, genes = "all"
    )},
    cores_8 = {
        helper_gene_mean_celltype(
            data, ncores = 8, genes = "all"
    )},
    cores_16 = {
        helper_gene_mean_celltype(
            data, ncores = 16, genes = "all"
    )},
    v2_cores_4 = {
        helper_gene_mean_celltype_v2(
            data, ncores = 4
    )},
    v2_cores_8 = {
        helper_gene_mean_celltype_v2(
            data, ncores = 8
    )},
    v2_cores_16 = {
        helper_gene_mean_celltype_v2(
            data, ncores = 16
    )},
    times = 10
)

saveRDS(res, paste0(dataset, "_big_bench.rds"))



res <- microbenchmark(
    cores_16 = {
        helper_gene_mean_celltype(
            data, ncores = 16, genes = "all"
    )},
    v2_cores_16 = {
        helper_gene_mean_celltype_v2(
            data, ncores = 16
    )},
    v3_cores_16 = {
        helper_gene_mean_celltype_v3(
            data, ncores = 16
    )},
    times = 10
)
