devtools::load_all()
library(microbenchmark)
library(profvis)

data_small <- readRDS("refactoring/example_scrnaseq.rds")

# old implementation w/ find_var gene vs v2
res <- microbenchmark(
    cores_1 = {
        helper_gene_mean_celltype(
            data, ncores = 1
    )},
    cores_2 = {
        helper_gene_mean_celltype(
            data, ncores = 2
    )},
    cores_4 = {
        helper_gene_mean_celltype(
            data, ncores = 4
    )},
    cores_8 = {
        helper_gene_mean_celltype(
            data, ncores = 8
    )},
    cores_16 = {
        helper_gene_mean_celltype(
            data, ncores = 16
    )},
    v2_cores_1 = {
        helper_gene_mean_celltype_v2(
            data, ncores = 1
    )},
    v2_cores_2 = {
        helper_gene_mean_celltype_v2(
            data, ncores = 2
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
    times = 4
)

# old implementation on all genes vs v2
res <- microbenchmark(
    cores_1 = {
        helper_gene_mean_celltype(
            data, ncores = 1, find_variable_genes = FALSE
    )},
    cores_2 = {
        helper_gene_mean_celltype(
            data, ncores = 2, find_variable_genes = FALSE
    )},
    cores_4 = {
        helper_gene_mean_celltype(
            data, ncores = 4, find_variable_genes = FALSE
    )},
    cores_8 = {
        helper_gene_mean_celltype(
            data, ncores = 8, find_variable_genes = FALSE
    )},
    cores_16 = {
        helper_gene_mean_celltype(
            data, ncores = 16, find_variable_genes = FALSE
    )},
    v2_cores_1 = {
        helper_gene_mean_celltype_v2(
            data, ncores = 1
    )},
    v2_cores_2 = {
        helper_gene_mean_celltype_v2(
            data, ncores = 2
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
    times = 4
)


# profile old
profvis({
    helper_gene_mean_celltype(
        data, ncores = 1
    )}
)

# profile new
profvis({
    helper_gene_mean_celltype_v2(
        data, ncores = 1
    )
    }
)


res <- microbenchmark(
    cores_16 = {
        helper_gene_mean_celltype(
            data, ncores = 16, find_variable_genes = FALSE
    )},
    v2_cores_16 = {
        helper_gene_mean_celltype_v2(
            data, ncores = 16
    )},
    v3_cores_16 = {
        helper_gene_mean_celltype_v3(
            data, ncores = 16
    )},
    times = 5
)


library(refactor)
devtools::load_all()
options('refactor.time' = TRUE)

feature_gene_mean_celltype <- helper_gene_mean_celltype(
    data, ncores = 1
)