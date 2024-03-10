
#' This function retrieves the hallmark gene set from the Molecular Signature 
#' Database (MSigDB) for a specified species. The function returns a list of 
#' gene sets where each list entry contains a vector of gene symbols in the gene set.
#' @noRd
get_geneset <- function( species = "Homo sapiens") {
  
    m_df <- msigdbr::msigdbr(species = species, category = "H")
    m_t2g <- m_df |>
        dplyr::select("gs_name", "entrez_gene") |>
        as.data.frame()
    geneset <- split(x = m_t2g$entrez_gene, f = m_t2g$gs_name)

    # convert the entrez id to gene symbol
    geneset <- lapply(geneset, function(x) {
        if (species == "Homo sapiens") {
            geneID <- ensembldb::select(EnsDb.Hsapiens.v79,
                keys = as.character(x),
                keytype = "ENTREZID", columns = c("SYMBOL", "GENEID")
            )
        } else {
            geneID <- ensembldb::select(EnsDb.Mmusculus.v79,
                keys = as.character(x),
                keytype = "ENTREZID", columns = c("SYMBOL", "GENEID")
            )
        }
        this_geneset <- unique(geneID$SYMBOL)
        this_geneset
    })

    return(geneset)
}



#' This function performs gene set enrichment analysis on a list object
#' containing expression data to generate pathway score. The type of pathway 
#' analysis method currently supported are ssgsea and aucell. The function 
#' returns a matrix of samples x features with the pathway scores.
#' @noRd
helper_pathway_gsva <- function(alldata, method = "aucell", geneset, ncores = 1) {
  
    if (method == "ssgsea") {
        # if the dataset has greater than 30000 cells
        #  then it is actually too large to be computed in one go in gsva
        #  split into multiple set
        if (ncol(alldata$data) > 30000) {
            index <- seq(1, ncol(alldata$data), by = 30000)
            index <- c(index, ncol( alldata$data ))
            topMatrixGSVA <- NULL

            for (i in c(2:length(index))) {
                start <- index[i - 1]
                finish <- index[i] - 1

                message("calculating ", start, " to ", finish, " cells")
                thesecell <- as.matrix(alldata$data[, start:finish])
                
                gsvaPar <- GSVA::ssgseaParam(thesecell, geneset ,
                                             minSize = 10, maxSize  = 999999)
                temp_topMatrixGSVA <- GSVA::gsva(gsvaPar,  verbose=TRUE )
      
                topMatrixGSVA <- cbind(topMatrixGSVA, temp_topMatrixGSVA)
            }
        } else { # if <30000, can directly be used as input into the GSVA function
          
           gsvaPar <- GSVA::ssgseaParam(as.matrix(alldata$data), geneset ,
                                        minSize = 10, maxSize  = 999999)
           topMatrixGSVA <- GSVA::gsva(gsvaPar,  verbose=TRUE )
      
        }


        X <- format_pathway(alldata, topMatrixGSVA, ncores)
    }

    if (method == "aucell") {
        
      capture.output( suppressMessages( { 
        
        cells_rankings <- AUCell::AUCell_buildRankings(
            alldata$data,
            nCores = ncores, plotStats = FALSE
        )


        # geneSets <- GSEABase::GeneSet(genes, setName="geneSet1") # alternative
        cells_AUC <- AUCell::AUCell_calcAUC(
            geneSets = geneset, rankings = cells_rankings,
            aucMaxRank = nrow(cells_rankings) * 0.05, nCores = ncores
        )

        cells_AUC <- AUCell::getAUC(cells_AUC)

        X <- format_pathway(alldata, cells_AUC, ncores)
        
      }))
    }


    return(X)
}





#' The helper_pathway_mean function use the mean expression of genes in a 
#' gene set as the basis of pathway score. The input is a list object 
#' containing the gene expression data, the gene set of interest. Frist 
#' compute the mean expression of each gene in the gene set. This is 
#' done on each cell. The mean expression values are then subtracted from 
#' the overall mean expression of all genes in each cell to control for the 
#' difference in library detection. The resulting gene set scores are passed 
#' to the format_pathway function to convert them into a sample x pathway 
#' feature matrix. 
#' @noRd
helper_pathway_mean <- function(alldata, geneset, ncores = 1) {
  
    BPparam <- generateBPParam(ncores)

    geneset_score_all <- BiocParallel::bplapply(geneset, function(x) {
      
        exprsMat_geneset <- alldata$data[rownames(alldata$data) %in% x, , drop=FALSE]

        geneset <- DelayedMatrixStats::colMeans2(
            DelayedArray::DelayedArray(exprsMat_geneset)
        ) - DelayedMatrixStats::colMeans2(
            DelayedArray::DelayedArray(alldata$data)
        )
        geneset <- as.matrix(geneset)

        geneset_score <- data.frame(geneset = geneset)
        geneset_score
    }, BPPARAM = BPparam)

    geneset_score_all <- as.data.frame(do.call(cbind, geneset_score_all))
    colnames(geneset_score_all) <- names(geneset)
    geneset_score_all <- t(geneset_score_all)
    colnames(geneset_score_all) <- colnames(alldata$data)
    
    X <- format_pathway(alldata, geneset_score_all, ncores)


    return(X)
}






#' Helper function to calculate proportion that a geneset is expressed in each
#' cell type in a sample
#'
#' The function first calculates the average expression level of the genes in 
#' the gene set across all cells. It then uses the third quantile of these expression 
#' levels as a threshold to determine whether a gene is considered "expressed" in a cell.
#' It then summarise the proportion of cells expressing the gene set for each cell type.
#' This computation is done for each sample in the dataset. 
#' 
#' @param data Data to run the calculation on.
#' @param this_geneset geneset to run analysis on
#' @return a data frame containing the proportion a gene is expressed in each
#'         cell type in a sample
#'
#' @importFrom stats quantile
#' @importFrom DelayedMatrixStats colMeans2
#' @importFrom DelayedArray DelayedArray
#' 
#' @noRd
individual_geneset_proportion_celltype <- function(alldata, this_geneset) {
    # first find the average expression of the genes across cells
    expression_level <- DelayedMatrixStats::colMeans2(DelayedArray::DelayedArray(
        alldata$data[rownames(alldata$data) %in% this_geneset, , drop=FALSE] 
    ))

    # find the third quantile as the threshold
    third_quantile <- as.numeric(quantile(expression_level, 0.75, na.rm=T))
    # if the expression is higher than the third quantile,
    # this gene is expressed
    expression_level <- as.numeric(expression_level >= third_quantile)
    alldata$expression_high <- expression_level

    # loop through all samples
    geneset_prop <- NULL
    
    # this_sample_name <- unique(alldata$sample)[1]
    for (this_sample_name in unique(alldata$sample)) {
      
        this_sample_celltype <- unname( alldata$celltype[alldata$sample == this_sample_name ] )
        this_sample_expression_high  <- unname( alldata$expression_high[alldata$sample == this_sample_name])
        
        # for each cell type, find the proportion of cells expressing this geneset
        this_geneset_prop <- table(
          celltype = this_sample_celltype,
          expression_high =  this_sample_expression_high
        )
        this_geneset_prop <- this_geneset_prop / rowSums(this_geneset_prop)
        this_geneset_prop <- as.data.frame(this_geneset_prop)
        this_geneset_prop <- this_geneset_prop[this_geneset_prop$expression_high == 1, ]

        if (nrow(this_geneset_prop) == 0) {
            this_geneset_prop <- data.frame(
                celltype = unique(this_sample_celltype),
                proportion_expressed = 0,
                sample = this_sample_name
            )
        } else {
            this_geneset_prop <- this_geneset_prop[, c(1, 3)]
            this_geneset_prop$sample <- this_sample_name
            colnames(this_geneset_prop) <- c(
                "celltype", "proportion_expressed", "sample"
            )
        }

        # check if any cell type in this patient is missing
        missing_celltype <- setdiff(
            unique(alldata$celltype), unique(this_geneset_prop$celltype)
        )

        if (length(missing_celltype) > 0) {
            for (this_missing_celltype in missing_celltype) {
                this_geneset_prop <- rbind(
                    this_geneset_prop,
                    data.frame(
                        celltype = this_missing_celltype,
                        proportion_expressed = 0,
                        sample = this_sample_name
                    )
                )
            }
        }
        this_geneset_prop <- this_geneset_prop[match(
            unique(alldata$celltype),
            this_geneset_prop$celltype
        ), ]

        geneset_prop <- rbind(geneset_prop, this_geneset_prop)
    }


    return(geneset_prop)
}


#' Calculates proportion that a geneset is expressed in each cell type in a
#' sample. The function first uses the individual_geneset_proportion_celltype
#' function to calculate the proportion of cells expressing each gene set 
#' in each sample. It then concatenates the gene set name to the cell type 
#' name in the output data frame. The output data frame contains one column 
#' for each gene set-cell type combination. 
#' @noRd
helper_pathway_prop <- function(alldata, geneset, ncores = 1) {
  
    BPparam <- generateBPParam(ncores)

    geneset_prop_df <- BiocParallel::bplapply(seq_along(geneset), function(i) {
        this_geneset <- geneset[[i]]
        geneset_prop <- individual_geneset_proportion_celltype(alldata, this_geneset)
        geneset_prop$geneset <- names(geneset)[i]
        geneset_prop
    }, BPPARAM = BPparam)


    geneset_prop_df <- do.call(rbind, geneset_prop_df)
    geneset_prop_df$celltype <- paste0(
        geneset_prop_df$geneset,
        "--",
        geneset_prop_df$celltype
    )
    geneset_prop_df <- geneset_prop_df[, c(
        "sample", "celltype",
        "proportion_expressed"
    )]
    geneset_prop_df <- geneset_prop_df |>
        pivot_wider(names_from = "celltype", values_from = "proportion_expressed")

    geneset_prop_df <- as.data.frame(geneset_prop_df)
    rownames(geneset_prop_df) <- geneset_prop_df$sample
    geneset_prop_df <- geneset_prop_df[, -1]
    geneset_prop_df <- geneset_prop_df[unique(alldata$sample), ]

    return(geneset_prop_df)
}





#' This function is used to convert the pathway score output from 
#' gene set analysis into a matrix of samples x features. It takes 
#' the object containing the pathway scores (in the format of pathway score x cell), 
#' aggregates the pathway scores by cell type and converts the data into a matrix of
#' samples x features, where the features are the pathways and the cell types they 
#' are associated with. 
#' @noRd
format_pathway <- function(alldata, topMatrixGSVA, ncores) {
    # aggregate the pathway score of each cell type
  
    result_list <- list( data = topMatrixGSVA , 
                         celltype = alldata$celltype, 
                         sample = alldata$sample)
    
    
    bulk_data <- bulk_sample_celltype(result_list , ncores)

    # convert to patient x (pathway a -- cell type a , pathway a -- cell type b)
    bulk_data <- reshape2::melt(t(as.data.frame(bulk_data)))
    
    celltype <- unlist(
        lapply(strsplit(as.character(bulk_data$Var1), "--"), `[`, 2)
    )
    bulk_data$Var2 <- paste0(bulk_data$Var2, "--", celltype)
    bulk_data$Var1 <- unlist(
        lapply(strsplit(as.character(bulk_data$Var1), "--"), `[`, 1)
    )
    bulk_data <- bulk_data |>
        pivot_wider(names_from = "Var2", values_from = "value")
    bulk_data[is.na(bulk_data)] <- 0

    bulk_data <- as.data.frame(bulk_data)
    rownames(bulk_data) <- bulk_data$Var1
    bulk_data <- bulk_data[, -1]

    bulk_data <- bulk_data[unique(alldata$sample), ]

    return(bulk_data)
}





#' Thie function calculates the mean expression level of each gene set
#' in a list of gene sets for each cell type in each sample, for spatial
#' transcriptomics data. It first uses the cell type composition and 
#' gene expression at each spot and uses linear regression to calculate the 
#' regression coefficient of each cell type. It then sum the regression 
#' coefficient of each cell type and returns the resulting matrix. 
#' 
#' @importFrom glue glue
#' @importFrom stats lm
#' 
#' @noRd
helper_pathway_mean_st <- function(alldata, geneset, ncores = 1) {
    BPparam <- generateBPParam(ncores)

    prob <- as.matrix(alldata$predictions)
    zero_celltype <- names(which(rowSums(prob) == 0))
    prob <- prob[!rownames(prob) %in% zero_celltype, ]

    rownames(prob) <- make.names(rownames(prob))

    celltype <- sort(rownames(prob))


    x <- 3

    allgeneset <- BiocParallel::bplapply(seq_along(geneset), function(x) {
        thisgeneset <- geneset[[x]]

        s <- unique(alldata$sample)[1]

        temp <- lapply(unique(alldata$sample), function(s) {
            index <- which(alldata$sample == s)
            gene <- intersect(rownames(alldata$data), thisgeneset)

            gene_count <- as.matrix(alldata$data[gene, index])
            gene_name <- rownames(gene_count)

            thisprob <- prob[, index]


            result_coef <- NULL

            i <- 1
            for (i in c(seq_len(nrow(gene_count)))) {
                thisgene <- data.frame(count = gene_count[i, ])
                thisgene <- cbind(thisgene, t(thisprob))

                for (thiscelltype in celltype) {
                    model <- lm(glue::glue("count ~ {thiscelltype}"), thisgene)
                    a <- summary(model)

                    if (nrow(a$coefficients) == 1) {
                        value <- data.frame(0)
                    } else {
                        value <- data.frame(a$coefficients[2, 2])
                    }

                    rownames(value) <- paste0(thiscelltype, "-", gene_name[i])
                    colnames(value) <- "val"

                    if (is.null(result_coef)) {
                        result_coef <- data.frame(value)
                    } else {
                        result_coef <- rbind(result_coef, value)
                    }
                }
            }

            result_coef
        })


        temp <- do.call(cbind, temp)
        colnames(temp) <- unique(alldata$sample)

        ct <- unlist(lapply(strsplit(rownames(temp), "-"), `[`, 1))

        final <- NULL
        i <- unique(celltype)[1]
        for (i in unique(celltype)) {
            thiscelltype <- which(ct == i)
            thiscelltype <- temp[thiscelltype, ]
            this_val <- colSums(thiscelltype)
            final <- rbind(final, this_val)
        }

        rownames(final) <- paste0(unique(celltype), "-", names(geneset)[x])
        final
    }, BPPARAM = BPparam)

    allgeneset <- do.call(rbind, allgeneset)

    colnames(allgeneset) <- unique(alldata$sample)

    allgeneset <- t(allgeneset)

    return(allgeneset)
}
