# helper: load an object shipped in CellChat (DBs, PPI) into a private env
.load_from_cellchat <- function(object_name) {
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("CellChat is not installed.", call. = FALSE)
  }
  
  e <- new.env(parent = emptyenv())
  utils::data(list = object_name, package = "CellChat", envir = e)
  if (!exists(object_name, envir = e, inherits = FALSE)) {
    stop("Failed to load '", object_name, "' from the CellChat package.", call. = FALSE)
  }
  get(object_name, envir = e, inherits = FALSE)
}

# helper: retrieve an exported function from CellChat without using CellChat::
.cellchat_fun <- function(fun_name) {
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    stop("CellChat is not installed.", call. = FALSE)
  }
  getExportedValue("CellChat", fun_name)
}

# helper function to run CCI
helper_CCI <- function(alldata, species, ncores = 1) {
  
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    cli::cli_abort(c(
      "CellChat is required to compute CCI features but is not installed.",
      "i" = "Install it from GitHub, then re-run this feature type."
    ))
  }
  
  BPparam <- generateBPParam(ncores)
  
  # Load DB + PPI once (exported to workers by BiocParallel)
  CellChatDB.human <- .load_from_cellchat("CellChatDB.human")
  CellChatDB.mouse <- .load_from_cellchat("CellChatDB.mouse")
  PPI.human        <- .load_from_cellchat("PPI.human")
  PPI.mouse        <- .load_from_cellchat("PPI.mouse")
  
  if (identical(species, "Homo sapiens")) {
    CellChatDB <- CellChatDB.human
    PPI <- PPI.human
  } else {
    CellChatDB <- CellChatDB.mouse
    PPI <- PPI.mouse
  }
  
 #  x <- unique(alldata$sample)[1]
  capture.output(
    suppressMessages(
      individual_cci <- BiocParallel::bplapply(
        unique(alldata$sample),
        function(x) {
          
          # Fetch CellChat functions *inside* the worker
          createCellChat <- .cellchat_fun("createCellChat")
          setIdent <- .cellchat_fun("setIdent")
          subsetData <- .cellchat_fun("subsetData")
          identifyOverExpressedGenes <- .cellchat_fun("identifyOverExpressedGenes")
          identifyOverExpressedInteractions <- .cellchat_fun("identifyOverExpressedInteractions")
          computeCommunProb <- .cellchat_fun("computeCommunProb")
          computeCommunProbPathway <- .cellchat_fun("computeCommunProbPathway")
          aggregateNet <- .cellchat_fun("aggregateNet")
          netVisual_bubble <- .cellchat_fun("netVisual_bubble")
          
          err <- try({
            this_sample_data <- alldata$data[, alldata$sample == x]
            colnames(this_sample_data) <- make.unique(colnames(this_sample_data))
            
            meta <- data.frame(labels = alldata$celltype[alldata$sample == x])
            rownames(meta) <- colnames(this_sample_data)
            
            cellchat <- createCellChat(
              object = this_sample_data,
              meta = meta,
              group.by = "labels"
            )
            
            cellchat <- setIdent(cellchat, ident.use = "labels")
            
            # set database / PPI
            cellchat@DB <- CellChatDB
            
            cellchat <- subsetData(cellchat)
            cellchat <- identifyOverExpressedGenes(cellchat)
            cellchat <- identifyOverExpressedInteractions(cellchat)
            cellchat <- computeCommunProb(cellchat)
            cellchat <- computeCommunProbPathway(cellchat)
            cellchat <- aggregateNet(cellchat)
            
            cellchat_score <- netVisual_bubble(cellchat, return.data = TRUE)$communication
            
            cellchat_score$feature <- paste0(
              cellchat_score$source, "->", cellchat_score$target,
              "--",
              cellchat_score$ligand, "->", cellchat_score$receptor
            )
            
            cellchat_score
          })
          
          if (inherits(err, "try-error")) {
            data.frame(
              source = "placeholder",
              target = "placeholder",
              ligand = "placeholder",
              receptor = "placeholder",
              prob = 0,
              pval = 0,
              interaction_name = "placeholder",
              interaction_name_2 = "placeholder",
              pathway_name = "placeholder",
              annotation = "placeholder",
              evidence = "placeholder",
              source.target = "placeholder",
              prob.original = 0,
              feature = "placeholder"
            )
          } else {
            err
          }
        },
        BPPARAM = BPparam
      )
    )
  )
  
  # gather sample x interaction probability matrix
  X <- NULL
  for (i in seq_along(individual_cci)) {
    temp <- individual_cci[[i]][, c("feature", "prob")]
    temp <- temp[!is.na(temp$prob), , drop = FALSE]
    temp <- temp[!duplicated(temp$feature), , drop = FALSE]
    
    if (is.null(X)) {
      X <- temp
    } else {
      X <- merge(X, temp, by = "feature", all = TRUE)
    }
  }
  
  rownames(X) <- X$feature
  X <- X[rownames(X) != "placeholder", , drop = FALSE]
  X <- X[, -1, drop = FALSE]
  colnames(X) <- unique(alldata$sample)
  X[is.na(X)] <- 0
  
  t(X)
}
