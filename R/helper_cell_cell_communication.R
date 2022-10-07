

# helper function to run CCI
helper_CCI <- function(data, species = "Homo sapiens", ncores = 1) {
  BPparam <- generateBPParam(ncores)

  # run CCI for all samples
  com_list <- BiocParallel::bplapply(unique(data$sample), function(x) {
    try({
      this_sample <- data[, data$sample == x]
      this_sample$celltype <- as.character(this_sample$celltype)
      cellchat_thissample <- cci_individual(this_sample, species)
      cellchat_thissample
    })
  }, BPPARAM = BPparam)

  names(com_list) <- unique(data$sample)

  com_list <- com_list[sapply(com_list, function(x) !inherits(x, "try-error"))]

  # gather the cell - cell interaction probability into sample x interaction probability matrix
  X <- NULL
  for (i in c(1:length(com_list))) {
    temp <- suppressMessages(CellChat::netVisual_bubble(com_list[[i]], return.data = TRUE))
    temp <- temp$communication
    temp$patient <- names(com_list)[i]

    temp$feature <- paste0(temp$source.target, "--", temp$interaction_name)
    temp <- temp[, c("feature", "patient", "prob")]

    err <- which(temp$prob == "Error in subsetCommunication_internal(net, LR, cells.level, slot.name = slot.name,  : \n  No significant signaling interactions are inferred based on the input!\n")
    if (length(err) > 0) {
      temp <- temp[-c(err), ]
    }

    temp$prob <- 10 * as.numeric(temp$prob)
    temp <- temp %>% tidyr::pivot_wider(names_from = feature, values_from = prob)
    X <- plyr::rbind.fill(X, temp)
  }

  X[is.na(X)] <- 0


  rownames(X) <- X$patient
  X <- X[, -1]
  X <- X[unique(data$sample), ]

  return(X)
}




# helper function to run CCI for each sample
cci_individual <- function(this_sample, species = "Homo sapiens") {
  data.input <- as.matrix(this_sample@assays$RNA@data)
  meta <- data.frame(
    labels = this_sample$celltype,
    row.names = names(this_sample$celltype)
  )

  cellchat <- CellChat::createCellChat(object = data.input, meta = meta, group.by = "labels")

  if (species == "Homo sapiens") {
    cellchat@DB <- CellChat::CellChatDB.human
    ppi <- CellChat::PPI.human
  } else {
    cellchat@DB <- CellChat::CellChatDB.mouse
    ppi <- CellChat::PPI.mouse
  }

  cellchat <- CellChat::subsetData(cellchat)

  # future::plan("multiprocess", workers =1)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::projectData(cellchat, ppi)

  #  future::plan("multiprocess", workers = 1)
  cellchat <- CellChat::computeCommunProb(cellchat, raw.use = TRUE)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)

  return(cellchat)
}
