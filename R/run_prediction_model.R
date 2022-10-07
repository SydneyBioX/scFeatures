#' run cross-validated classification
#'
#' This function takes a feature matrix in the form of samples by features and performs cross-validated classification.
#'
#' @param infile Path to the input file
#'
#' @return A matrix of the infile
#'
#' @import ClassifyR
#' @importFrom dplyr %>%
#' @importFrom caret confusionMatrix
#'
#' @export
run_classification <- function(X, y, model = "randomforest", ncores = 1) {
  BPparam <- generateBPParam(ncores)

  X <- as.matrix(t(X))

  y <- as.factor(y)

  if (model == "randomforest") {
    trainParams <-  ClassifyR::TrainParams("randomForest", trees = 100)
    predictParams <- ClassifyR::PredictParams("randomForest")
  }


  if (model == "svm") {
    trainParams <- ClassifyR::TrainParams("SVM", kernel = "linear")
    predictParams <- ClassifyR::PredictParams("SVM")
  }

  if (model == "dlda") {
    trainParams <- ClassifyR::TrainParams("DLDA")
    predictParams <- ClassifyR::PredictParams("DLDA")
  }

  params <- list(trainParams, predictParams)

  result <- ClassifyR::runTests(X, y,
    datasetName = "scfeatures",
    classificationName = model,
    permutations = 20, folds = 3, seed = 2018,
    params = params, verbose = 1,
    parallelParams = BPparam
  )



  temp <- BiocParallel::bplapply(result@rankedFeatures, function(x) {
    importance <- NULL
    for (i in c(1:length(x))) {
      this <- x[[i]]
      this <- data.frame(this)
      this$importance <- 1:nrow(this)
      importance <- rbind(importance, this)
    }
    importance
  }, BPPARAM = BPparam)

  importance <- do.call(rbind, temp)

  importance <- importance %>%
    dplyr::group_by(this) %>%
    dplyr::summarize(Mean = mean(importance, na.rm = TRUE))

  importance <- importance[order(importance$Mean), ]
  colnames(importance) <- c("feature", "relative ranking")



  classLevels <- levels(ClassifyR::actualClasses(result))
  samples <- lapply(result@predictions, function(sample) factor(sample[, "sample"], levels = ClassifyR::sampleNames(result)))
  predictedClasses <- lapply(result@predictions, function(sample) factor(sample[, "class"], levels = classLevels))
  actualClasses <- lapply(result@predictions, function(sample) {
    factor(ClassifyR::actualClasses(result)[match(sample[, "sample"], ClassifyR::sampleNames(result))],
      levels = classLevels, ordered = TRUE
    )
  })


  score <- NULL
  for (i in c(1:20)) {
    predicted <- predictedClasses[[i]]
    true <- actualClasses[[i]]

    performance <- caret::confusionMatrix(predicted, reference = true)
    performance <- performance[["byClass"]]

    temp <- suppressMessages(reshape2::melt(data.frame(performance)))

    # if multi-class classification
    if (length(classLevels) > 2) {
      temp$class <- levels(predictedClasses[[1]])
    } else {
      temp$class <- names(performance)
      temp$variable <- NULL
    }
    temp$repeats <- i
    score <- rbind(score, temp)
  }




  return(list(importance = importance, performance = score))
}
