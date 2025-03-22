










#'Functional_classification
#'
#' @param folds The number of folds to use for cross - validation. Default is 5, which means the data will be split into 5 subsets for cross - validation purposes.
#' @param data2 Matrix of preprocessed single-cell metabolomics data
#' @param label_col A vector indicating the group labels (e.g., 0 and 1) for each sample in the data, used to distinguish different groups for the statistical test.
#' @param method The classification method to be used. Valid options include "AdaBoost", "Bagging", "Decision Trees", etc.
#'
#' @return A list containing two elements. The first element is a data frame with the true and predicted labels for each fold of the cross - validation. The second element is a data frame with the top 10 most important features (as determined by the selected classification method) for each fold.
#' @export
#' @examples
#'\donttest{
#' }
#'
#'
#'
#'
#'
#'
classification <- function(data2, label_col, folds = 5, method) {

  library(adabag)
  library(C50)
  library(pROC)
  library(kknn)
  library(MASS)
  library(e1071)
  library(AUC)
  library(multiROC)
  library(caret)
  library(mlbench)
  library(dummies)
  library(randomForest)


  labels <- as.factor(data2[, label_col])
  data <- apply(data2[, -c(1:4)], 2, as.numeric)
  data4 <- data2[, 1:4]

  data[data == 0] <- NA
  data <- as.data.frame(data)


  if (method == "AdaBoost") {
    set.seed(10)
    final_res <- list()
    data_label <- as.data.frame(cbind(labels, data))
    colnames(data_label)[1] <- "labels"
    test.fold <- split(sample(1:length(labels)), 1:folds)

    for (mm in 1:folds) {
      test <- test.fold[[mm]]
      train_df <- data_label[-test, ]
      test_df <- data_label[test, ]
      train_df$labels <- as.factor(train_df$labels)
      test_df$labels <- as.factor(test_df$labels)

      boost_model <- boosting(labels ~ ., data = train_df)
      boost_pred <- predict(boost_model, test_df)

      true_label <- dummies::dummy(test_df$labels, sep = ".")
      true_label <- data.frame(true_label)
      colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
      colnames(true_label) <- paste0(colnames(true_label), "_true")

      pred_label <- boost_pred$prob
      colnames(pred_label) <- gsub("true", "pred_SVM", colnames(true_label))

      final_df_2 <- cbind(true_label, pred_label)
      final_df2_t <- t(final_df_2)
      final_df2_t_la <- rownames(final_df2_t)
      final_df2_t_c <- as.data.frame(cbind(final_df2_t_la, final_df2_t))
      colnames(final_df2_t_c)[-1] <- paste0("Fold_", mm, "_", colnames(final_df2_t_c)[-1])

      import_f <- boost_model$importance
      import_s <- import_f[order(import_f, decreasing = T)][1:10]
      import_a <- as.data.frame(cbind(row.names(t(t(import_s))), t(t(import_s))))
      colnames(import_a) <- c("Feature", paste0("Importance_Fold_", mm))

      if (mm == 1) {
        final_df_m <- final_df2_t_c
        imports <- import_a
      } else {
        final_df_m <- merge(final_df_m, final_df2_t_c, by = "final_df2_t_la", all = TRUE)
        imports <- merge(imports, import_a, by = "Feature", all = TRUE)
      }
    }
    final_res[[1]] <- final_df_m
    final_res[[2]] <- imports
    return(final_res)
  } else if (method == "Bagging") {
    set.seed(10)
    final_res <- list()
    data_label <- cbind(labels, data)
    x <- data_label
    y <- factor(data_label[,1])
    colnames(x)[1] <- "labels"
    test.fold <- split(sample(1:length(labels)), 1:folds)

    for (mm in 1:folds) {
      test <- test.fold[[mm]]
      train_df <- x[-test, ]
      test_df <- x[test, ]

      boost_model <- bagging(labels~., data=train_df)
      boost_pred <- predict.bagging(boost_model, test_df)

      true_label <- dummies::dummy(test_df$labels, sep = ".")
      true_label <- data.frame(true_label)
      colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
      colnames(true_label) <- paste0(colnames(true_label), "_true")

      pred_label <- boost_pred$prob
      colnames(pred_label) <- gsub("true","pred_SVM",colnames(true_label))

      final_df_2 <- cbind(true_label, pred_label)
      final_df2_t <- t(final_df_2)
      final_df2_t_la <- rownames(final_df2_t)
      final_df2_t_c <- as.data.frame(cbind(final_df2_t_la, final_df2_t))

      import_f <- boost_model$importance
      import_s <- import_f[order(import_f, decreasing = T)][1:10]
      import_a <- as.data.frame(cbind(row.names(t(t(import_s))), t(t(import_s))))
      colnames(import_a) <- c("V1", paste0("V_", mm))

      if (mm == 1) {
        final_df_m <- final_df2_t_c
        imports <- import_a
      } else {
        final_df_m <- merge(final_df_m, final_df2_t_c, by = "final_df2_t_la", all = T)
        imports <- merge(imports, import_a, by = "V1", all = T)
      }
    }
    final_res[[1]] <- final_df_m
    final_res[[2]] <- imports
    return(final_res)
  } else if (method == "Decision Trees") {
    set.seed(10)
    final_res <- list()
    data_label <- cbind(labels, data)
    x <- data_label
    y <- factor(data_label[,1])
    colnames(x)[1] <- "labels"
    test.fold <- split(sample(1:length(labels)), 1:folds)

    for (mm in 1:folds) {
      test <- test.fold[[mm]]
      train_df <- x[-test, ]
      test_df <- x[test, ]

      credit_boost <- C5.0(labels~., data = train_df, trials = 10, probability = TRUE)
      credit_boost_pred <- predict(credit_boost, test_df, probability = TRUE)

      true_label <- dummies::dummy(test_df$labels, sep = ".")
      true_label <- data.frame(true_label)
      colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
      colnames(true_label) <- paste0(colnames(true_label), "_true")

      pred_label <- dummies::dummy(credit_boost_pred, sep = ".")
      pred_label <- data.frame(pred_label)
      colnames(pred_label) <- gsub(".*?\\.", "", colnames(pred_label))
      colnames(pred_label) <- paste0(colnames(pred_label), "_pred_SVM")

      final_df_2 <- cbind(true_label, pred_label)
      final_df2_t <- t(final_df_2)
      final_df2_t_la <- rownames(final_df2_t)
      final_df2_t_c <- as.data.frame(cbind(final_df2_t_la,final_df2_t))

      import_f <- as.data.frame(varImp(credit_boost))
      import_fn <- import_f[,1]
      import_s <- import_fn[order(import_fn, decreasing = T)][1:10]
      import_a <- as.data.frame(cbind(row.names(t(t(import_f))),t(t(import_s))))
      colnames(import_a) <- c("V1", paste0("V_", mm))

      if (mm == 1) {
        final_df_m <- final_df2_t_c
        imports <- import_a
      } else {
        final_df_m <- merge(final_df_m, final_df2_t_c, by = "final_df2_t_la", all = T)
        imports <- merge(imports, import_a, by = "V1", all = T)
      }
    }
    final_res[[1]] <- final_df_m
    final_res[[2]] <- imports
    return(final_res)
  } else if (method == "K-Nearest Neighbor") {
    set.seed(10)
    final_res <- list()
    data_label <- cbind(labels, data)
    x <- data_label
    y <- factor(data_label[,1])
    colnames(x)[1] <- "labels"
    test.fold <- split(sample(1:length(labels)), 1:folds)

    for (mm in 1:folds) {
      test <- test.fold[[mm]]
      train_df <- x[-test, ]
      test_df <- x[test, ]

      iris.kknn <- kknn(labels~., train_df, test_df, distance = 1, kernel = "triangular")

      true_label <- dummies::dummy(test_df$labels, sep = ".")
      true_label <- data.frame(true_label)
      colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
      colnames(true_label) <- paste0(colnames(true_label), "_true")

      pred_label <- iris.kknn$prob
      colnames(pred_label) <- gsub("true","pred_SVM",colnames(true_label))

      final_df_2 <- cbind(true_label, pred_label)
      final_df2_t <- t(final_df_2)
      final_df2_t_la <- rownames(final_df2_t)
      final_df2_t_c <- as.data.frame(cbind(final_df2_t_la,final_df2_t))

      if (mm == 1) {
        final_df_m <- final_df2_t_c
        imports <- 1:10
      } else {
        final_df_m <- merge(final_df_m, final_df2_t_c, by = "final_df2_t_la", all = T)
        imports <- 1:10
      }
    }
    final_res[[1]] <- final_df_m
    final_res[[2]] <- imports
    return(final_res)
  }  else if (method == "Naive Bayes") {
    set.seed(10)
    final_res <- list()
    data_label <- cbind(as.numeric(labels), data)
    x <- data_label
    y <- factor(data_label[,1])
    colnames(x)[1] <- "labels"
    test.fold <- split(sample(1:length(labels)), 1:folds)

    for (mm in 1:folds) {
      test <- test.fold[[mm]]
      train_df <- x[-test, ]
      test_df <- x[test, ]

      model <- naiveBayes(labels ~ ., data = train_df, laplace = 3)
      pred <- predict(model, test_df)

      true_label <- dummies::dummy(test_df$labels, sep = ".")
      true_label <- data.frame(true_label)
      colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
      colnames(true_label) <- paste0(colnames(true_label), "_true")

      pred_label <- dummies::dummy(pred, sep = ".")
      pred_label <- data.frame(pred_label)
      colnames(pred_label) <- gsub(".*?\\.", "", colnames(pred_label))
      colnames(pred_label) <- paste0(colnames(pred_label), "_pred_SVM")

      final_df_2 <- cbind(true_label, pred_label)
      final_df2_t <- t(final_df_2)
      final_df2_t_la <- rownames(final_df2_t)
      final_df2_t_c <- as.data.frame(cbind(final_df2_t_la,final_df2_t))

      if (mm == 1) {
        final_df_m <- final_df2_t_c
        imports <- 1:10
      } else {
        final_df_m <- merge(final_df_m, final_df2_t_c, by = "final_df2_t_la", all = T)
        imports <- 1:10
      }
    }
    final_res[[1]] <- final_df_m
    final_res[[2]] <- imports
    return(final_res)
  } else if (method == "Linear Discriminat Analysis") {
    set.seed(10)
    final_res <- list()
    data_label <- cbind(labels, data)
    x <- data_label
    y <- factor(data_label[, 1])
    colnames(x)[1] <- "labels"
    test.fold <- split(sample(1:length(labels)), 1:folds) #ignore warning
    for (mm in 1:folds) {
      test <- test.fold[[mm]]
      train_df <- x[-test, ]
      test_df <- x[test, ]
      ld <- lda(labels~., train_df)
      z_pred <- predict(ld, test_df)
      true_label <- dummies::dummy(test_df$labels, sep = ".")
      true_label <- data.frame(true_label)
      colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
      colnames(true_label) <- paste0(colnames(true_label), "_true")
      pred_label <- z_pred$posterior
      colnames(pred_label) <- gsub("true", "pred_SVM", colnames(true_label))
      pred_label_s <- t(sapply(1:dim(pred_label)[1], function(i) ifelse(pred_label[i, ] == max(pred_label[i, ]), 1, 0)))
      final_df_2 <- cbind(true_label, pred_label_s)
      final_df2_t <- t(final_df_2)
      final_df2_t_la <- rownames(final_df2_t)
      final_df2_t_c <- as.data.frame(cbind(final_df2_t_la, final_df2_t))
      import_f <- ld$scaling[, 1]
      import_s <- import_f[order(import_f, decreasing = T)][1:10]
      import_a <- as.data.frame(cbind(row.names(t(t(import_s))), t(t(import_s))))
      colnames(import_a) <- c("V1", paste0("V_", mm))
      if (mm == 1) {
        final_df_m <- final_df2_t_c
        imports <- import_a
      } else {
        final_df_m <- merge(final_df_m, final_df2_t_c, by = "final_df2_t_la", all = T)
        imports <- merge(imports, import_a, by = "V1", all = T)
      }
    }
    final_res[[1]] <- final_df_m
    final_res[[2]] <- imports
    return(final_res)
  } else if (method == "Random Forest") {
    set.seed(10)
    final_res <- list()
    labels <- factor(labels)
    data_label <- cbind(labels, data)
    x <- data_label
    colnames(x)[1] <- "labels"
    test.fold <- split(sample(1:length(labels)), 1:folds)

    for (mm in 1:folds) {
      test <- test.fold[[mm]]
      train_df <- x[-test, ]
      test_df <- x[test, ]

      iris.rf <- randomForest(labels ~ ., data = train_df, importance = TRUE, proximity = TRUE)
      pre_iris <- predict(iris.rf, test_df)

      true_label <- dummies::dummy(test_df$labels, sep = ".")
      true_label <- data.frame(true_label)
      colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
      colnames(true_label) <- paste0(colnames(true_label), "_true")

      pred_label <- dummies::dummy(pre_iris, sep = ".")
      pred_label <- data.frame(pred_label)
      colnames(pred_label) <- gsub(".*?\\.", "", colnames(pred_label))
      colnames(pred_label) <- paste0(colnames(pred_label), "_pred_SVM")

      final_df_2 <- cbind(true_label, pred_label)
      final_df2_t <- t(final_df_2)
      final_df2_t_la <- rownames(final_df2_t)
      final_df2_t_c <- as.data.frame(cbind(final_df2_t_la, final_df2_t))

      import_f <- iris.rf$importance[, 4]
      import_s <- import_f[order(import_f, decreasing = T)][1:10]
      import_a <- as.data.frame(cbind(row.names(t(t(import_s))), t(t(import_s))))
      colnames(import_a) <- c("V1", paste0("V_", mm))

      if (mm == 1) {
        final_df_m <- final_df2_t_c
        imports <- import_a
      } else {
        final_df_m <- merge(final_df_m, final_df2_t_c, by = "final_df2_t_la", all = T)
        imports <- merge(imports, import_a, by = "V1", all = T)
      }
    }
    final_res[[1]] <- final_df_m
    final_res[[2]] <- imports
    return(final_res)
  } else if (method == "Support Vector Machine") {
    set.seed(10)
    final_res <- list()
    data_label <- cbind(as.numeric(labels), data)
    x <- data_label
    y <- factor(data_label[, 1])
    colnames(x)[1] <- "labels"
    test.fold <- split(sample(1:length(labels)), 1:folds)

    for (mm in 1:folds) {
      test <- test.fold[[mm]]
      train_df <- x[-test, ]
      test_df <- x[test, ]

      svmmodel <- svm(labels ~ ., data = train_df, kernel = "radial", type = "C-classification", probability = TRUE)
      svm_pred <- predict(svmmodel, test_df, probability = TRUE)
      svm_pred <- attr(svm_pred, "probabilities")
      colnames(svm_pred) <- paste0(colnames(svm_pred), "_pred_SVM")

      true_label <- dummies::dummy(test_df$labels, sep = ".")
      true_label <- data.frame(true_label)
      colnames(true_label) <- gsub(".*?\\.", "", colnames(true_label))
      colnames(true_label) <- paste0(colnames(true_label), "_true")

      final_df_2 <- cbind(true_label, svm_pred)
      final_df2_t <- t(final_df_2)
      final_df2_t_la <- rownames(final_df2_t)
      final_df2_t_c <- as.data.frame(cbind(final_df2_t_la, final_df2_t))

      w <- t(svmmodel$coefs) %*% svmmodel$SV
      rankingCriteria <- w * w
      ranking <- rev(sort(rankingCriteria[1, ]))
      import_f <- t(t(ranking))
      import_s <- import_f[order(import_f, decreasing = T)][1:10]
      import_a <- as.data.frame(cbind(row.names(t(t(import_f))), t(t(import_s))))
      colnames(import_a) <- c("V1", paste0("V_", mm))

      if (mm == 1) {
        final_df_m <- final_df2_t_c
        imports <- import_a
      } else {
        final_df_m <- merge(final_df_m, final_df2_t_c, by = "final_df2_t_la", all = T)
        imports <- merge(imports, import_a, by = "V1", all = T)
      }
    }
    final_res[[1]] <- final_df_m
    final_res[[2]] <- imports
    return(final_res)
  } else {
    stop(paste("Invalid method:", method))
  }
}

















