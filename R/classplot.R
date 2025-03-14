








#'Classification_plots
#'
#' @param final_res_output
#'
#' @return AUCplot and ROCplot
#' @export
#'
#' @examples
#' \donttest{
#' }
classification_plots <- function(final_res_output) {

  library(magrittr)

  RP_curve <- function(data, force_diag = TRUE) {
    group_names <- colnames(data)[grepl("_true", colnames(data))] %>%
      gsub("(.*)_true", "\\1", .)
    method_names <- colnames(data)[grepl("_pred.*", colnames(data))] %>%
      gsub(".*_pred_(.*)", "\\1", .) %>% unique

    y_true <- data[, grepl("_true", colnames(data))]
    colnames(y_true) <- gsub("_true", "", colnames(y_true))
    y_true <- y_true[, match(group_names, colnames(y_true))]

    res_sp <- list()
    res_se <- list()
    res_auc <- list()

    for (i in seq_along(method_names)) {
      res_sp[[i]] <- list()
      res_se[[i]] <- list()
      res_auc[[i]] <- list()
      method <- method_names[i]

      y_pred <- data[, grepl(method, colnames(data))]
      colnames(y_pred) <- gsub("_pred.*", "", colnames(y_pred))
      y_pred <- y_pred[, match(group_names, colnames(y_pred))]

      for (j in seq_along(group_names)) {
        y_true_vec <- as.numeric(y_true[, j])
        y_pred_vec <- as.numeric(y_pred[, j])

        valid_idx <- which(!is.na(y_true_vec) & !is.na(y_pred_vec))
        y_true_vec <- y_true_vec[valid_idx]
        y_pred_vec <- y_pred_vec[valid_idx]

        if (length(y_true_vec) > 1) {
          roc_res <- cal_confus(y_true_vec, y_pred_vec, force_diag = force_diag)
          res_sp[[i]][[j]] <- roc_res$TPR
          res_se[[i]][[j]] <- roc_res$PPV
          res_auc[[i]][[j]] <- cal_auc(X = roc_res$TPR, Y = roc_res$PPV)
        } else {
          res_sp[[i]][[j]] <- NA
          res_se[[i]][[j]] <- NA
          res_auc[[i]][[j]] <- NA
        }
      }

      names(res_sp[[i]]) <- group_names
      names(res_se[[i]]) <- group_names
      names(res_auc[[i]]) <- group_names

      all_sp <- res_sp[[i]] %>% unlist %>% unique %>% sort(decreasing = TRUE)
      all_se <- rep(0, length(all_sp))

      for (j in seq_along(group_names)) {
        all_se <- all_se + approx(res_sp[[i]][[j]], res_se[[i]][[j]],
                                  all_sp, yleft = 1, yright = 0)$y
      }
      all_se <- all_se / length(group_names)
      res_sp[[i]]$macro <- all_sp
      res_se[[i]]$macro <- all_se
      res_auc[[i]]$macro <- cal_auc(X = res_sp[[i]]$macro, Y = res_se[[i]]$macro)

      y_true_vec_bin <- as.numeric(as.vector(as.matrix(y_true)))
      y_pred_vec_bin <- as.numeric(as.vector(as.matrix(y_pred)))
      roc_res_bin <- cal_confus(y_true_vec_bin, y_pred_vec_bin)
      res_sp[[i]]$micro <- roc_res_bin$TPR
      res_se[[i]]$micro <- roc_res_bin$PPV
      res_auc[[i]]$micro <- cal_auc(X = res_sp[[i]]$micro, Y = res_se[[i]]$micro)
    }

    names(res_sp) <- method_names
    names(res_se) <- method_names
    names(res_auc) <- method_names

    return(list(TPR = res_sp, PPV = res_se, PR = res_auc,
                Methods = method_names, Groups = group_names))
  }

  class_data <- function() {
    class.data <- list()

    final_df_m <- final_res_output[[1]]

    final_df_m[is.na(final_df_m)] <- 0
    rownames(final_df_m) <- final_df_m[, 1]
    final_df_m_t <- as.data.frame(t(final_df_m[, -1]))
    final_df_m_t <- apply(final_df_m_t, 2, as.numeric)

    roc_res <- multi_roc(final_df_m_t, force_diag = TRUE)
    if (is.null(roc_res)) {
      stop("ROC computation failed. Check the input data.")
    }
    plot_roc_df <- plot_roc_data(roc_res)
    class.data$plot_roc_df <- plot_roc_df

    if (!is.numeric(final_df_m_t)) {
      stop("Input data for RP_curve must be numeric.")
    }
    RPdata <- RP_curve(final_df_m_t, force_diag = TRUE)
    if (is.null(RPdata)) {
      stop("PR computation failed. Check the input data.")
    }
    class.data$RPdata <- RPdata

    return(class.data)
  }


  par(mfrow = c(1, 2))


  class.data <- class_data()
  plot_roc_df <- class.data$plot_roc_df
  plot_roc_df_mic <- subset(plot_roc_df, Group == "Micro")
  AUC_mic <- round(unique(plot_roc_df_mic$AUC), 3)

  plot(x = c(0, 1), y = c(0, 1), col = "grey", lty = 2, type = "l",
       xlim = c(0, 1), ylim = c(0, 1), lwd = 1,
       xlab = "1 - Specificity", ylab = "Sensitivity")

  lines(x = 1 - plot_roc_df_mic$Specificity, y = plot_roc_df_mic$Sensitivity,
        col = "#C16E71", lwd = 2)

  legend("bottomright", legend = c("Random Guess", "ROC Curve"),
         col = c("grey", "#C16E71"), lty = c(2, 1), lwd = c(1, 2), cex = 0.8)


  RPdata <- class.data$RPdata

  plot(x = c(0, 1), y = c(1, 0), col = "grey", lty = 2, type = "l",
       xlim = c(0, 1), ylim = c(0, 1), lwd = 1,
       xlab = "Recall", ylab = "Precision")

  lines(x = RPdata$TPR$SVM$micro, y = RPdata$PPV$SVM$micro,
        col = "#ABC8E5", lwd = 2)

  legend("bottomleft", legend = c("Baseline", "PR Curve"),
         col = c("grey", "#ABC8E5"), lty = c(2, 1), lwd = c(1, 2), cex = 0.8)


  par(mfrow = c(1, 1))


  AUC_R <- round(unique(plot_roc_df_mic$AUC), 3)
  AUC_ROC_text <- paste0("The AUC Value of ROC Curve is ", AUC_R)
  print(AUC_ROC_text)


  AUC_P <- round(RPdata$PR$SVM$macro, 3)
  AUC_PR_text <- paste0("The AUC Value of PR Curve is ", AUC_P)
  print(AUC_PR_text)
}









