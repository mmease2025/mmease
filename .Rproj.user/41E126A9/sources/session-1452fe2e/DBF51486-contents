








#'Classification_plots
#'
#' @param final_res_output
#'
#' @return Receiver operating characteristic (ROC) curve and area under the curve (AUC) value
#' @export
#'
#' @examples
#' \donttest{
#' }
classification_plots <- function(final_res_output) {

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

    return(class.data)
  }

  class.data <- class_data()
  plot_roc_df <- class.data$plot_roc_df
  plot_roc_df_mic <- subset(plot_roc_df, Group == "Micro")
  AUC_mic <- round(unique(plot_roc_df_mic$AUC), 3)

  # 设置较小的边距
  par(mar = c(4, 4, 2, 2))

  plot(x = c(0, 1), y = c(0, 1), col = "grey", lty = 2, type = "l",
       xlim = c(0, 1), ylim = c(0, 1), lwd = 1,
       xlab = "1 - Specificity", ylab = "Sensitivity",
       cex.lab = 0.8, cex.axis = 0.8)  # 调整标签和坐标轴字体大小

  lines(x = 1 - plot_roc_df_mic$Specificity, y = plot_roc_df_mic$Sensitivity,
        col = "#C16E71", lwd = 2)

  legend("bottomright", legend = c("Random Guess", "ROC Curve"),
         col = c("grey", "#C16E71"), lty = c(2, 1), lwd = c(1, 2), cex = 0.6)  # 调整图例字体大小

  # 调整位置使 AUC 值更居中显示在图的上方区域
  x_text <- mean(par("usr")[1:2])
  y_text <- max(par("usr")[3:4]) * 0.9
  text(x = x_text, y = y_text, labels = paste0("AUC = ", AUC_mic), pos = 3, cex = 0.8)

  AUC_R <- round(unique(plot_roc_df_mic$AUC), 3)
  AUC_ROC_text <- paste0("The AUC Value of ROC Curve is ", AUC_R)
  print(AUC_ROC_text)
}









