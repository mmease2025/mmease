#' Transformation
#'
#' @param data Imputed single-cell metabolomics data matrix
#' @param method G-log, Log2, Log10
#' @param lambda lambda is a parameter used in the G - log transformation method within the transformation function, which is involved in the calculation of the g_log_transform function to adjust the result of the G - log transformation as log((data + sqrt(data^2 + lambda^2)) / 2).
#'
#' @return Transformed single-cell metabolomics data matrix
#' @export
#'
#' @examples
#'
#' imputed_data <- cbind(matrix(rnorm(40), nrow = 10), matrix(rnorm(60), nrow = 10))
#' transformed_data <- transformation(imputed_data, method = "G-log")
#'
transformation <- function(data, method = "G-log", lambda = 1) {
  # method: "G-log", "Log2", "Log10"

  data2 <- data[, -c(1:4)]
  data4 <- data[, 1:4]

  data2[data2 == 0] <- NA

  g_log_transform <- function(data, lambda) {
    log((data + sqrt(data^2 + lambda^2)) / 2)
  }

  #train_data_t <- t(data[, -(1:4)])

  if (method == "G-log") {
    transformed_data <- g_log_transform(t(data2), lambda = lambda)
  } else if (method == "Log2") {
    transformed_data <- log2(t(data2))
  } else if (method == "Log10") {
    transformed_data <- log10(t(data2))
  } else {
    stop("Unknown transform method, please select 'G-log', 'Log2' or 'Log10'")
  }

  return(cbind(data4,t(transformed_data)))
}



















