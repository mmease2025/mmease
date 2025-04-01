#'Filtering
#'
#' @param data Single-cell metabolomics data matrix
#' @param percentage Percent of Missing Values
#'
#' @return Filtered single-cell metabolomics data matrix
#' @export
#'
#' @examples
#'
#' set.seed(123)
#' data <- matrix(rnorm(100), ncol = 10)
#'
#' data <- cbind(matrix(sample(1:10, 40, replace = TRUE), ncol = 4), data)
#'
#' filted_data <- filtering(data, 0.2)
#'
filtering <- function(data, percentage) {



  if (percentage < 0 || percentage > 1) {
    stop("The percentage parameter must be between 0 and 1")
  }

  data2 <- data[, -c(1:4)]
  data4 <- data[, 1:4]

  data2[data2 == 0] <- NA

  col_f_missing <- apply(data2, 2, function(x) mean(is.na(x)))
  data2_f_missing <- data2[, col_f_missing <= percentage, drop = FALSE]

  return(cbind(data4, data2_f_missing))
}
