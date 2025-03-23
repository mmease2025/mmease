

#' Imputation
#'
#' @param data Filtered single-cell metabolomics data matrix
#' @param method knn,1/5 of minimum positive value
#' @param k k represents the number of nearest neighbors considered for imputing missing values in the KNN algorithm. A small k may make the imputation unstable due to outliers, while a large k may introduce bias by including dissimilar samples.
#' @param rowmax rowmax is a proportion indicating the maximum allowable proportion of missing values per row. Rows with a higher proportion won't be used for calculating neighbors, avoiding the influence of rows with excessive missing values.
#' @param colmax colmax is a proportion defining the maximum allowable proportion of missing values per column. Columns exceeding this proportion won't be imputed, preventing inaccurate imputation for columns with too many missing values.
#' @param maxp maxp is an integer that limits the maximum number of features considered when calculating neighbors. A small maxp may overlook important features, while a large one increases computational complexity and may introduce noise.
#'
#' @return Imputed single-cell metabolomics data matrix
#' @export
#'
#' @examples
#'
#' set.seed(123)
#' filted_data <- matrix(rnorm(100), ncol = 10)
#'
#' filted_data <- cbind(matrix(sample(1:10, 40, replace = TRUE), ncol = 4), filted_data)
#'
#' imputed_data <- imputation(filted_data, method = "1/5 of minimum positive value")
#'
#' imputed_data <- imputation(filted_data, method = "KNN")
#'
imputation <- function(data, method = "KNN", k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500) {
  # method: "KNN" "back"

  data2 <- data[, -c(1:4)]
  data4 <- data[, 1:4]

  data2[data2 == 0] <- NA

  back_impute <- function(x) {
    filterdata <- t(x)
    x <- filterdata
    x <- x[apply(x, 1, function(y) !all(is.na(y))), ]
    filterdata <- x
    filterdata[is.na(filterdata)] <- min(filterdata[filterdata > 0], na.rm = TRUE) / 5
    return(t(filterdata))
  }

  if (method == "KNN") {
    library(impute)
    filter_train_data <- as.matrix(t(data2))
    data.imputed <- impute.knn(filter_train_data, k = k, rowmax = rowmax, colmax = colmax, maxp = maxp, rng.seed = 1024)
    imput_m <- t(data.imputed$data)
  } else if (method == "1/5 of minimum positive value") {
    imput_m <- back_impute(data2)
  } else {
    stop("For unknown interpolation method, please select 'KNN' or '1/5 of minimum positive value'")
  }

  return(cbind(data4, imput_m))
}
