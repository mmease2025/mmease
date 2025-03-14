#' Normalization
#'
#' @param data Transformed single-cell metabolomics data matrix
#' @param method 'Auto Scaling', 'Cyclic Loess', 'Mean', 'Median', 'MSTUS' , 'SIS'
#' @param sis_column sis_column is a parameter in the normalization function that is specifically used when the normalization method is set to "SIS". It indicates the column number in the data matrix that will be used as the reference vector for the internal standard normalization carried out by the SIS function.
#'
#' @return normalized single-cell metabolomics data matrix
#' @export
#'
#' @examples
#'
#' transformed_data <- cbind(matrix(rnorm(40), nrow = 10), matrix(rnorm(60), nrow = 10))
#' normalized_data <- normalization(transformed_data, method = "MSTUS")
#'
normalization <- function(data, method = "Auto Scaling", sis_column = 2) {
  # method: "Auto Scaling", "Cyclic Loess", "Mean", "Median", "MSTUS", "SIS"

  library(metabolomics)

  data2 <- data[, -c(1:4)]
  data4 <- data[, 1:4]

  data2[data2 == 0] <- NA

  auto_scaling <- function(data) {
    centered.data <- data - apply(data, 1, mean)
    scaling.auto <- apply(data, 1, sd)
    return(centered.data / scaling.auto)
  }

  cyclic_loess <- function(data) {
    result <- tryCatch({
      fastlo(data)
    }, error = function(e) {
      message("Error in fastlo: ", e$message)
      NULL
    })
    return(result)
  }

  mean_normalization <- function(data) {
    inputdata <- data.frame(as.factor(rep("sample", ncol(data))), t(data))
    norm_mean <- Normalise(inputdata, method = "mean")
    return(t(norm_mean$output[,-1]))
  }

  median_normalization <- function(data) {
    inputdata <- data.frame(as.factor(rep("sample", ncol(data))), t(data))
    norm_med <- Normalise(inputdata, method = "median")
    return(t(norm_med$output[,-1]))
  }

  mstuds_normalization <- function(data) {
    data_sum <- matrix(colSums(data), nrow = 1)
    area.uni <- matrix(rep(1, nrow(data)), ncol = 1) %*% data_sum
    return(data / area.uni)
  }

  SIS <- function(data, nc) {
    library(metabolomics)
    norm_is <- Normalise(data, method = "is", refvec=data[, nc[1]])
    return(norm_is$output)
  }

  if (method == "Auto Scaling") {
    return(cbind(data4,auto_scaling(data2)))
  } else if (method == "Cyclic Loess") {
    return(cbind(data4,cyclic_loess(data2)))
  } else if (method == "Mean") {
    return(cbind(data4,mean_normalization(data2)))
  } else if (method == "Median") {
    return(cbind(data4,median_normalization(data2)))
  } else if (method == "MSTUS") {
    return(cbind(data4,mstuds_normalization(data2)))
  } else if (method == "SIS") {
    return(cbind(data4,(SIS(data2, sis_column))))
  } else {
    stop("Unknown normalization method, please select 'Auto Scaling', 'Cyclic Loess', 'Mean', 'Median', 'MSTUS' or 'SIS'")
  }
}





