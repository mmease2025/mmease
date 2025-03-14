#' Batch_correction
#'
#' @param data Normalized single-cell metabolomics data matrix
#' @param method ComBat,Limma
#'
#' @return Corrected single-cell metabolomics data matrix
#' @export
#'
#' @examples
#'
#' normalized_data <- cbind(matrix(rnorm(40), nrow = 10), matrix(rnorm(60), nrow = 10))
#'
#' batch_info <- factor(sample(1:2, 10, replace = TRUE))
#' normalized_data <- cbind(normalized_data, batch_info)
#'
#' corrected_data <- batch_correction(t(normalized_data), method = "ComBat")
#'
batch_correction <- function(data, method = "ComBat") {
  # method: "ComBat", "Limma"

  batch <- factor(data[, 4])
  data2 <- data[, -c(1:4)]
  data4 <- data[, 1:4]

  data2[data2 == 0] <- NA


  combat_correction <- function(data2, batch) {
    requireNamespace("sva", quietly = TRUE)
    combat_data <- tryCatch({
      sva::ComBat(dat = data2, batch = batch, par.prior = TRUE, prior.plots = FALSE, mean.only = FALSE)
    }, error = function(e) {
      message("ComBat error: ", e$message)
      return(NULL)
    })
    return(combat_data)
  }

  limma_correction <- function(data2, batch) {
    requireNamespace("limma", quietly = TRUE)
    design <- stats::model.matrix(~batch)
    combat_data <- tryCatch({
      fit <- limma::lmFit(data2, design)
      corrected_data <- limma::removeBatchEffect(data2, batch = batch, design = design)
      return(corrected_data)
    }, error = function(e) {
      message("Limma error: ", e$message)
      return(NULL)
    })
    return(combat_data)
  }

  batch <- factor(data[, 4])
  if (method == "ComBat") {
    return(as.matrix(t(rbind(t(data4),combat_correction(t(data2), batch)))))
  } else if (method == "Limma") {
    return(as.matrix(t(rbind(t(data4),limma_correction(t(data2), batch)))))
  } else {
    stop("Unknown batch effect correction method, select 'ComBat' or 'Limma'")
  }
}




