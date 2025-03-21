% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformation.R
\name{transformation}
\alias{transformation}
\title{Transformation}
\usage{
transformation(data, method = "G-log", lambda = 1)
}
\arguments{
\item{data}{Imputed single-cell metabolomics data matrix}

\item{method}{G-log, Log2, Log10}

\item{lambda}{lambda is a parameter used in the G - log transformation method within the transformation function, which is involved in the calculation of the g_log_transform function to adjust the result of the G - log transformation as log((data + sqrt(data^2 + lambda^2)) / 2).}
}
\value{
Transformed single-cell metabolomics data matrix
}
\description{
Transformation
}
\examples{

imputed_data <- cbind(matrix(rnorm(40), nrow = 10), matrix(rnorm(60), nrow = 10))
transformed_data <- transformation(imputed_data, method = "G-log")

}
