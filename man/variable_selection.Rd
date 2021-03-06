% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variable_selection.R
\name{variable_selection}
\alias{variable_selection}
\title{This function allows the user to select the most relevant variables thanks to
the estimation of their selection frequencies obtained by the stability
selection approach.}
\usage{
variable_selection(Y, X, square_root_inv_hat_Sigma, nb_repli = 1000,
  parallel = FALSE, nb.cores = 1)
}
\arguments{
\item{Y}{a response matrix}

\item{X}{a matrix of covariables}

\item{square_root_inv_hat_Sigma}{Estimation of the inverse of the square root of the covariance matrix
of each row of the residuals matrix obtained by the whitening function.}

\item{nb_repli}{numerical, number of replications in the stability selection}

\item{parallel}{logical, if TRUE then a parallelized version of the code is used}

\item{nb.cores}{numerical, number of cores used}
}
\value{
A data frame containing the selection frequencies of the different variables obtained
by the stability selection, the corresponding level in the design matrix and the associated
column of the observations matrix.
}
\description{
This function allows the user to select the most relevant variables thanks to
the estimation of their selection frequencies obtained by the stability
selection approach.
}
\examples{
data("copals_camera")
Y <- scale(Y[, 1:50])
X <- model.matrix(~ group + 0)
residuals <- lm(as.matrix(Y) ~ X - 1)$residuals
S12_inv <- whitening(residuals, "AR1", pAR = 1, qMA = 0)
Frequencies <- variable_selection(
  Y = Y, X = X,
  square_root_inv_hat_Sigma = S12_inv,
  nb_repli = 10, nb.cores = 1, parallel = FALSE
)
}
