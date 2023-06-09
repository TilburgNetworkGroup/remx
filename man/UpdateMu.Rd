% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{UpdateMu}
\alias{UpdateMu}
\title{UpdateMu}
\usage{
UpdateMu(p, K, TauMu, Sigma, beta)
}
\arguments{
\item{p}{an integer representing the number of random effects in the model}

\item{K}{is the number of groups in the multilevel analysis}

\item{TauMu}{is a matrix with dimensions p x p representing the covariance matrix of the prior of Mu}

\item{Sigma}{is a p x p matrix representing the covariance matrix of the random effects}

\item{beta}{is a p x K matrix representing one MCMC sample of the p random effects for each of the K groups}
}
\value{
a matrix with dimensions 1 x p representing one MCMC sample of the random-effect mean vector
}
\description{
A function to update the vector of random-effect means
}
