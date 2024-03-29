% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/latent_hellinger.R
\name{latent_hellinger}
\alias{latent_hellinger}
\title{Hellinger distance for the latent variables of the estimation process.}
\usage{
latent_hellinger(
  object,
  burnin = 50,
  plot_div = TRUE,
  scales.free = FALSE,
  all.obs = TRUE,
  obs = 1
)
}
\arguments{
\item{object}{This is an object of the class "bqr", produced by a call
to the bayesQR function.}

\item{burnin}{Initial part of the chain, which is to be discarded. Default
value is 50.}

\item{plot_div}{If TRUE, the function prints the plot with all probabilities.
Default is set to TRUE.}

\item{scales.free}{If FALSE, the default, then all plots will use the same
y scale. If TRUE, for each tau the plot will use the best possible scale
in order to visualize the probability information for all observations.}

\item{all.obs}{if TRUE, calculates KL divergence for all observations}

\item{obs}{if \code{all.obs} is FALSE, specifies the observation to
calculate the KL divergence.}
}
\value{
Prints a plot of the posterior probability of being an outlier for
 all observations and returns a data.frame with all values of the
 probabilities.
}
\description{
Returns the Hellinger distance for the latent variables, which are
 used in estimation process, and it could indicate possible outlying
 observations.
}
