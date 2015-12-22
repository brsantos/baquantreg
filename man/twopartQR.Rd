% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/twopartQR.R
\name{twopartQR}
\alias{twopartQR}
\title{Two part model using Bayesian quantile regression model}
\usage{
twopartQR(formula, tau = 0.5, data, itNum, thin = 1, betaValue,
  sigmaValue = 1, gammaValue, sigmaGamma = 0.5, link = 1,
  priorVar = 100, refresh = 100)
}
\arguments{
\item{formula}{a formula object, with the response on the left of a ~
operator, and the terms, separated by + operators, on the right.}

\item{tau}{Quantile of interest.}

\item{data}{a data.frame from which to find the variables defined in the
formula}

\item{itNum}{Number of iterations.}

\item{thin}{Thinning parameter.}

\item{betaValue}{Initial values for the parameter beta for the continuous part.}

\item{sigmaValue}{Initial value for the scale parameter.}

\item{gammaValue}{Initial value for the parameter gamma of the discrete
part.}

\item{sigmaGamma}{Tuning parameter for the Metropolis-Hastings step.}

\item{link}{Integer defining the link function used for the probability
model. Default is 1 for the logit link function.}

\item{priorVar}{Value that multiplies a identity matrix in the elicition
process of the prior variance of the regression parameters.}

\item{refresh}{Interval between printing a message during the iteration
process. Default is set to 100.}
}
\value{
A list with the chains of all parameters of interest.
}
\description{
This function estimates a two part model using a Bayesian quantile
regression model to describe the continous part of the conditional
distribution. The response variable is assumed to follow a mixed
discrete-continuous distribution.
}
\examples{
set.seed(1)
}
\references{
Santos and Bolfarine (2015) - Bayesian quantile regression
analysis for continuous data with a discrete component at zero.
\emph{Preprint}. \url{http://arxiv.org/abs/1511.05925}
}
