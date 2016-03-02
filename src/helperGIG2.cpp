#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double logPosteriorV(arma::vec aux, double theta, double psi2, double sigma,
                     arma::vec vSample, double newV, arma::mat C, int indice){

  double output;
  arma::vec newzSample = vSample;
  newzSample(indice) = newV;

  arma::vec u = aux - theta*newzSample;
  double value = arma::as_scalar(u.t() * diagmat(1/sqrt(newzSample)) * C *
                                 diagmat(1/sqrt(newzSample)) * u);
  return - 0.5*newV - 0.5*value/(psi2*sigma)  - newV/sigma;
}

double mtM(arma::vec aux, double theta, double psi2, double sigma,
               arma::vec vSample, double curV, int indice, arma::mat C,
               double tuneV, int k){

  double output;
  int nNeg, nNegRef;
  NumericVector prop(k), xRef(k), weights(k), weightsRef(k);
  bool flagProp = true;
  for (int i = 0; i < k; i++){
    prop[i] = rnorm(1, curV, tuneV)[0];
    while(prop[i] < 0) prop[i] = rnorm(1, curV, tuneV)[0];
    weights[i] = exp(logPosteriorV(aux, theta, psi2, sigma, vSample, prop[i],
                                   C, indice));
  }
  double y = RcppArmadillo::sample(prop, 1, false, weights)[0];
  for (int ii = 0; ii < k-1; ii++){
    xRef[ii] = rnorm(1, y, tuneV)[0];
    while(xRef[ii] < 0) xRef[ii] = rnorm(1, y, tuneV)[0];
    weightsRef[ii] = exp(logPosteriorV(aux, theta, psi2, sigma, vSample, xRef[ii],
                                       C, indice));
  }
  xRef[k-1] = curV;
  weightsRef[k-1] = exp(logPosteriorV(aux, theta, psi2, sigma, vSample, xRef[k-1],
                                      C, indice));

  double probA = sum(weights)/sum(weightsRef);
  if (runif(1)[0] < probA) output = y;
  else output = curV;
  return output;
}