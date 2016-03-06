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
  return - 0.5*log(newV) - 0.5*value/(psi2*sigma)  - newV/sigma;
}

double logPosteriorV2(arma::vec aux, double theta, double psi2, double sigma,
                     arma::vec vSample, arma::mat C){

  double output;

  arma::vec u = aux - theta*vSample;
  double value = arma::as_scalar(u.t() * diagmat(1/sqrt(vSample)) * C *
                                 diagmat(1/sqrt(vSample)) * u);
  return - 0.5*sum(log(vSample)) - 0.5*value/(psi2*sigma)  -
    (1/sigma)*sum(vSample);
}

double mtM(arma::vec aux, double theta, double psi2, double sigma,
               arma::vec vSample, double curV, int indice, arma::mat C,
               double tuneV, int k){

  double output;
  int nNeg, nNegRef;
  NumericVector prop(k), xRef(k), weights(k), weightsRef(k);
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

NumericVector calcLogWeights(NumericVector weights){
  double a0 = max(weights);
  NumericVector output = exp(weights - a0)/sum(exp(weights - a0));
  return output;
}

NumericVector mtM2(arma::vec aux, double theta, double psi2, double sigma,
           arma::vec vSample, int n, arma::mat C, double tuneV, int k){

  NumericVector output;
  int nNeg, nNegRef;
  IntegerVector indices = seq_len(k);
  NumericMatrix prop(n, k), xRef(n, k);
  NumericVector weights(k), weightsRef(k);
  for (int i = 0; i < k; i++){
    prop(_,i) = rexp(n, tuneV);
    weights[i] = logPosteriorV2(aux, theta, psi2, sigma, prop(_,i), C) +
      sum(log(dexp(prop(_,i), tuneV)));
  }

  NumericVector weights2 = calcLogWeights(weights);

  int y = RcppArmadillo::sample(indices, 1, false, weights2)[0];
  for (int ii = 0; ii < k-1; ii++){
    xRef(_,ii) = rexp(n, tuneV);
    weights[ii] = logPosteriorV2(aux, theta, psi2, sigma, xRef(_,ii), C) +
      sum(log(dexp(xRef(_,ii), tuneV)));
  }
  xRef(_,k-1) = NumericVector(vSample.begin(), vSample.end());
  weightsRef[k-1] = logPosteriorV2(aux, theta, psi2, sigma, xRef(_,k-1), C) +
    sum(log(dexp(xRef(_,k-1), tuneV)));

  NumericVector weightsRef2 = calcLogWeights(weightsRef);

  double probA = sum(weights2)/sum(weightsRef2);
  if (runif(1)[0] < probA) output = prop( _, y-1);
  else output = NumericVector(vSample.begin(), vSample.end());
  return output;
}