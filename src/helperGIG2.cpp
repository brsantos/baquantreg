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


NumericVector calcLogWeights(NumericVector weights){
  double a0 = max(weights);
  NumericVector output = exp(weights - a0)/sum(exp(weights - a0));
  return output;
}

double mtM(arma::vec aux, double theta, double psi2, double sigma,
               arma::vec vSample, double curV, int indice, arma::mat C,
               double tuneV, int k){

  double output;
  int nNeg, nNegRef;
  NumericVector prop(k), xRef(k), weights(k), weights2(k), weightsRef(k);
  prop = rexp(k, tuneV);
  for (int i = 0; i < k; i++){
    weights[i] = logPosteriorV(aux, theta, psi2, sigma, vSample, prop[i],
                                   C, indice);
  }

  weights += log(dexp(prop, tuneV));
  weights2 = calcLogWeights(weights);

  double y = RcppArmadillo::sample(prop, 1, false, weights2)[0];
  for (int ii = 0; ii < k-1; ii++){
    xRef[ii] = rexp(1, tuneV)[0];
    weightsRef[ii] = logPosteriorV(aux, theta, psi2, sigma, vSample, xRef[ii],
                                       C, indice);
  }
  xRef[k-1] = curV;
  weightsRef[k-1] = logPosteriorV(aux, theta, psi2, sigma, vSample, xRef[k-1],
                                      C, indice);

  weightsRef += log(dexp(xRef, tuneV));

  double probA = (max(weights) + log(sum(exp(weights - max(weights))))) -
    (max(weightsRef) + sum(exp(weightsRef - max(weightsRef))));
  if (runif(1)[0] < exp(probA)) output = y;
  else output = curV;
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
    weightsRef[ii] = logPosteriorV2(aux, theta, psi2, sigma, xRef(_,ii), C) +
      sum(log(dexp(xRef(_,ii), tuneV)));
  }
  xRef(_,k-1) = NumericVector(vSample.begin(), vSample.end());
  weightsRef[k-1] = logPosteriorV2(aux, theta, psi2, sigma, xRef(_,k-1), C) +
    sum(log(dexp(xRef(_,k-1), tuneV)));

  NumericVector weightsRef2 = calcLogWeights(weightsRef);

  double probA = (max(weights) + log(sum(exp(weights - max(weights))))) -
    (max(weightsRef) + sum(exp(weightsRef - max(weightsRef))));

  if (runif(1)[0] < exp(probA)) output = prop( _, y-1);
  else output = NumericVector(vSample.begin(), vSample.end());
  return output;
}
