#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double logPosteriorV(double v, double delta2, double gama2){
  if (delta2 < 0) delta2 = 0.01;
  if (gama2 < 0) gama2 = 0.01;
  return -0.5*log(v) -0.5*(delta2/v + gama2*v);
}

double mtMHRcpp(double x, double delta2, double gama2, double tuneV, int k){
  double output;
  int nNeg, nNegRef;
  NumericVector prop(k), xRef(k), weights(k), weightsRef(k);
  bool flagProp = true;
  for (int i = 0; i < k; i++){
    prop[i] = rnorm(1, x, tuneV)[0];
    while(prop[i] < 0) prop[i] = rnorm(1, x, tuneV)[0];
    weights[i] = exp(logPosteriorV(prop[i], delta2=delta2, gama2=gama2));
  }
  double y = RcppArmadillo::sample(prop, 1, false, weights)[0];
  for (int ii = 0; ii < k-1; ii++){
    xRef[ii] = rnorm(1, y, tuneV)[0];
    while(xRef[ii] < 0) xRef[ii] = rnorm(1, y, tuneV)[0];
    weightsRef[ii] = exp(logPosteriorV(xRef[ii], delta2=delta2, gama2=gama2));
  }
  xRef[k-1] = x;
  weightsRef[k-1] = exp(logPosteriorV(xRef[k-1], delta2=delta2, gama2=gama2));

  double probA = sum(weights)/sum(weightsRef);
  if (runif(1)[0] < probA) output = y;
  else output = x;
  return output;
}

