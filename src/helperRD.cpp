#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double rinvgammaRcpp(double shape, double scale){
  double out = pow(rgamma(1, shape, 1/scale)[0],-1);
  return out;
}

arma::colvec mvrnormRcpp(arma::colvec mu, arma::mat Sigma){
  int p = mu.size();
  arma::colvec eV(p), X(p);
  arma::mat eVec(p,p);

  eig_sym(eV, eVec, Sigma);

  X = rnorm(p);

  arma::mat produto = eVec * diagmat(sqrt(eV));

  arma::colvec output;

  output =  mu + produto * X;
  return output;
}

double Flaplace (double predictor, double tau, double sigma){
  double output;
  if (predictor > 0) output = tau * exp(-((1-tau)/sigma)*predictor);
  else output = 1 - (1-tau)*exp((tau/sigma)*predictor);
  return output;
}
