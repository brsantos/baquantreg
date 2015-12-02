#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double rinvgammaRcpp(double shape, double scale){
  double out = pow(rgamma(1, shape, 1/scale)[0],-1);
  return out;
}

colvec mvrnormRcpp(colvec mu, mat Sigma){
  int p = mu.size();
  colvec eV(p), X(p);
  mat eVec(p,p);

  eig_sym(eV, eVec, Sigma);

  X = rnorm(p);

  mat produto = eVec * diagmat(sqrt(eV));

  colvec output;

  output =  mu + produto * X;
  return output;
}

double Flaplace (double predictor, double tau, double sigma){
  double output;
  if (predictor > 0) output = tau * exp(-((1-tau)/sigma)*predictor);
  else output = 1 - (1-tau)*exp((tau/sigma)*predictor);
  return output;
}
