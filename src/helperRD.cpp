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

double rinvgauss_rcpp(double mu, double lambda){
  double Y = pow(rnorm(1, 0.0, 1.0)[0], 2.0);
  double dispersion = 1/lambda;
  double Yphi = Y * dispersion * mu;
  double X1, out, firstroot = runif(1)[0];
  if (Yphi > 5e+5) X1 = 1/Yphi;
  else X1 = 1 + Yphi/2 * (1 - pow(1 + 4/Yphi, 0.5));
  if (firstroot < 1/(1 + X1))  out = mu * X1;
  else out = mu * (1 / X1);
  return 1/out;
}

