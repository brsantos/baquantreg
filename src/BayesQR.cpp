#include <RcppArmadillo.h>
#include <RcppGSL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "helper1.h"
#include "helperGIG.h"
#include "helperRD.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

// [[Rcpp::export]]
List BayesQRrcpp(double tau, arma::colvec y, arma::mat X, int itNum, int thin,
                 arma::colvec betaValue, double sigmaValue, int refresh){

  RNGScope scope;

  int n, p;

  n = X.n_rows;
  p = X.n_cols;

  // Hiperparâmetros da priori normal
  arma::mat B0(p,p);
  arma::colvec b0(p);

  B0 = 100 * B0.eye(p,p);
  b0.fill(0);

  NumericMatrix betaSample(itNum, p);
  NumericVector sigmaSample(itNum);

  arma::colvec zSample(n);

  double theta = (1-2*tau)/(tau*(1-tau));
  double psi2 = 2/(tau*(1-tau));

  zSample.fill(1);

  int n0 = 3;
  double s0 = 0.1;

  arma::mat sigmaMinusOne(p,p);

  arma::mat xAux(n,p);

  arma::colvec mu(p);

  NumericVector termsSum(n);

  double delta2, gama2, aux, lambda= 0.5, meanModel, sdModel, sTilde, zValue;

  int nTilde;

  NumericVector InfPar(1), LimSup(1); InfPar[0] = R_NegInf; LimSup[0] = 0;

  arma::mat diagU, diagU1;

  IntegerVector seqRefresh = seq(1, itNum/refresh)*(refresh);

  for(int k = 1; k < itNum; k++) {
    for(int j = 0; j < thin; j++) {

      if(is_true(any(k+1 == seqRefresh))){
        Rcout << "Iteração = " << k+1 << std::endl;
      }

      diagU = diagmat(zSample)*(psi2*sigmaValue);
      diagU1 = diagmat(zSample);

      sigmaMinusOne = (X.t() * diagU.i() * X) + B0.i();
      arma::mat sigma = sigmaMinusOne.i();

      mu =  sigma * (B0.i()*b0 + (1/(psi2 * sigmaValue)) *
        (X.t()*diagU1.i() * (y - theta*zSample)));

      betaValue = mvrnormRcpp(mu, sigma);

      gama2 = 2/sigmaValue + pow(theta,2.0)/(psi2*sigmaValue);

      for(int o = 0; o < n; o++){
        aux = arma::as_scalar(X.row(o) * betaValue);
        delta2 = pow(y[o] - aux,2.0)/(psi2*sigmaValue);

        zSample[o] = rgigRcpp(delta2, gama2, lambda);

        termsSum[o] = pow(y[o] - aux - theta*zSample[o],2.0)/zSample[o];
      }

      nTilde = n0 + 3*n;
      sTilde =  s0 + 2*sum(zSample) + sum(termsSum)/psi2;
      sigmaValue = rinvgammaRcpp(nTilde/2,sTilde/2);
    }

    for(int jj = 0; jj < p; jj++) betaSample(k,jj) = betaValue[jj];
    sigmaSample[k] = sigmaValue;
  }

  return Rcpp::List::create(
    Rcpp::Named("BetaSample") = betaSample,
    Rcpp::Named("SigmaSample") = sigmaSample);
}