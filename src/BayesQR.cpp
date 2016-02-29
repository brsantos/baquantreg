#include <RcppArmadillo.h>
#include <RcppGSL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "helper1.h"
#include "helperGIG.h"
#include "helperRD.h"
#include "helperTN.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

// [[Rcpp::export]]
List BayesQR(double tau, arma::colvec y, arma::mat X, int itNum, int thin,
                 arma::colvec betaValue, double sigmaValue, double priorVar,
                 int refresh, bool quiet, bool tobit){

  RNGScope scope;

  int n, p;

  n = X.n_rows;
  p = X.n_cols;

  arma::mat B0(p,p);
  arma::colvec b0(p);

  B0 = priorVar * B0.eye(p,p);
  b0.fill(0);

  NumericVector sigmaSample(itNum);

  arma::colvec zSample(n);

  double theta = (1-2*tau)/(tau*(1-tau));
  double psi2 = 2/(tau*(1-tau));

  zSample.fill(1);

  int n0 = 3;
  double s0 = 0.1;

  arma::mat sigmaMinusOne(p,p), vSample(itNum, n), betaSample(itNum, p);

  arma::mat xAux(n,p);

  arma::colvec mu(p), aux(n), delta2(n);

  double gama2, lambda = 0.5, meanModel, sdModel, sTilde, zValue, termsSum;

  int nTilde;

  NumericVector InfPar(1), LimSup(1); InfPar[0] = R_NegInf; LimSup[0] = 0;

  arma::mat diagU, diagU1;

  IntegerVector seqRefresh = seq(1, itNum/refresh)*(refresh);

  arma::colvec yS = y;

  for(int k = 1; k < itNum; k++){
    for(int j = 0; j < thin; j++) {

      if(!quiet){
        if(is_true(any(k+1 == seqRefresh))){
          Rcout << "Iteration = " << k+1 << std::endl;
        }
      }

      if(tobit){
        for (int aa = 0; aa < n; aa++){
          meanModel = aux[aa] + theta*zSample(aa);
          sdModel = sqrt(psi2*sigmaValue*zSample(aa));
          yS[aa] = rnorm_trunc(meanModel, sdModel, InfPar[0], LimSup[0]);
        }
      }

      diagU = diagmat(zSample)*(psi2*sigmaValue);
      diagU1 = diagmat(zSample);

      sigmaMinusOne = (X.t() * diagU.i() * X) + B0.i();
      arma::mat sigma = sigmaMinusOne.i();

      mu =  sigma * (B0.i()*b0 + (1/(psi2 * sigmaValue)) *
        (X.t()*diagU1.i() * (yS - theta*zSample)));

      betaValue = mvrnormRcpp(mu, sigma);

      aux = X * betaValue;
      gama2 = 2/sigmaValue + pow(theta,2.0)/(psi2*sigmaValue);

      delta2 = diagvec((1/(psi2*sigmaValue)) * diagmat(yS - aux) *
        diagmat(yS - aux));

      for(int o = 0; o < n; o++){
        delta2[o] = std::max(delta2[o], 1e-10);
        zSample[o] = rgigRcpp(delta2[o], gama2, lambda);
      }

      termsSum = arma::as_scalar((yS - aux - theta*zSample).t() *
        diagmat(zSample).i() * (yS - aux - theta*zSample));

      nTilde = n0 + 3*n;
      sTilde =  s0 + 2*sum(zSample) + termsSum/psi2;
      sigmaValue = rinvgammaRcpp(nTilde/2,sTilde/2);
    }

    betaSample.row(k) = betaValue.t();
    sigmaSample[k] = sigmaValue;
    vSample.row(k) = zSample.t();
  }

  return Rcpp::List::create(
    Rcpp::Named("BetaSample") = betaSample,
    Rcpp::Named("SigmaSample") = sigmaSample,
    Rcpp::Named("vSample") = vSample);
}