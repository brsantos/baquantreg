#include <RcppArmadillo.h>
#include <RcppGSL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "helper1.h"
#include "helperGIG.h"
#include "helperRD.h"
#include "helperKappa.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

// [[Rcpp::export]]
List spBayesQR(double tau, arma::colvec y, arma::mat X, int itNum, int thin,
                    arma::colvec betaValue, double sigmaValue,
                    arma::vec spCoord1, arma::vec spCoord2, double kappa,
                    double tuneP, double nugget, double priorVar){

 RNGScope scope;

 int n = X.n_rows;
 int p = X.n_cols;

 double theta, psi2, s0, n0, gama2, aux, aux2, quantEst, lambda,
          meanModel, sdModel, nTilde, sTilde, zValue, delta2, Cii,
          alphaValue;

 NumericVector sigmaSample(itNum), termsSum(n), InfPar(1), LimSup(1);

 arma::colvec b0(p), zSample(n), mu(p), resVec(n);

 arma::mat B0(p,p), SigmaMinusOne(p,p), diagU, Sigma,
  betaSample(itNum, p), vSample(itNum, n);

 /* ************************************************* */

 B0 = priorVar * B0.eye(p,p);
 b0.fill(0);

 theta = (1-2*tau)/(tau*(1-tau));
 psi2 = 2/(tau*(1-tau));

 zSample.fill(1);

 n0 = 3.0;
 s0 = 0.1;

 lambda = 0.5;
 InfPar[0] = R_NegInf; LimSup[0] = 0;
 /* *********************** */

arma::mat covMat(n, n, arma::fill::zeros), covMatInv(n, n, arma::fill::zeros),
  auxCov(n, n, arma::fill::zeros), cholCov(n,n);

arma::colvec kappaSample(itNum, arma::fill::ones);

kappaSample[0] = kappa;

double kappa1value = kappa;

auxCov = auxCov.eye();

/* *********************** */

IntegerVector seqRefresh = seq(1, 100)*(itNum/100);

 for(int k = 1; k < itNum; k++) {
    for(int j = 0; j < thin; j++) {

      if(is_true(any(k == seqRefresh))){
        Rcout << "iter = " << k << std::endl;
      }

      for(int aa = 0; aa < n; aa++)
        for(int bb = 0; bb < n; bb++)
          covMat(aa,bb) = exp(-kappa1value *
            (pow(spCoord1[aa] - spCoord1[bb],2) +
            pow(spCoord2[aa] - spCoord2[bb],2)));

      diagU = diagmat(sqrt(zSample))*(sqrt(psi2));
      cholCov = chol(covMat + nugget*auxCov.eye());
      covMatInv = cholCov.i() * cholCov.i().t();

      SigmaMinusOne = ((1/sigmaValue) * X.t() * diagU.i().t() *
        covMatInv * diagU.i() * X) + B0.i();
      Sigma = SigmaMinusOne.i();

      mu = Sigma * (B0.i()*b0 +
              (1/(sigmaValue))*(X.t() * diagU.i().t() * covMatInv *
              diagU.i() * (y - theta*zSample)));

      betaValue = mvrnormRcpp(mu, Sigma);

      resVec = y - theta*zSample - X * betaValue;
      nTilde = n0 + 3*n;
      sTilde =  arma::as_scalar(s0 + 2*sum(zSample) + resVec.t() *
        diagU.i().t() * covMatInv * diagU.i() * resVec);

      sigmaValue = rinvgammaRcpp(nTilde/2,sTilde/2);

      for(int o = 0; o < n; o++){
        Cii = arma::as_scalar(covMatInv(o,o));
        quantEst = arma::as_scalar(X.row(o) * betaValue);
        aux = y[o] - quantEst;
        aux2 = (arma::as_scalar(covMatInv.row(o) * resVec) - resVec[o]*Cii)*aux;
        delta2 = std::max((pow(aux,2.0)*Cii + aux2)/(psi2*sigmaValue), 1e-10);
        gama2 = std::max((1/sigmaValue)*(2 + pow(theta,2.0) *
          (Cii + arma::as_scalar(covMatInv.row(o)*zSample) - Cii*zSample[o])/psi2),
          1e-10);
        zSample[o] = rgigRcpp(delta2, gama2, lambda);
      }

      kappa1value = mhKappa(kappa1value, X, y, betaValue, sigmaValue,
                            zSample, tau, psi2, theta, spCoord1, spCoord2,
                            tuneP, nugget);
    }

    betaSample.row(k) = betaValue.t();
    sigmaSample[k] = sigmaValue;
    vSample.row(k) = zSample.t();
 }

 return Rcpp::List::create(
      Rcpp::Named("BetaSample") = betaSample,
      Rcpp::Named("SigmaSample") = sigmaSample,
      Rcpp::Named("vSample") = vSample,
      Rcpp::Named("kappa1") = kappaSample);
}