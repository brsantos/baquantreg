#include <RcppArmadillo.h>
#include <RcppGSL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "helper1.h"
#include "helperGIG.h"
#include "helperGIG2.h"
#include "helperRD.h"
#include "helperKappa.h"
#include "helperAlpha.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

// [[Rcpp::export]]
List spBayesQR(double tau, arma::colvec y, arma::mat X, int itNum, int thin,
                    arma::colvec betaValue, double sigmaValue,
                    arma::vec spCoord1, arma::vec spCoord2, double kappa1value,
                    double tuneP, double alphaValue, double tuneA,
                    double priorVar, int refresh, bool quiet, double jitter,
                    bool includeAlpha, double tuneV, int kMT){

 RNGScope scope;

 int n = X.n_rows;
 int p = X.n_cols;

 double theta, psi2, s0, n0, lambda, nTilde, sTilde, delta2;

 NumericVector sigmaSample(itNum), termsSum(n), InfPar(1), LimSup(1);

 arma::colvec b0(p), zSample(n), mu(p), resVec(n), aux(n),
  Cii(n), gamma2(n);

 arma::mat B0(p,p), SigmaMinusOne(p,p), diagU, Sigma,
  betaSample(itNum, p), vSample(itNum, n), aux2(n, 1);

 /* ************************************************* */

 B0 = priorVar * B0.eye(p,p);
 b0.fill(0);

 theta = (1-2*tau)/(tau*(1-tau));
 psi2 = 2/(tau*(1-tau));

 zSample.fill(1);

 n0 = 3.0;
 s0 = 0.1;

 lambda = 0.5;
 /* *********************** */

arma::mat covMat(n, n, arma::fill::zeros),
          covMatInv(n, n, arma::fill::zeros),
          covMatAux(n, n, arma::fill::zeros),
          covMat2(n, n, arma::fill::zeros), cholCov(n,n);

arma::colvec kappaSample(itNum, arma::fill::ones),
  alphaSample(itNum, arma::fill::zeros);

kappaSample[0] = kappa1value;
if (includeAlpha) alphaSample[0] = alphaValue;
else alphaValue = 0.0;

IntegerVector seqRefresh = seq(1, itNum/refresh)*(refresh);

 for(int k = 1; k < itNum; k++) {
    for(int j = 0; j < thin; j++) {

      if(!quiet){
        if(is_true(any(k+1 == seqRefresh))){
          Rcout << "Iteration = " << k+1 << std::endl;
        }
      }

      for(int aa = 0; aa < n; aa++)
        for(int bb = 0; bb < n; bb++)
          covMat(aa,bb) = exp(-kappa1value *
            (pow(spCoord1[aa] - spCoord1[bb],2) +
            pow(spCoord2[aa] - spCoord2[bb],2)));

      diagU = diagmat(1/sqrt(sigmaValue*psi2*zSample));
      covMatAux.diag().fill(alphaValue + jitter);

      cholCov = chol((1-alphaValue)*covMat + covMatAux).i();
      covMatInv = cholCov * cholCov.t();

      SigmaMinusOne = X.t() * diagU * covMatInv * diagU * X + B0.i();
      Sigma = SigmaMinusOne.i();

      mu = Sigma * (B0.i()*b0 + (X.t() * diagU * covMatInv *
        diagU * (y - theta*zSample)));

      betaValue = mvrnormRcpp(mu, Sigma);

      resVec = y - theta*zSample - X * betaValue;
      nTilde = n0 + 3*n;
      sTilde =  arma::as_scalar(s0 + 2*sum(zSample) + (sigmaValue*psi2) *
        resVec.t() * diagU * covMatInv * diagU * resVec);

      sigmaValue = rinvgammaRcpp(nTilde/2,sTilde/2);

      aux = 1/(psi2*sigmaValue)*(y - X * betaValue);
      // aux2 = covMatInv * resVec;

      Cii = covMatInv.diag();

      // Rcout << "Cheguei aqui" << std::endl;

      for(int o = 0; o < n; o++){
//         delta2 = std::max((aux[o]*aux2(o,0) +
//           aux[o]*Cii[o]*theta*zSample[o])/(psi2*sigmaValue), 0.01);
        delta2 = (aux[o]*aux2(o,0) + aux[o]*Cii[o]*theta*zSample[o]);
//         Rcout << o << std::endl;
//         Rcout << "delta2 = " << delta2 << std::endl;
//         Rcout << "gamma2 = " << gamma2[o] << std::endl;
//         Rcout << "zSample[o] = " << zSample[o] << std::endl;
        zSample[o] = mtM(y - X * betaValue, theta, psi2,
                              sigmaValue, zSample, zSample(o), o, covMatInv,
                              tuneV, kMT);
      }

      kappa1value = mhKappa(kappa1value, spCoord1, spCoord2, resVec, diagU,
                            covMat, covMatInv,
                            tuneP, alphaValue, jitter);

      if (includeAlpha){
        alphaValue = mhAlpha(alphaValue, resVec,
                             diagU, covMat, tuneA, jitter);
      }
    }

    betaSample.row(k) = betaValue.t();
    sigmaSample[k] = sigmaValue;
    kappaSample[k] = kappa1value;
    alphaSample[k] = alphaValue;
    vSample.row(k) = zSample.t();
 }

 return Rcpp::List::create(
      Rcpp::Named("BetaSample") = betaSample,
      Rcpp::Named("SigmaSample") = sigmaSample,
      Rcpp::Named("vSample") = vSample,
      Rcpp::Named("kappa1") = kappaSample,
      Rcpp::Named("alphaSample") = alphaSample);
}