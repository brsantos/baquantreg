#include <RcppArmadillo.h>
#include <RcppGSL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "helper1.h"
#include "helperGIG.h"
#include "helperRD.h"

using namespace Rcpp;   // inline does that for us already

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

// [[Rcpp::export]]
List tpBayesQR(double tau, arma::colvec y, arma::mat X, int itNum, int thin,
              arma::vec betaValue, double sigmaValue, arma::vec vSampleInit,
              arma::vec gammaValue, double sigmaGamma,
              int link, double priorVar, int refresh, bool quiet){

   RNGScope scope;

   int n, n2, n2_aux, p;

   n = X.n_rows;
   p = X.n_cols;

   arma::colvec varInd(n, arma::fill::zeros);
   for (int dd = 0; dd < n; dd++){
     if (y[dd]==0) varInd[dd] = 1;
   }

   n2 = n - sum(varInd);
   arma::colvec y2(n2);
   arma::mat X2(n2,p);
   n2_aux = 0;
   for (int aa = 0; aa < n; aa++){
     if (varInd[aa]==0){
       y2[n2_aux] = y[aa];
       X2.row(n2_aux) = X.row(aa);
       n2_aux++;
     }
   }

   double theta, psi2, s0, n0, gama2, lambda, nTilde, sTilde, zValue,
    postNew, postAnt, probAcceptance, acceptRate, accept, termsSum;

   NumericVector InfPar(1), LimSup(1);

   arma::colvec b0(p), zSample(n2), mu(p), newgammaValue(p),
    aux(n2), sigmaSample(itNum/thin), delta2(n2);

   arma::mat betaSample(itNum/thin, p, arma::fill::zeros),
   vSample(itNum/thin, n2, arma::fill::zeros),
   gammaSample(itNum/thin, p, arma::fill::zeros), B0(p,p, arma::fill::zeros),
    sigmaMinusOne(p,p,arma::fill::zeros), diagU, diagU1,
    sigmaMat(p,p,arma::fill::zeros), matIdsigmaGamma(p,p,arma::fill::zeros);

   // Hyperparameter of the normal prior distribution
   B0 = priorVar * B0.eye(p,p);
   b0.fill(0);

   /* Constants that depend only on tau */
   theta = (1-2*tau)/(tau*(1-tau));
   psi2 = 2/(tau*(1-tau));

   /* Initializing values. */
   zSample = vSampleInit;

   n0 = 3;
   s0 = 0.1;

   lambda = 0.5;
   InfPar[0] = R_NegInf; LimSup[0] = 0;
   /* *********************** */

   accept = 0.0;

   matIdsigmaGamma.diag().fill(sigmaGamma);

   IntegerVector seqRefresh = seq(1, itNum/refresh)*(refresh);

   for(int k = 1; k < itNum; k++) {
      for(int j = 0; j < thin; j++) {

        if(!quiet){
          if(is_true(any(k+1 == seqRefresh))){
            Rcout << "Iteration = " << k+1 << std::endl;
          }
        }

        // ***************************************** //
        // Updating the point mass part of the model //
        newgammaValue = mvrnormRcpp(gammaValue, matIdsigmaGamma);

        postAnt = postDensMHstep(gammaValue, link, X, varInd, b0, B0);
        postNew = postDensMHstep(newgammaValue, link, X, varInd, b0, B0);

        probAcceptance = std::min(0.0,postNew-postAnt);
        if(log(runif(1)[0])<probAcceptance){
          gammaValue = newgammaValue;
          accept++;
        }

        diagU = diagmat(1/zSample)/(psi2*sigmaValue);
        diagU1 = diagmat(1/zSample);

        sigmaMinusOne = (X2.t() * diagU * X2) + B0.i();
        sigmaMat = sigmaMinusOne.i();

        mu = sigmaMat * (B0.i()*b0 + (1/(psi2 * sigmaValue)) *
          (X2.t()*diagU1 * (y2 - theta*zSample)));

        betaValue = mvrnormRcpp(mu, sigmaMat);

        aux = X2 * betaValue;
        gama2 = 2/sigmaValue + pow(theta,2.0)/(psi2*sigmaValue);

        delta2 = diagvec((1/(psi2*sigmaValue)) * diagmat(y2 - aux) *
          diagmat(y2 - aux));

        for(int o = 0; o < n2; o++){
          delta2[o] = std::max(delta2[o], 1e-10);
          zSample[o] = rgigRcpp(delta2[o], gama2, lambda);
        }

        termsSum = arma::as_scalar((y2 - aux - theta*zSample).t() *
          diagmat(zSample).i() * (y2 - aux - theta*zSample));

        nTilde = n0 + 3*n2;
        sTilde =  s0 + 2*sum(zSample) + termsSum/psi2;

        sigmaValue = rinvgammaRcpp(nTilde/2,sTilde/2);
      }

      betaSample.row(k) = betaValue.t();
      gammaSample.row(k) = gammaValue.t();
      vSample.row(k) = zSample.t();
      sigmaSample[k] = sigmaValue;
   }

   acceptRate = accept/itNum;

   return Rcpp::List::create(
        Rcpp::Named("BetaSample") = betaSample,
        Rcpp::Named("GammaSample") = gammaSample,
        Rcpp::Named("acceptRate") = acceptRate,
        Rcpp::Named("SigmaSample") = sigmaSample,
        Rcpp::Named("vSample") = vSample);
}