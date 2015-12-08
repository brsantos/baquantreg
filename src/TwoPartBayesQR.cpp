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
                      arma::colvec betaValue, double sigmaValue,
                      arma::colvec gammaValue, double sigmaGamma,
                      int link, double priorVar, int refresh){

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

   double theta, psi2, s0, n0, gama2, lambda, meanModel, sdModel, nTilde,
    sTilde, zValue, delta2, postNew, postAnt, probAcceptance, acceptRate,
    accept, thresholdLik, thetaElip, thetaElip_min, thetaElip_max;

   NumericMatrix betaSample(itNum/thin, p), gammaSample(itNum/thin, p);
   NumericVector sigmaSample(itNum/thin), termsSum(n), InfPar(1), LimSup(1);

   arma::colvec b0(p), zSample(n2), mu(p), newgammaValue(p),
    newEllipseValue(p), vetorUm(p), aux(n2), vetorZeros(p, arma::fill::zeros);

   arma::mat B0(p,p), sigmaMinusOne(p,p), diagU, diagU1, sigmaMat(p,p),
    matIdsigmaGamma(p,p);

   // Hyperparameter of the normal prior distribution
   B0 = priorVar * B0.eye(p,p);
   b0.fill(0);

   /* Constants that depend only on tau */
   theta = (1-2*tau)/(tau*(1-tau));
   psi2 = 2/(tau*(1-tau));

   /* Initializing values. */
   zSample.fill(1);

   n0 = 3;
   s0 = 0.1;

   lambda = 0.5;
   InfPar[0] = R_NegInf; LimSup[0] = 0;
   /* *********************** */

   accept = 0.0;

   vetorUm.fill(1);
   matIdsigmaGamma = sigmaGamma * diagmat(vetorUm);

   IntegerVector seqRefresh = seq(1, itNum/refresh)*(refresh);

   for(int k = 1; k < itNum; k++) {
      for(int j = 0; j < thin; j++) {

        if(is_true(any(k+1 == seqRefresh))){
          Rcout << "Iteração = " << k+1 << std::endl;
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

        diagU = diagmat(zSample)*(psi2*sigmaValue);
        diagU1 = diagmat(zSample);

        sigmaMinusOne = (X2.t() * diagU.i() * X2) + B0.i();
        sigmaMat = sigmaMinusOne.i();

        mu = sigmaMat * (B0.i()*b0 + (1/(psi2 * sigmaValue)) *
          (X2.t()*diagU1.i() * (y2 - theta*zSample)));

        betaValue = mvrnormRcpp(mu, sigmaMat);

        aux = X2 * betaValue;

        for(int o = 0; o < n2; o++){
          delta2 = pow(y2[o] - aux[o],2.0)/(psi2*sigmaValue);
          gama2 = 2/sigmaValue + pow(theta,2.0)/(psi2*sigmaValue);
          zSample[o] = rgigRcpp(delta2, gama2, lambda);

          termsSum[o] = pow(y2[o] - aux[o] - theta*zSample[o],2.0)/zSample[o];
        }

        nTilde = n0 + 3*n2;
        sTilde =  s0 + 2*sum(zSample) + sum(termsSum)/psi2;

        sigmaValue = rinvgammaRcpp(nTilde/2,sTilde/2);
      }

      for(int jj = 0; jj < p; jj++){
        betaSample(k,jj) = betaValue[jj];
        gammaSample(k,jj) = gammaValue[jj];
      }
      sigmaSample[k] = sigmaValue;
   }

   acceptRate = accept/itNum;

   return Rcpp::List::create(
        Rcpp::Named("BetaSample") = betaSample,
        Rcpp::Named("gammaSample") = gammaSample,
        Rcpp::Named("acceptRate") = acceptRate,
        Rcpp::Named("sigmaSample") = sigmaSample);
}