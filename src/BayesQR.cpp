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
                 arma::colvec betaValue, double sigmaValue,
                 arma::vec vSampleInit, double priorVar,
                 NumericVector hyperSigma, int refresh,
                 bool sigmaSampling,
                 bool quiet, bool tobit, bool recordLat,
                 int blocksV, bool stopOrdering, int numOrdered){

  RNGScope scope;

  int n, p;

  n = X.n_rows;
  p = X.n_cols;

  arma::mat B0(p,p);
  arma::colvec b0(p);

  B0 = priorVar * B0.eye(p,p);
  b0.fill(0);

  NumericVector sigmaSample(itNum);

  arma::colvec zSample = vSampleInit;

  double theta = (1-2*tau)/(tau*(1-tau));
  double psi2 = 2/(tau*(1-tau));

  double n0 = hyperSigma[0];
  double s0 = hyperSigma[1];

  arma::mat sigmaMinusOne(p,p), vSample(itNum, n), betaSample(itNum, p);

  arma::mat xAux(n,p);

  arma::colvec mu(p), aux(n), delta2(n);

  double gama2, lambda = 0.5, meanModel, sdModel, sTilde, termsSum;

  double nTilde;

  NumericVector InfPar(1), LimSup(1); InfPar[0] = R_NegInf; LimSup[0] = 0;

  arma::mat diagU, diagU1;

  IntegerVector seqRefresh = seq(1, itNum/refresh)*(refresh);

  arma::colvec yS = y;

  arma::uvec sortRes(n);

  for(int k = 1; k < itNum; k++){
    for(int j = 0; j < thin; j++) {

      if(!quiet){
        if(is_true(any(k+1 == seqRefresh))){
          Rcout << "Iteration = " << k+1 << std::endl;
        }
      }

      if(tobit){
        for (int aa = 0; aa < n; aa++){
          if (yS[aa] <= 0){
            meanModel = aux[aa] + theta*zSample(aa);
            sdModel = sqrt(psi2*sigmaValue*zSample(aa));
            yS[aa] = rnorm_trunc(meanModel, sdModel, InfPar[0], LimSup[0]);
          }
        }
      }

      diagU = diagmat(1/zSample)*(1/(psi2*sigmaValue));
      diagU1 = diagmat(1/zSample);

      sigmaMinusOne = (X.t() * diagU * X) + B0.i();
      arma::mat sigma = sigmaMinusOne.i();

      mu =  sigma * (B0.i()*b0 + (1/(psi2 * sigmaValue)) *
        (X.t()*diagU1 * (yS - theta*zSample)));

      betaValue = mvrnormRcpp(mu, sigma);

      aux = X * betaValue;
      gama2 = 2/sigmaValue + pow(theta,2.0)/(psi2*sigmaValue);

      delta2 = diagvec((1/(psi2*sigmaValue)) * diagmat(yS - aux) *
        diagmat(yS - aux));

      if(blocksV == 0){
        for(int o = 0; o < n; o++){
          delta2[o] = std::max(delta2[o], 1e-8);
          zSample[o] = rgigRcpp(delta2[o], gama2, lambda);
        }
      }
      else{
        NumericVector tempVec2 = floor(seq_len(blocksV) * n/blocksV);
        IntegerVector tempVec = as<IntegerVector>(tempVec2);

        // Reordering terms.
        if (!stopOrdering) sortRes = arma::sort_index(delta2);
        else if (stopOrdering && k < numOrdered) sortRes = arma::sort_index(delta2);
        // X = X.rows(sortRes);
        // yS = yS(sortRes);
        // delta2 = delta2(sortRes);
        // aux = aux(sortRes);

        for(int kk = 0; kk < blocksV; kk++){
          double zSampleAux;
          double delta_aux = 0.0;
          int oo, beg_span;
          if (kk == 0){
            oo = 0;
            beg_span = 0;
          }
          else{
            oo = tempVec[kk - 1];
            beg_span = tempVec[kk - 1];
          }

          int range = tempVec[kk] - beg_span;

          while (oo < tempVec[kk]){
            delta_aux = delta_aux + std::max(delta2(sortRes(oo)), 1e-8);
            oo++;
          }

          zSampleAux = rgigRcpp(delta_aux/range, gama2, lambda);
          zSample.elem(sortRes(arma::span(beg_span, tempVec[kk]-1))).fill(zSampleAux);
          // zSample(sortRes(arma::span(beg_span, tempVec2[kk] - 1))) = arma::ones<arma::vec>(range) * zSampleAux;
        }
      }

      if (sigmaSampling){
        termsSum = arma::as_scalar((yS - aux - theta*zSample).t() *
          diagmat(zSample).i() * (yS - aux - theta*zSample));

        nTilde = n0 + 3*n;
        sTilde =  s0 + 2*sum(zSample) + termsSum/psi2;
        sigmaValue = rinvgammaRcpp(nTilde/2,sTilde/2);
      }
    }

    betaSample.row(k) = betaValue.t();
    sigmaSample[k] = sigmaValue;
    vSample.row(k) = zSample.t();
  }

  Rcpp::List output;

  if (recordLat){
    output = Rcpp::List::create(
      Rcpp::Named("BetaSample") = betaSample,
      Rcpp::Named("SigmaSample") = sigmaSample,
      Rcpp::Named("vSample") = vSample);
  }
  else{
    output = Rcpp::List::create(
      Rcpp::Named("BetaSample") = betaSample,
      Rcpp::Named("SigmaSample") = sigmaSample);
  }

  return output;
}