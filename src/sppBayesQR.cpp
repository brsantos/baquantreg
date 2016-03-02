#include <RcppArmadillo.h>
#include <RcppGSL.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "helper1.h"
#include "helperGIG.h"
#include "helperRD.h"
#include "helperKappa.h"
#include "helperAlpha.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

// [[Rcpp::export]]
List sppBayesQR(double tau, arma::colvec y, arma::mat X, int itNum,
                   int thin, arma::colvec betaValue, double sigmaValue,
                   arma::vec spCoord1, arma::vec spCoord2, double kappa1value,
                   double tuneP, arma::uvec indices, int m,
                   double alphaValue, double tuneA, double priorVar,
                   bool quiet, int refresh, double jitter){

   RNGScope scope;

   int n = X.n_rows;
   int p = X.n_cols;

   double theta, psi2, s0, n0, lambda, nTilde, sTilde, delta2;

   NumericVector sigmaSample(itNum), termsSum(n);

   arma::colvec b0(p), zSample(n), mu(p), resVec(n), Cii(n),
    gamma2(n), aux(n), aux2(n);

   arma::mat B0(p,p), SigmaMinusOne(p,p), diagU, Sigma,
    betaSample(itNum, p), vSample(itNum, n), sigmaDot(n, n);

   B0 = priorVar * B0.eye(p,p);
   b0.fill(0);

   theta = (1-2*tau)/(tau*(1-tau));
   psi2 = 2/(tau*(1-tau));

   zSample.fill(1);

   n0 = 3.0;
   s0 = 0.1;

   lambda = 0.5;

  arma::mat covMatAux(n, m, arma::fill::zeros), covMat(n, n, arma::fill::zeros),
    covMat2(m, m, arma::fill::zeros), covMatInv(m, m, arma::fill::zeros),
    auxCov(m, m, arma::fill::zeros), cholCov(m,m), cholCov2(m,m), linPred(n, p),
    matAux(n,n, arma::fill::zeros), matM(m, m, arma::fill::zeros),
    matM2(m, m, arma::fill::zeros), CovCov(n, n, arma::fill::zeros);

  arma::colvec kappaSample(itNum), alphaSample(itNum);

  kappaSample[0] = kappa1value;
  alphaSample[0] = alphaValue;

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

        covMat2 = covMat.submat(indices, indices);
        covMatAux = covMat.cols(indices);

        auxCov.diag().fill(alphaValue + jitter);

        // Rcout << "Primeira matriz" << std::endl;
        cholCov = arma::chol((1-alphaValue)*covMat2 + auxCov).i();
        covMatInv =  cholCov * cholCov.t();
        matAux = covMatAux*covMatInv*covMatAux.t();
        diagU = diagmat(sqrt(1/zSample));
        sigmaDot = diagmat(alphaValue +
          (1-alphaValue)*(covMat.diag() - matAux.diag())).i();

//         sigmaDot.i().print("SigmaDot = ");
//
//         Rcout << "Segunda matriz" << std::endl;

//         cholCov2 =  arma::chol((1-alphaValue)*covMat2 +
//           covMatAux.t()*sigmaDot*covMatAux).i();
//         matM = cholCov2 * cholCov2.t();

        cholCov2 = (1-alphaValue)*covMat2 + covMatAux.t()*sigmaDot*covMatAux;
        // cholCov2.row(1).print("Primeira linha");

        matM = arma::chol(cholCov2 + auxCov).i();
        matM2 = matM * matM.t();

        CovCov = sigmaDot - sigmaDot * covMatAux * matM2 * covMatAux.t() *
          sigmaDot;

        SigmaMinusOne = ((1/(sigmaValue*psi2)) * X.t() * diagU * CovCov *
          diagU * X) + B0.i();
        Sigma = SigmaMinusOne.i();

        mu = Sigma * ((1/(sigmaValue*psi2))*(X.t() * diagU * CovCov *
                diagU * (y - theta*zSample)));

        betaValue = mvrnormRcpp(mu, Sigma);

        resVec = y - theta*zSample - X * betaValue;

        nTilde = n0 + 3*n;
        sTilde =  arma::as_scalar(s0 + 2*sum(zSample) + resVec.t() * diagU *
          CovCov * diagU * resVec);

        sigmaValue = rinvgammaRcpp(nTilde/2,sTilde/2);

        Cii = CovCov.diag();

        aux = y - X * betaValue;
        aux2 = (aux.t() * CovCov * resVec).t();

        gamma2 = 2/sigmaValue + (pow(theta,2.0)/(psi2*sigmaValue))*Cii;

        for(int o = 0; o < n; o++){
          delta2 = std::max((aux2[o] +
            aux[o]*Cii[o]*theta*zSample[o])/(psi2*sigmaValue), 1e-10);
          zSample[o] = rgigRcpp(delta2, gamma2[o], lambda);
        }

        kappa1value = mhKappa2(kappa1value, spCoord1, spCoord2, resVec, diagU,
                               covMat, CovCov,
                               tuneP, alphaValue, jitter, indices, m);

        alphaValue = mhAlpha2(alphaValue, resVec, diagU, covMat, covMat2,
                              covMatAux, tuneA, jitter, indices, m);

      }

      betaSample.row(k) = betaValue.t();
      vSample.row(k) = zSample.t();
      sigmaSample[k] = sigmaValue;
      kappaSample[k] = kappa1value;
      alphaSample[k] = alphaValue;
   }

   return Rcpp::List::create(
        Rcpp::Named("BetaSample") = betaSample,
        Rcpp::Named("SigmaSample") = sigmaSample,
        Rcpp::Named("vSample") = vSample,
        Rcpp::Named("kappa1") = kappaSample,
        Rcpp::Named("alphaSample") = alphaSample);
}
