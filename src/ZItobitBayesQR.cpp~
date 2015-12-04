#include <RcppArmadillo.h>
#include "helper1.h"
#include "helperGIG.h"
#include "helperTN.h"
#include "helperRD.h"

using namespace Rcpp;   // inline does that for us already

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

// [[Rcpp::export]]
List ziTobitBayesQR(double tau, arma::colvec y, arma::mat X, int itNum, int thin,
                      arma::colvec betaValue, double sigmaValue, arma::colvec betaZeroValue,
                      double sigmaBetaZero, int link){

   RNGScope scope;

   int n, n1, n2, l, m, p, cont_n2, cont_l, cont_m;
   cont_n2 = 0;
   cont_l = 0;
   cont_m = 0;

   n = X.n_rows;
   p = X.n_cols;

   arma::colvec varInd(n, fill::zeros);
   for (int dd = 0; dd < n; dd++){
     if (y[dd]==0){
       varInd[dd] = 1.0;
     }
   }

   n1 = sum(varInd);
   n2 = n - n1;

   double fObsLap, probCens;
   arma::colvec censInd(n, fill::zeros);

   double theta, psi2, s0, n0, gama2, lambda, meanModel, sdModel, nTilde, sTilde, zValue,
            postNew, postAnt, probAcceptance, acceptRate, accept, termsSum;

   acceptRate = 0;

   arma::mat betaSample(itNum/thin, p, arma::fill::zeros), betaZeroSample(itNum/thin, p, arma::fill::zeros);
   arma::colvec sigmaSample(itNum/thin, arma::fill::zeros);
   NumericVector InfPar(1), LimSup(1);

   arma::colvec b0(p), mu(p), newBetaZeroValue(p), vetorUm(p, arma::fill::ones), delta2(n);

   arma::mat B0(p,p), sigmaMinusOne(p,p), diagU, sigmaMat(p,p), matIdsigmaBetaZero(p,p);

   // Hiperparâmetros da priori normal
   B0 = 100 * B0.eye(p,p);
   b0.fill(0);

   /* Constantes que dependem somente do tau estimado */
   theta = (1-2*tau)/(tau*(1-tau));
   psi2 = 2/(tau*(1-tau));

   /* Inicializando valores. */
   n0 = 3;
   s0 = 0.1;

   lambda = 0.5;
   InfPar[0] = R_NegInf; LimSup[0] = 0;
   /* *********************** */

   accept = 0.0;

   matIdsigmaBetaZero = sigmaBetaZero * diagmat(vetorUm);

   IntegerVector seqRefresh = seq(1, itNum/100)*100;

   arma::colvec zSample(n, arma::fill::ones), aux(n), aux2;

   sigmaSample[0] = sigmaValue;

   arma::mat matrizIndCens(itNum/thin, n, arma::fill::zeros), matrizIndNCens(itNum/thin, n, arma::fill::zeros);

   arma::uvec ind_n2(n2, arma::fill::zeros), ind_l(n1, arma::fill::zeros), ind_m(n1, arma::fill::zeros);

   double probZero;

   // Laço que vai fazer a atualização dos estados das cadeias.
   for(int k = 1; k < itNum; k++){
      for(int j = 0; j < thin; j++){

        if(is_true(any(k+1 == seqRefresh))){
          Rcout <<  "k = " << k + 1 << std::endl;
        }

        ind_m.fill(0); ind_l.fill(0);
        censInd.fill(0);
        cont_m = 0; cont_l = 0;

        for (int ee = 0; ee < n; ee++){
          if (y[ee]==0){
            fObsLap = Flaplace(arma::as_scalar(X.row(ee) * betaValue), tau, sigmaValue);
            probZero = 1-predicProb(betaZeroValue, link, X.row(ee));
            probCens = (probZero * fObsLap)/(1-probZero + probZero*fObsLap);
            if (runif(1)[0] < probCens){
              censInd[ee] = 1.0;
              matrizIndCens(k, ee) = 1.0;
              ind_m[cont_m] = ee;
              cont_m++;
            }
            else{
              ind_l[cont_l] = ee;
              matrizIndNCens(k, ee) = 1.0;
              cont_l++;
            }
          }
          else if (k==1){
            ind_n2[cont_n2] = ee;
            cont_n2++;
          }
        }

        m = sum(censInd);
        l = n1 - m;

        arma::uvec ind_v(n2 + m, arma::fill::zeros);

        if (m > 0) {
          arma::uvec ind_m2 = ind_m.subvec(0, cont_m - 1);
          ind_v.subvec(0, m - 1) = ind_m2;
          ind_v.subvec(m, n2+m-1) = ind_n2;
        }
        else {
          ind_v = ind_n2;
        }

        arma::colvec yS = y;
        arma::mat X2 = X.rows(ind_v);
        arma::colvec zSample2(n2+m, arma::fill::ones);

        aux = X * betaValue;
        for (int aa = 0; aa < n; aa++){
           if (censInd(aa)==1.0){
              meanModel = aux[aa] + theta*zSample(aa);
              sdModel = sqrt(psi2*sigmaValue*zSample(aa));
              yS[aa] = rnorm_trunc(meanModel, sdModel, InfPar[0], LimSup[0]);
           }
        }

        arma::colvec y2 = yS(ind_v);

        delta2 = diagvec((1/(psi2*sigmaValue)) * diagmat(yS - aux) * diagmat(yS - aux));

        gama2 = 2/sigmaValue + (theta*theta)/(psi2*sigmaValue);

        for(int o = 0; o < n; o++){
          if (varInd(o)==censInd(o)){
            delta2[o] = std::max(delta2[o], 1e-10);
            gama2 = std::max(gama2, 1e-10);
            zSample[o] = rgigRcpp(delta2[o], gama2, lambda);
          }
        }

        // Atualizando a cadeia de sigma.

        zSample2 = zSample(ind_v);
        aux2 = aux.rows(ind_v);
        termsSum = arma::as_scalar((y2 - aux2 - theta*zSample2).t() * diagmat(zSample2).i() * (y2 - aux2 - theta*zSample2));

        nTilde = n0 + 3*(n2+m);
        sTilde =  s0 + 2*sum(zSample2) + termsSum/psi2;

        sigmaValue = rinvgammaRcpp(nTilde/2,sTilde/2);

        // Atualizando beta do modelo quantílico
        diagU = diagmat(zSample2)*(psi2*sigmaValue);

        sigmaMinusOne = (X2.t() * diagU.i() * X2) + B0.i();

        sigmaMat = sigmaMinusOne.i();

        mu = sigmaMat * (B0.i()*b0 +
          (1/(psi2 * sigmaValue))*(X2.t() * diagmat(zSample2).i() * (y2 - theta*zSample2)));

        betaValue = mvrnormRcpp(mu, sigmaMat);

        // ********************************************* //
        // Estimação da parte da massa pontual do modelo //

        newBetaZeroValue = mvrnormRcpp(betaZeroValue, matIdsigmaBetaZero);

        postAnt = postDensMHstep(betaZeroValue, link, X, varInd - censInd, b0, B0);
        postNew = postDensMHstep(newBetaZeroValue, link, X, varInd - censInd, b0, B0);

        probAcceptance = std::min(0.0,postNew-postAnt);

        if(log(runif(1)[0])<probAcceptance){
          betaZeroValue = newBetaZeroValue;
          accept++;
        }

        // ********************************************* //

      }


      // Armazenando valores finais considerando parâmetro thin
      for(int jj = 0; jj < p; jj++){
        betaSample(k,jj) = betaValue[jj];
        betaZeroSample(k,jj) = betaZeroValue[jj];
      }

      sigmaSample[k] = sigmaValue;
   }

   acceptRate = accept/itNum;

   return List::create(
        Named("BetaSample") = betaSample,
        Named("BetaZeroSample") = betaZeroSample,
        Named("acceptRate") = acceptRate,
        Named("tau") = tau,
        Named("SigmaSample") = sigmaSample,
        Named("matrizIndCens") = matrizIndCens);
}
