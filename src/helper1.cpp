#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

double postDensMHstep(arma::colvec parameters, int link, arma::mat X, arma::colvec varInd,
                      arma::colvec priorBeta, arma::mat priorSigma){

  int n = X.n_rows;

  double priorWeight, likelihoodInf, postDens;
  arma::colvec valLik(n), aux(n);
  Rcpp::NumericVector aux1(1);

  priorWeight = -0.5 * as_scalar((parameters.t() - priorBeta.t()) * priorSigma.i() * (parameters - priorBeta));

  aux = X * parameters;

  if (link==1){
    for (int j=0; j < n; j++){
      if(varInd(j)==1) valLik[j] = log(exp(aux[j])/(1 + exp(aux[j])));
      else valLik[j] = -log(1 + exp(aux[j]));
    }
  }
  else{
    for (int j=0; j < n; j++){
      aux1[0] = as_scalar(X.row(j) * parameters);
      if(varInd(j)==1) valLik[j] = log(pnorm(aux1, 0.0, 1.0)[0]);
      else valLik[j] = log(1-pnorm(aux1, 0.0, 1.0)[0]);
    }
  }

  likelihoodInf = sum(valLik);
  postDens = likelihoodInf + priorWeight;
  return postDens;
}

double predicProb(arma::colvec parameters, int link, arma::rowvec x){
  double prob;
  Rcpp::NumericVector aux(1);
  if (link==1) prob = exp(as_scalar(x * parameters))/(1 + exp(as_scalar(x * parameters)));
  else{
    aux[0] = as_scalar(x * parameters);
    prob = pnorm(aux, 0.0, 1.0)[0];
  }
  return prob;
}