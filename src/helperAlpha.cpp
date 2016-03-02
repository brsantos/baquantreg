#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double logLikelihoodAlpha (double alpha, arma::mat aux, arma::mat diagU,
                           arma::mat covMat, double jitter){

  double output;
  double valDet;
  double sign;

  int n = covMat.n_rows;
  arma::mat covMatAux(n,n, arma::fill::zeros);
  covMatAux.diag().fill(alpha + jitter);

  arma::mat cholCov = chol((1-alpha)*covMat + covMatAux).i();
  arma::mat covMatInv = cholCov * cholCov.t();

  arma::log_det(valDet, sign, (1-alpha)*covMat);

  output = -0.5*valDet-0.5*(as_scalar(aux.t()*diagU*covMatInv*diagU*aux));

  return output;
}

// [[Rcpp::export]]
double logLikelihoodAlpha2 (double alpha, arma::mat aux, arma::mat diagU,
                            arma::mat covMat, arma::mat covMat2,
                            arma::mat covMatAux, double jitter,
                            arma::uvec indices, int m){

  double output;
  double valDet;
  double sign;

  arma::mat auxCov(m, m, arma::fill::zeros);

  auxCov.diag().fill(alpha + jitter);

  arma::mat cholCov = arma::chol((1-alpha)*covMat2 + auxCov).i();
  arma::mat covMatInv =  cholCov * cholCov.t();
  arma::mat matAux = covMatAux*covMatInv*covMatAux.t();
  arma::mat sigmaDot = diagmat(alpha +
    (1-alpha)*(covMat.diag() - matAux.diag())).i();

  arma::mat cholCov2 = (1-alpha)*covMat2 +
    covMatAux.t()*sigmaDot*covMatAux;

  arma::mat matM = arma::chol(cholCov2 + auxCov).i();
  arma::mat matM2 = matM * matM.t();

  arma::mat CovCov = sigmaDot -
    sigmaDot * covMatAux * matM2 * covMatAux.t() * sigmaDot;

  arma::log_det(valDet, sign, matAux + sigmaDot.i());

  output = -0.5*valDet-0.5*(as_scalar(aux.t()*diagU*covMatInv*diagU*aux));

  return output;
}

// [[Rcpp::export]]
double mhAlpha(double alpha, arma::mat aux, arma::mat diagU,
               arma::mat covMat, double tuneA, double jitter){

  double postCur, postProp, alphaProp, new_alpha;

  double p1, p2;
  double precision = tuneA;
  p1 = alpha * precision;
  p2 = (1-alpha) * precision;
  // alphaProp = rbeta(1, p1, p2)[0];
  alphaProp = runif(1)[0]*0.8 + 0.1;

  postCur = logLikelihoodAlpha(alpha, aux, diagU, covMat, jitter);
  postProp = logLikelihoodAlpha(alphaProp, aux, diagU, covMat, jitter);

  double logAccepProb = postProp - postCur;

  if(log(runif(1)[0]) < logAccepProb) new_alpha = alphaProp;
  else new_alpha = alpha;

  return new_alpha;
}

// [[Rcpp::export]]
double mhAlpha2(double alpha, arma::mat aux, arma::mat diagU,
                arma::mat covMat, arma::mat covMat2,
                arma::mat covMatAux, double tuneA, double jitter,
                arma::uvec indices, int m){

  double postCur, postProp, alphaProp, new_alpha;

  double p1, p2;
  double precision = tuneA;
  p1 = alpha * precision;
  p2 = (1-alpha) * precision;
  // alphaProp = rbeta(1, p1, p2)[0];
  alphaProp = runif(1)[0]*0.8 + 0.1;

  postCur = logLikelihoodAlpha2(alpha, aux, diagU, covMat, covMat2, covMatAux,
                                jitter, indices, m);
  postProp = logLikelihoodAlpha2(alphaProp, aux, diagU, covMat, covMat2,
                                 covMatAux, jitter, indices, m);

  double logAccepProb = postProp/postCur;

  if(runif(1)[0] < logAccepProb) new_alpha = alphaProp;
  else new_alpha = alpha;

  return new_alpha;
}