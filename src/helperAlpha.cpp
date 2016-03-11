#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double logLikelihoodAlpha (double alpha, arma::mat aux, arma::mat diagU,
                           arma::mat covMat, double jitter){

  double output;
  double valDet;
  double sign;

  int n = covMat.n_rows;
  arma::mat covMatAux(n,n, arma::fill::zeros);
  covMatAux.diag().fill(alpha + jitter);

  arma::mat cholCov = chol((1-alpha)*covMat + covMatAux, "lower");
  arma::mat covMatInv = solve(trimatl(cholCov),diagU*aux);

  arma::log_det(valDet, sign, (1-alpha)*covMat+covMatAux);

  output = -0.5*valDet-0.5*(as_scalar(covMatInv.t()*covMatInv));

  return output;
}

double logLikelihoodAlpha2 (double alpha, arma::mat aux, arma::mat diagU,
                            arma::mat covMat, arma::mat covMat2,
                            arma::mat covMatAux, double jitter,
                            arma::uvec indices, int m){

  double output;
  double valDet;
  double sign;

  arma::mat auxCov(m, m, arma::fill::zeros);
  auxCov.diag().fill(alpha + jitter);

  arma::mat cholCov = arma::chol((1-alpha)*covMat2 + auxCov, "lower");
  arma::mat covMatInv =  solve(trimatl(cholCov),covMatAux.t());
  arma::mat matAux = covMatInv.t()*covMatInv;
  arma::mat sigmaDot = diagmat(1/(alpha +
    (1-alpha)*(covMat.diag() - matAux.diag())));

  arma::mat cholCov2 = (1-alpha)*covMat2 + covMatAux.t()*sigmaDot*covMatAux;

  arma::mat matM = arma::chol(cholCov2 + auxCov, "lower");
  arma::mat matM2 = solve(trimatl(matM), covMatAux.t());
  arma::mat matM3 = matM2.t() * matM2;

  arma::mat CovCov = sigmaDot - sigmaDot * matM3 * sigmaDot;

  arma::log_det(valDet, sign, matAux + diagmat(1/sigmaDot.diag()));

  output = -0.5*valDet-0.5*(as_scalar(aux.t()*diagU*CovCov*diagU*aux));

  return output;
}

double mhAlpha(double alpha, arma::mat aux, arma::mat diagU,
               arma::mat covMat, double tuneA, double jitter){

  double postCur, postProp, alphaProp, new_alpha;

  double p1, p2;
  double precision = tuneA;
  p1 = alpha * precision;
  p2 = (1-alpha) * precision;
  // alphaProp = rbeta(1, p1, p2)[0];
  alphaProp = runif(1)[0];

  postCur = logLikelihoodAlpha(alpha, aux, diagU, covMat, jitter);
  postProp = logLikelihoodAlpha(alphaProp, aux, diagU, covMat, jitter);

  double logAccepProb = postProp - postCur;

  if(log(runif(1)[0]) < logAccepProb) new_alpha = alphaProp;
  else new_alpha = alpha;

  return new_alpha;
}

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
  alphaProp = runif(1)[0];

  postCur = logLikelihoodAlpha2(alpha, aux, diagU, covMat, covMat2, covMatAux,
                                jitter, indices, m);
  postProp = logLikelihoodAlpha2(alphaProp, aux, diagU, covMat, covMat2,
                                 covMatAux, jitter, indices, m);

  double logAccepProb = postProp - postCur;

  if(log(runif(1)[0]) < logAccepProb) new_alpha = alphaProp;
  else new_alpha = alpha;

  return new_alpha;
}
