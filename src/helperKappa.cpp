#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double logPriorKappa (double value){
  NumericVector aux(1); aux[0] = value;
  double output = log(dgamma(aux, 1.5, 2.0)[0]);
  return output;
}

// [[Rcpp::export]]
double logLikelihoodKappa (double kappa, arma::mat aux, arma::mat diagU,
                           arma::mat covMat, arma::mat covMatInv,
                           arma::vec spCoord1, arma::vec spCoord2,
                           double alpha, double jitter, bool newkappa){

  double output;
  double valDet;
  double sign;

  int n = aux.n_rows;

  if(newkappa){
    arma::mat covMat(n, n, arma::fill::zeros);

    for (int aa = 0; aa < n; aa++)
      for (int bb = 0; bb < n; bb++)
        covMat(aa,bb) = exp(-kappa * (pow(spCoord1(aa) - spCoord1(bb),2) +
          pow(spCoord2(aa) - spCoord2(bb),2)));

    arma::mat auxCov(n,n, arma::fill::zeros);
    auxCov.diag().fill(alpha+jitter);

    arma::mat R = chol((1-alpha)*covMat + auxCov).i();
    arma::mat covMatInv = R * R.t();

    arma::log_det(valDet, sign, covMat);

    output = -0.5*valDet -0.5*(as_scalar(aux.t()*diagU*covMatInv*diagU*aux));
  }
  else {
    arma::log_det(valDet, sign, (1-alpha)*covMat);

    output = -0.5 * valDet -
      0.5 * (as_scalar(aux.t() * diagU * covMatInv * diagU * aux));
  }

  return output;
}

double logLikelihoodKappa2 (double kappa, arma::mat aux, arma::mat diagU,
                            arma::mat covMat, arma::mat covMatInv,
                            arma::vec spCoord1, arma::vec spCoord2,
                            double alpha, double jitter, bool newkappa,
                            arma::uvec indices, int m){

  double output;
  double valDet;
  double sign;

  int n = aux.n_rows;

  if(newkappa){
    arma::mat covMat(n, n, arma::fill::zeros);
    arma::mat auxCov(m, m, arma::fill::zeros);

    for (int aa = 0; aa < n; aa++)
      for (int bb = 0; bb < n; bb++)
        covMat(aa,bb) = exp(-kappa * (pow(spCoord1(aa) - spCoord1(bb),2) +
          pow(spCoord2(aa) - spCoord2(bb),2)));

    arma::mat covMat2 = covMat.submat(indices, indices);
    arma::mat covMatAux = covMat.cols(indices);

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

    output = -0.5*valDet -0.5*(as_scalar(aux.t()*diagU*CovCov*diagU*aux));
  }
  else {
    arma::log_det(valDet, sign, (1-alpha)*covMat);

    output = -0.5 * valDet -
      0.5 * (as_scalar(aux.t() * diagU * covMatInv * diagU * aux));
  }

  return output;
}

// [[Rcpp::export]]
double mhKappa(double kappa, arma::vec spCoord1, arma::vec spCoord2,
               arma::mat aux, arma::mat diagU,
               arma::mat covMat, arma::mat covMatInv, double tuneParam,
               double alpha, double jitter){

  double postCur, postProp, kappaProp, densCur, densProp, new_kappa;

  kappaProp = rexp(1, tuneParam)[0];

  postCur = logLikelihoodKappa(kappa, aux, diagU, covMat, covMatInv,
                               spCoord1, spCoord2, alpha, jitter, false) +
                                 logPriorKappa(kappa);
  postProp = logLikelihoodKappa(kappaProp, aux, diagU, covMat, covMatInv,
                                   spCoord1, spCoord2, alpha, jitter, true) +
              logPriorKappa(kappaProp);

  NumericVector xx(1); xx[0] = kappa;
  NumericVector yy(1); yy[0] = kappaProp;

  densCur = dexp(xx, tuneParam)[0];
  densProp = dexp(yy, tuneParam)[0];

  double logAccepProb = postProp/postCur;

  if(runif(1)[0] < logAccepProb) new_kappa = kappaProp;
  else new_kappa = kappa;

  return new_kappa;
}

double mhKappa2(double kappa, arma::vec spCoord1, arma::vec spCoord2,
               arma::mat aux, arma::mat diagU,
               arma::mat covMat, arma::mat covMatInv, double tuneParam,
               double alpha, double jitter, arma::uvec indices, int m){

  double postCur, postProp, kappaProp, densCur, densProp, new_kappa;

  kappaProp = rexp(1, tuneParam)[0];

  postCur = logLikelihoodKappa2(kappa, aux, diagU, covMat, covMatInv,
                               spCoord1, spCoord2, alpha, jitter, false,
                               indices, m) +
                                 logPriorKappa(kappa);
  postProp = logLikelihoodKappa2(kappaProp, aux, diagU, covMat, covMatInv,
                                spCoord1, spCoord2, alpha, jitter, true,
                                indices, m) +
                                  logPriorKappa(kappaProp);

  NumericVector xx(1); xx[0] = kappa;
  NumericVector yy(1); yy[0] = kappaProp;

  densCur = dexp(xx, tuneParam)[0];
  densProp = dexp(yy, tuneParam)[0];

  double logAccepProb = exp(postProp-postCur)*densCur/densProp;

  if(runif(1)[0] < logAccepProb) new_kappa = kappaProp;
  else new_kappa = kappa;

  return new_kappa;
}