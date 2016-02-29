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

  densCur = dexp(xx, kappaProp)[0];
  densProp = dexp(yy, kappa)[0];

  double logAccepProb = exp(postProp-postCur)*densCur/densProp;

  if(runif(1)[0] < logAccepProb) new_kappa = kappaProp;
  else new_kappa = kappa;

  return new_kappa;
}

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

// double logLikelihood2Kappa (double kappa, arma::mat aux, arma::mat diagU,
//                        arma::mat CovCov){
//
//   double output;
//
//   double valDetC, signDetC, valDetM, signDetM;
//
//   log_det(valDetC, signDetC, CovCov);
//   log_det(valDetM, signDetM, matM);
//
//   output = -0.5 * ((-as_scalar(sum(log(sigmaDot))) - valDetC + valDetM)) -
//     0.5*(as_scalar(aux.t() * diagU.t() * CovCov * diagU *aux));
//
//   return output;
// }


// double logPosterior2 (double value, arma::mat X, arma::vec y, arma::vec beta,
//                       double sigma, arma::vec zSample, double tau, double psi2,
//                       double theta, arma::vec spCoord1, arma::vec spCoord2,
//                       arma::uvec indices, int m, double alpha){
//
//   double output = logLikelihood2(value, X, y, beta, sigma, zSample, tau, psi2,
//                                  theta, spCoord1, spCoord2, indices, m,
//                                  alpha) + logPriorKappa(value);
//   return output;
// }

// double mhKappa2(double value, arma::mat X, arma::vec y, arma::vec beta,
//                 double sigma, arma::vec zSample, double tau, double psi2,
//                 double theta, arma::vec spCoord1, arma::vec spCoord2,
//                 double tuneParam, arma::uvec indices, int m, double alpha){
//
//   double postCur, postProp, kappaProposal, densCur, densProp, new_kappa;
//
//   kappaProposal = rexp(1, tuneParam)[0];
//
//   postCur = logPosterior2(value, X, y, beta, sigma, zSample, tau, psi2,
//                           theta, spCoord1, spCoord2, indices, m, alpha);
//   postProp = logPosterior2(kappaProposal, X, y, beta, sigma, zSample,
//                            tau, psi2, theta, spCoord1, spCoord2, indices,
//                            m, alpha);
//
//   NumericVector xx(1); xx[0] = value;
//   NumericVector yy(1); yy[0] = kappaProposal;
//
//   densCur = dexp(xx, kappaProposal)[0];
//   densProp = dexp(yy, value)[0];
//
//   double accepProb = std::min(1.0, (postProp*densCur)/(postCur*densProp));
//
//   if(runif(1)[0] < accepProb) new_kappa = kappaProposal;
//   else new_kappa = value;
//
//   return new_kappa;
// }


// double mhAlpha2(double kappa, arma::mat X, arma::vec y, arma::vec beta,
//                 double sigma, arma::vec zSample, double tau,
//                 double psi2, double theta, arma::vec spCoord1,
//                 arma::vec spCoord2, double value,
//                 arma::uvec indices, int m, double tuneA){
//
//   double postCur, postProp, alphaProp, new_alpha;
//
//   double p1, p2;
//   double precision = tuneA;
//   p1 = value * precision;
//   p2 = (1-value) * precision;
//   alphaProp = rbeta(1, p1, p2)[0];
//
//   postCur = logLikelihood2(kappa, X, y, beta, sigma, zSample, tau, psi2,
//                            theta, spCoord1, spCoord2, indices, m, value);
//   postProp = logLikelihood2(kappa, X, y, beta, sigma, zSample, tau, psi2,
//                             theta, spCoord1, spCoord2, indices, m, alphaProp);
//
//   double accepProb = std::min(1.0, (postProp)/(postCur));
//
//   if(runif(1)[0] < accepProb) new_alpha = alphaProp;
//   else new_alpha = value;
//
//   return new_alpha;
// }