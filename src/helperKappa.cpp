#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double logPriorKappa (double value){
  NumericVector aux(1); aux[0] = value;
  double output = log(dgamma(aux, 1.5, 2.0)[0]);
  return output;
}

double logLikelihood (double value, arma::mat X, arma::vec y, arma::vec beta,
                      double sigma, arma::vec zSample, double tau, double psi2,
                      double theta, arma::vec spCoord1, arma::vec spCoord2,
                      double nugget){

  double output;

  arma::mat aux = y - X*beta - theta*zSample;
  int n = X.n_rows;
  arma::mat covMat(n, n, arma::fill::zeros);

  for (int aa = 0; aa < n; aa++)
    for (int bb = 0; bb < n; bb++)
      covMat(aa,bb) = exp(-value * (pow(spCoord1(aa) - spCoord1(bb),2) +
        pow(spCoord2(aa) - spCoord2(bb),2)));

  arma::mat diagU = diagmat(sqrt(zSample))*(sqrt(psi2));

  arma::mat auxCov(n,n, arma::fill::zeros);

  arma::mat R = chol(covMat + nugget * auxCov.eye()).i();

  arma::mat covMatInv = R * R.t();

  output = -0.5 * (as_scalar(aux.t() * diagU.i().t() * covMatInv *
    diagU.i() *aux));

  return output;
}

double logLikelihood2 (double value, arma::mat X, arma::vec y,
                       arma::vec beta, double sigma, arma::vec zSample,
                       double tau, double psi2, double theta,
                       arma::vec spCoord1, arma::vec spCoord2,
                       arma::uvec indices, int m, double nugget){

  double output;

  arma::mat aux = y - X*beta - theta*zSample;
  int n = X.n_rows;
  arma::mat covMat(n, n, arma::fill::zeros), covMat2(m, m, arma::fill::zeros),
    covMatAux(n, m, arma::fill::zeros);

  for(int aa = 0; aa < n; aa++)
    for(int bb = 0; bb < n; bb++)
      covMat(aa,bb) = exp(-value * (pow(spCoord1(aa) - spCoord1(bb),2) +
        pow(spCoord2(aa) - spCoord2(bb),2)));

  arma::mat diagU = diagmat(sqrt(1/zSample)) * (sqrt(1/psi2));

  arma::mat auxCov(m, m, arma::fill::zeros);

  // Selecting the terms of the knots.
  covMat2 = covMat.submat(indices, indices);

  covMatAux = covMat.cols(indices);

  arma::mat R = chol(covMat2 + nugget * auxCov.eye()).i();

  arma::mat covMatInv = R * R.t();

  output = -0.5*(as_scalar(aux.t() * diagU.t() * covMatAux *
    covMatInv * covMatAux.t() * diagU *aux));

  return output;
}

double logLikelihood3 (double kappa, arma::mat X, arma::vec y, arma::vec beta,
                       double sigma, arma::vec zSample, double tau,
                       double psi2, double theta, arma::vec spCoord1,
                       arma::vec spCoord2, arma::uvec indices, int m,
                       double nugget){

  double output;

  arma::mat aux = y - X * beta - theta * zSample;
  int n = X.n_rows;
  arma::mat covMat(n, n, arma::fill::zeros), covMat2(m, m, arma::fill::zeros),
    covMatAux(n, m, arma::fill::zeros);

  for(int aa = 0; aa < n; aa++)
    for(int bb = 0; bb < n; bb++)
      covMat(aa,bb) = exp(-kappa * (pow(spCoord1(aa) - spCoord1(bb),2) +
        pow(spCoord2(aa) - spCoord2(bb),2)));

  arma::mat diagU = diagmat(sqrt(1/zSample)) * (sqrt(1/psi2));

  arma::mat auxCov(m,m, arma::fill::zeros);

  // Selecting the terms of the knots.
  covMat2 = covMat.submat(indices, indices);

  covMatAux = covMat.cols(indices);

  arma::mat R = chol(covMat2 + nugget * auxCov.eye()).i();

  arma::mat covMatInv = R * R.t();

  output = -0.5*(as_scalar(aux.t()* diagU.t() * covMatAux * covMatInv *
    covMatAux.t() * diagU *aux));

  return output;
}

double logLikelihood4 (double kappa, arma::mat X, arma::vec y,
                       arma::vec beta, double sigma, arma::vec zSample,
                       double tau, double psi2, double theta,
                       arma::vec spCoord1, arma::vec spCoord2,
                       double nugget){

  double output;

  arma::mat aux = y - X * beta - theta * zSample;
  int n = X.n_rows;
  arma::mat covMat(n, n, arma::fill::zeros), auxCov(n, n, arma::fill::zeros);

  for(int aa = 0; aa < n; aa++)
    for(int bb = 0; bb < n; bb++)
      covMat(aa,bb) = exp(-kappa * (pow(spCoord1(aa) - spCoord1(bb),2) +
        pow(spCoord2(aa) - spCoord2(bb),2)));

  arma::mat diagU = diagmat(sqrt(1/zSample)) * (sqrt(1/psi2));

  arma::mat R = chol(covMat + nugget * auxCov.eye()).i();

  arma::mat covMatInv = R * R.t();

  output = -0.5*(as_scalar(aux.t() * diagU.t() * covMatInv * diagU *aux));

  return output;
}

double logPosterior (double value, arma::mat X, arma::vec y, arma::vec beta,
                     double sigma, arma::vec zSample, double tau, double psi2,
                     double theta, arma::vec spCoord1, arma::vec spCoord2,
                     double nugget){

  double output = logLikelihood(value, X, y, beta, sigma, zSample, tau,
                                psi2, theta, spCoord1, spCoord2, nugget) +
                                  logPriorKappa(value);
  return output;
}

double logPosterior2 (double value, arma::mat X, arma::vec y, arma::vec beta,
                      double sigma, arma::vec zSample, double tau, double psi2,
                      double theta, arma::vec spCoord1, arma::vec spCoord2,
                      arma::uvec indices, int m, double nugget){

  double output = logLikelihood2(value, X, y, beta, sigma, zSample, tau, psi2,
                                 theta, spCoord1, spCoord2, indices, m,
                                 nugget) + logPriorKappa(value);
  return output;
}

double mhKappa(double value, arma::mat X, arma::vec y, arma::vec beta,
               double sigma, arma::vec zSample, double tau, double psi2,
               double theta, arma::vec spCoord1, arma::vec spCoord2,
               double tuneParam, double nugget){

  double postCur, postProp, kappaProposal, densCur, densProp, new_kappa;

  kappaProposal = rexp(1, tuneParam)[0];

  postCur = logPosterior(value, X, y, beta, sigma, zSample, tau, psi2,
                         theta, spCoord1, spCoord2, nugget);
  postProp = logPosterior(kappaProposal, X, y, beta, sigma, zSample, tau,
                          psi2, theta, spCoord1, spCoord2, nugget);

  NumericVector xx(1); xx[0] = value;
  NumericVector yy(1); yy[0] = kappaProposal;

  densCur = dexp(xx, kappaProposal)[0];
  densProp = dexp(yy, value)[0];

  double accepProb = std::min(1.0, (postProp*densCur)/(postCur * densProp));

  if(runif(1)[0] < accepProb) new_kappa = kappaProposal;
  else new_kappa = value;

  return new_kappa;
}

double mhKappa2(double value, arma::mat X, arma::vec y, arma::vec beta,
                double sigma, arma::vec zSample, double tau, double psi2,
                double theta, arma::vec spCoord1, arma::vec spCoord2,
                double tuneParam, arma::uvec indices, int m, double nugget){

  double postCur, postProp, kappaProposal, densCur, densProp, new_kappa;

  kappaProposal = rexp(1, tuneParam)[0];

  postCur = logPosterior2(value, X, y, beta, sigma, zSample, tau, psi2,
                          theta, spCoord1, spCoord2, indices, m, nugget);
  postProp = logPosterior2(kappaProposal, X, y, beta, sigma, zSample,
                           tau, psi2, theta, spCoord1, spCoord2, indices,
                           m, nugget);

  NumericVector xx(1); xx[0] = value;
  NumericVector yy(1); yy[0] = kappaProposal;

  densCur = dexp(xx, kappaProposal)[0];
  densProp = dexp(yy, value)[0];

  double accepProb = std::min(1.0, (postProp*densCur)/(postCur*densProp));

  if(runif(1)[0] < accepProb) new_kappa = kappaProposal;
  else new_kappa = value;

  return new_kappa;
}