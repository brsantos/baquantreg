#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double logPosteriorV(arma::vec aux, double theta, double psi2, double sigma,
                     arma::vec vSample, double newV, arma::mat C, int indice){

  double output;
  arma::vec newzSample = vSample;
  newzSample(indice) = newV;

  arma::vec u = aux - theta*newzSample;
  double value = arma::as_scalar(u.t() * diagmat(1/sqrt(newzSample)) * C *
                                 diagmat(1/sqrt(newzSample)) * u);
  return - 0.5*log(newV) - 0.5*value/(psi2*sigma)  - newV/sigma;
}

NumericVector calcLogWeights(NumericVector weights){
  double a0 = max(weights);
  NumericVector output = exp(weights - a0)/sum(exp(weights - a0));
  return output;
}

double mtM(arma::vec aux, double theta, double psi2, double sigma,
           arma::vec vSample, double curV, int indice, arma::mat C,
           double tuneV, int k){

  double output;
  int nNeg, nNegRef;
  NumericVector prop(k), xRef(k), weights(k), weights2(k), weightsRef(k);
  prop = rexp(k, tuneV);
  for (int i = 0; i < k; i++){
    weights[i] = logPosteriorV(aux, theta, psi2, sigma, vSample, prop[i],
                               C, indice);
  }

  weights += log(dexp(prop, tuneV));
  weights2 = calcLogWeights(weights);

  double y = RcppArmadillo::sample(prop, 1, false, weights2)[0];
  for (int ii = 0; ii < k-1; ii++){
    xRef[ii] = rexp(1, tuneV)[0];
    weightsRef[ii] = logPosteriorV(aux, theta, psi2, sigma, vSample, xRef[ii],
                                   C, indice);
  }
  xRef[k-1] = curV;
  weightsRef[k-1] = logPosteriorV(aux, theta, psi2, sigma, vSample, xRef[k-1],
                                  C, indice);

  weightsRef += log(dexp(xRef, tuneV));

  double probA = (max(weights) + log(sum(exp(weights - max(weights))))) -
    (max(weightsRef) + sum(exp(weightsRef - max(weightsRef))));
  if (runif(1)[0] < exp(probA)) output = y;
  else output = curV;
  return output;
}

double logPriorKappa (double value){
  NumericVector aux(1); aux[0] = value;
  double output = log(dgamma(aux, 1, 0.5)[0]);
  return output;
}

NumericVector logPriorKappa2 (NumericVector value, double shape, double rate){
  return log(dgamma(value, shape, rate));
}

double logLikelihoodKappa (double kappa, arma::mat aux, arma::mat diagU,
                           arma::mat covMat, arma::mat covMatInv,
                           arma::mat matDist,
                           double alpha, double jitter, bool newkappa){

  double output;
  double valDet;
  double sign;

  int n = aux.n_rows;

  arma::mat auxCov(n,n, arma::fill::zeros);
  auxCov.diag().fill(alpha+jitter);

  if(newkappa){
    arma::mat covMatNew(n, n, arma::fill::zeros);

    covMatNew = exp(-kappa*matDist);

    arma::mat R = chol((1-alpha)*covMatNew + auxCov, "lower");
    arma::mat covMatInvNew = solve(trimatl(R),diagU*aux);

    arma::log_det(valDet, sign, (1-alpha)*covMatNew + auxCov);

    output = -0.5*valDet -0.5*(as_scalar(covMatInvNew.t()*covMatInvNew));
  }
  else {
    arma::log_det(valDet, sign, (1-alpha)*covMat + auxCov);
    output = -0.5*valDet -0.5*(as_scalar(aux.t()*diagU*covMatInv*diagU*aux));
  }
  return output;
}

double logLikelihoodKappa2 (double kappa, arma::mat aux, arma::mat diagU,
                            arma::mat covMat, arma::mat covMatInv,
                            arma::mat matDist,
                            double alpha, double jitter, bool newkappa,
                            arma::uvec indices, int m){

  double output;
  double valDet;
  double sign;

  int n = aux.n_rows;

  if(newkappa){

    arma::mat auxCov(m, m, arma::fill::zeros);
    auxCov.diag().fill(alpha + jitter);
    arma::mat covMatNew(n, n, arma::fill::zeros);

    covMatNew = exp(-kappa*matDist);

    arma::mat covMat2 = covMatNew.submat(indices, indices);
    arma::mat covMatAux = covMatNew.cols(indices);

    arma::mat cholCov = arma::chol((1-alpha)*covMat2 + auxCov, "lower");
    arma::mat covMatInvNew =  solve(trimatl(cholCov), covMatAux.t());
    arma::mat matAux = covMatInvNew.t()*covMatInvNew;
    arma::mat sigmaDot = diagmat(1/(alpha +
      (1-alpha)*(covMatNew.diag() - matAux.diag())));

    arma::mat cholCov2 = (1-alpha)*covMat2 + covMatAux.t()*sigmaDot*covMatAux;

    arma::mat matM = arma::chol(cholCov2 + auxCov, "lower");
    arma::mat matM2 = solve(trimatl(matM),covMatAux.t());
    arma::mat matM3 = matM2.t() * matM2;

    arma::mat CovCov = sigmaDot - sigmaDot * matM3 * sigmaDot;

    arma::log_det(valDet, sign, matAux + diagmat(1/sigmaDot.diag()));
    output = -0.5*valDet -0.5*(as_scalar(aux.t()*diagU*CovCov*diagU*aux));
  }
  else {
    arma::log_det(valDet, sign, covMat);

    output = -0.5*valDet -0.5*(as_scalar(aux.t()*diagU*covMatInv*diagU*aux));
  }

  return output;
}


double discKappa2(NumericVector lambda, NumericVector lambdaPrior,
                  arma::mat matDist,
                  arma::mat aux, arma::mat diagU,
                  arma::mat covMat, arma::mat covMatInv,
                  double alpha, double jitter, arma::uvec indices, int m){

  double postCur, postProp, kappaProp, densCur, densProp, new_kappa;

  int sizeLambda = lambda.size();

  NumericVector logPost(sizeLambda);

  for(int k = 0; k < sizeLambda; k++){
    logPost[k] = logLikelihoodKappa2(lambda[k], aux, diagU, covMat, covMatInv,
                                     matDist, alpha, jitter, true,
                                     indices, m) + lambdaPrior[k];
  }

  NumericVector weights = calcLogWeights(logPost);

  double output = RcppArmadillo::sample(lambda, 1, false, weights)[0];

  return output;
}