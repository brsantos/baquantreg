double logPosteriorV(arma::vec aux, double theta, double psi2, double sigma,
                     arma::vec vSample, double newV, arma::mat C, int indice);

double mtM(arma::vec aux, double theta, double psi2, double sigma,
           arma::vec vSample, double curV, int indice, arma::mat C,
           double tuneV, int k);


double logPriorKappa (double value);

Rcpp::NumericVector logPriorKappa2 (Rcpp::NumericVector value, double shape, double rate);

double logLikelihoodKappa (double kappa, arma::mat aux, arma::mat diagU,
                           arma::mat covMat, arma::mat covMatInv,
                           arma::mat matDist,
                           double alpha, double jitter, bool newkappa);

double logLikelihoodKappa2 (double kappa, arma::mat aux, arma::mat diagU,
                            arma::mat covMat, arma::mat covMatInv,
                            arma::mat matDist,
                            double alpha, double jitter, bool newkappa,
                            arma::uvec indices, int m);

double discKappa2(Rcpp::NumericVector lambda, Rcpp::NumericVector lambdaPrior,
                  arma::mat matDist,
                  arma::mat aux, arma::mat diagU,
                  arma::mat covMat, arma::mat covMatInv,
                  double alpha, double jitter, arma::uvec indices, int m);
