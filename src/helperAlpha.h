double logLikelihoodAlpha (double alpha, arma::mat aux, arma::mat diagU,
                           arma::mat covMat, double jitter);


double logLikelihoodAlpha2 (double alpha, arma::mat aux, arma::mat diagU,
                            arma::mat covMat, arma::mat covMat2,
                            arma::mat covMatAux, double jitter,
                            arma::uvec indices, int m);

double mhAlpha(double alpha, arma::mat aux, arma::mat diagU,
               arma::mat covMat, double tuneA, double jitter);

double mhAlpha2(double alpha, arma::mat aux, arma::mat diagU,
                arma::mat covMat, arma::mat covMat2,
                arma::mat covMatAux, double tuneA, double jitter,
                arma::uvec indices, int m);
