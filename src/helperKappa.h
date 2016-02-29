double logPriorKappa (double value);

double logLikelihoodKappa (double kappa, arma::mat aux, arma::mat diagU,
                           arma::mat covMat, arma::mat covMatInv);

double logLikelihoodKappaNew (double kappa, arma::mat aux, arma::mat diagU,
                              arma::vec spCoord1, arma::vec spCoord2,
                              double alpha, double jitter);

double mhKappa(double kappa, arma::vec spCoord1, arma::vec spCoord2,
               arma::mat aux, arma::mat diagU,
               arma::mat covMat, arma::mat covMatInv, double tuneParam,
               double alpha, double jitter);

double logLikelihoodAlpha (double alpha, arma::mat aux, arma::mat diagU,
                           arma::mat covMat, double jitter);

double mhAlpha(double alpha, arma::mat aux, arma::mat diagU,
               arma::mat covMat, double tuneA, double jitter);


