double logPriorKappa (double value);

double logLikelihoodKappa (double kappa, arma::mat aux, arma::mat diagU,
                           arma::mat covMat, arma::mat covMatInv,
                           arma::vec spCoord1, arma::vec spCoord2,
                           double alpha, double jitter, bool newkappa);

double logLikelihoodKappa2 (double kappa, arma::mat aux, arma::mat diagU,
                            arma::mat covMat, arma::mat covMatInv,
                            arma::vec spCoord1, arma::vec spCoord2,
                            double alpha, double jitter, bool newkappa,
                            arma::uvec indices, int m);

double mhKappa(double kappa, arma::vec spCoord1, arma::vec spCoord2,
               arma::mat aux, arma::mat diagU,
               arma::mat covMat, arma::mat covMatInv, double tuneParam,
               double alpha, double jitter);

double mhKappa2(double kappa, arma::vec spCoord1, arma::vec spCoord2,
               arma::mat aux, arma::mat diagU,
               arma::mat covMat, arma::mat covMatInv, double tuneParam,
               double alpha, double jitter, arma::uvec indices, int m);

