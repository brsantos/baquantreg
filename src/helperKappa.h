double logPriorKappa (double value);

double logLikelihood (double value, arma::mat X, arma::vec y, arma::vec beta,
                      double sigma, arma::vec zSample, double tau, double psi2,
                      double theta, arma::vec spCoord1, arma::vec spCoord2,
                      double nugget);

double logLikelihood2 (double value, arma::mat X, arma::vec y,
                       arma::vec beta, double sigma, arma::vec zSample,
                       double tau, double psi2, double theta,
                       arma::vec spCoord1, arma::vec spCoord2,
                       arma::uvec indices, int m, double nugget);

double logLikelihood3 (double kappa, arma::mat X, arma::vec y, arma::vec beta,
                       double sigma, arma::vec zSample, double tau,
                       double psi2, double theta, arma::vec spCoord1,
                       arma::vec spCoord2, arma::uvec indices, int m,
                       double nugget);

double logLikelihood4 (double kappa, arma::mat X, arma::vec y,
                       arma::vec beta, double sigma, arma::vec zSample,
                       double tau, double psi2, double theta,
                       arma::vec spCoord1, arma::vec spCoord2,
                       double nugget);

double logPosterior (double value, arma::mat X, arma::vec y, arma::vec beta,
                     double sigma, arma::vec zSample, double tau, double psi2,
                     double theta, arma::vec spCoord1, arma::vec spCoord2,
                     double nugget);

double logPosterior2 (double value, arma::mat X, arma::vec y, arma::vec beta,
                      double sigma, arma::vec zSample, double tau, double psi2,
                      double theta, arma::vec spCoord1, arma::vec spCoord2,
                      arma::uvec indices, int m, double nugget);

double mhKappa(double value, arma::mat X, arma::vec y, arma::vec beta,
               double sigma, arma::vec zSample, double tau, double psi2,
               double theta, arma::vec spCoord1, arma::vec spCoord2,
               double tuneParam, double nugget);

double mhKappa2(double value, arma::mat X, arma::vec y, arma::vec beta,
                double sigma, arma::vec zSample, double tau, double psi2,
                double theta, arma::vec spCoord1, arma::vec spCoord2,
                double tuneParam, arma::uvec indices, int m, double nugget);
