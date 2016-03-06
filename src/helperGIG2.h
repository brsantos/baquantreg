double logPosteriorV(double v, double delta2, double gama2);

double mtM(arma::vec aux, double theta, double psi2, double sigma,
                arma::vec vSample, double curV, int indice, arma::mat C,
                double tuneV, int k);

double logPosteriorV2(arma::vec aux, double theta, double psi2, double sigma,
                      arma::vec vSample, arma::mat C);

Rcpp::NumericVector mtM2(arma::vec aux, double theta, double psi2, double sigma,
                   arma::vec vSample, int n, arma::mat C, double tuneV, int k);