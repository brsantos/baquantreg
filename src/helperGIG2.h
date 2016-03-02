double logPosteriorV(double v, double delta2, double gama2);

double mtM(arma::vec aux, double theta, double psi2, double sigma,
                arma::vec vSample, double curV, int indice, arma::mat C,
                double tuneV, int k);
