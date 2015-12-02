using namespace Rcpp;
using namespace arma;

double rinvgammaRcpp(double shape, double scale);

colvec mvrnormRcpp(colvec mu, mat Sigma);

double Flaplace (double predictor, double tau, double sigma);
