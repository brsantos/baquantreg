double postDensMHstep(arma::colvec parameters, int link, arma::mat X, arma::colvec varInd,
                      arma::colvec priorBeta, arma::mat priorSigma);

double predicProb(arma::colvec parameters, int link, arma::rowvec x);
