#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat checkPoints(arma::colvec gridy1, arma::colvec gridy2,
                      arma::mat directions, arma::mat orthBasis,
                      arma::mat fullestimates, arma::colvec xvalue,
                      bool splines, arma::colvec addterm){

  int n1 = gridy1.size();
  int n2 = gridy2.size();

  int ndirections = directions.n_rows;
  int p = fullestimates.n_rows;

  arma::mat effDir = fullestimates.row(p - 1);
  arma::mat estimates;
  if(p == 2) estimates  = fullestimates.row(0);
  else estimates = fullestimates.rows(0, p - 2);

  arma::mat checking = arma::zeros(n1, n2);

  for(int k = 0; k < n1; k++){
    for(int j = 0; j < n2; j++){

      int partcontour = 0;
      int count = 0;
      while (partcontour == 0 && count < ndirections){
        double lefthand = directions(count, 0) * gridy1(k) +
          directions(count, 1) * gridy2(j);
        double righthand = effDir(0, count) * (orthBasis(count, 0) * gridy1(k) +
                                  orthBasis(count, 1) * gridy2(j)) +
                                  arma::as_scalar(xvalue.t() *
                                  estimates.col(count));
        if (splines) righthand = righthand + addterm(count);
        if (lefthand < righthand) partcontour = 1;
        else if (count == ndirections - 1){
          checking(k, j) = 1;
        }
        count++;
      }
    }
  }

  return checking;
}
