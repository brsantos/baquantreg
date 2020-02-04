#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat checkPoints_4d(arma::colvec gridy1, arma::colvec gridy2,
                           arma::colvec gridy3, arma::colvec gridy4,
                           arma::mat directions,
                           arma::mat orthBasis1, arma::mat orthBasis2,
                           arma::mat orthBasis3,
                           arma::mat fullestimates, arma::colvec xvalue,
                           bool splines, arma::colvec addterm){

  int n1 = gridy1.size();
  int n2 = gridy2.size();
  int n3 = gridy3.size();
  int n4 = gridy4.size();

  int ndirections = directions.n_rows;
  int p = fullestimates.n_rows;

  arma::mat effDir1 = fullestimates.row(p - 3);
  arma::mat effDir2 = fullestimates.row(p - 2);
  arma::mat effDir3 = fullestimates.row(p - 1);
  arma::mat estimates;
  if(p == 2) estimates = fullestimates.row(0);
  else estimates = fullestimates.rows(0, p - 4);

  arma::mat checking;
  int number_points = 0;

  for(int k = 0; k < n1; k++){
    for(int j = 0; j < n2; j++){
      for(int l = 0; l < n3; l++){
        for(int m = 0; m < n4; m++){

          int partcontour = 0;
          int count = 0;
          while (partcontour == 0 && count < ndirections){
            double lefthand = directions(count, 0) * gridy1(k) +
              directions(count, 1) * gridy2(j) +
              directions(count, 2) * gridy3(l) +
              directions(count, 3) * gridy4(m);
            double righthand = effDir1(0, count) *
              (orthBasis1(count, 0) * gridy1(k) +
              orthBasis1(count, 1) * gridy2(j) +
              orthBasis1(count, 2) * gridy3(l) +
              orthBasis1(count, 2) * gridy4(m)) +
              effDir2(0, count) *
              (orthBasis2(count, 0) * gridy1(k) +
              orthBasis2(count, 1) * gridy2(j) +
              orthBasis2(count, 2) * gridy3(l) +
              orthBasis2(count, 2) * gridy4(m)) +
              effDir3(0, count) *
              (orthBasis3(count, 0) * gridy1(k) +
              orthBasis3(count, 1) * gridy2(j) +
              orthBasis3(count, 2) * gridy3(l) +
              orthBasis3(count, 2) * gridy4(m)) +
              arma::as_scalar(xvalue.t() * estimates.col(count));
            if (splines) righthand = righthand + addterm(count);
            if (lefthand < righthand) partcontour = 1;
            else if (count == ndirections - 1){
              arma::rowvec coordinates =
                {gridy1(k), gridy2(j), gridy3(l), gridy4(m)};
              checking.insert_rows(number_points, coordinates);
              number_points++;
            }
            count++;
          }
        }
      }
    }
  }

  return checking;
}
