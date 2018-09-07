#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat checkPoints(arma::colvec gridy1, arma::colvec gridy2, arma::vec effDir,
                                     arma::mat directions, arma::mat orthBasis,
                                     arma::mat estimates, arma::colvec xvalue) {

  int n1 = gridy1.size();
  int n2 = gridy2.size();

  int ndirections = directions.n_rows;

  arma::mat checking = arma::zeros(n1, n2);

  Rcout << "The number of directions is = " << ndirections  <<  std::endl;

  for(int k = 0; k < n1; k++){
    Rcout << "Max count = " << std::endl;
    for(int j = 0; j < n2; j++){

      int partcontour = 0;
      int count = 0;
      while (partcontour == 0 && count < ndirections){
        double lefthand = directions(count, 0) * gridy1(k) + directions(count, 1) * gridy2(j);
        double righthand = effDir(count) * (orthBasis(count, 0) * gridy1(k) + orthBasis(count, 1) * gridy2(j)) +
          arma::as_scalar(xvalue.t() * estimates.col(count));
        if (lefthand < righthand) partcontour = 1;
        else if (count == ndirections - 1){
          Rcout << "Got here" << std::endl;
          checking(k, j) = 1;
        }
        // Rcout << "count = " << count << std::endl;
        count++;
      }
      Rcout << count << " ";
    }
    Rcout << std::endl;
  }

  return checking;
}

/*** R
directionPoint <- 8
angles <- (0:(directionPoint-1))*2*pi/directionPoint
vectorDir <- cbind(cos(angles), sin(angles))

orthBasis <- t(apply(vectorDir, 1, function(a){
  u_1 <- c(1,0)

  A <- cbind(a, u_1)
  x.qr <- qr.Q(qr(A))[, 2]
}))

xvalue <- c(1, 1, 0, 1)
beta <- matrix(rnorm(length(xvalue)*length(angles), sd = 0.5), ncol = length(angles))

effDir <- rnorm(directionPoint)

checkPoints(seq(-1, 1, length = 10), seq(-1, 1, length = 10), effDir,
            vectorDir, orthBasis, beta, xvalue)
*/
