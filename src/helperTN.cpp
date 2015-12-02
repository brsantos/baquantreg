#include <Rcpp.h>
using namespace Rcpp;

double rtnormRcpp(NumericVector a, NumericVector b, double mean, double sigma){
  double out(1);
  NumericVector unif(1);
  NumericVector quant(1);

  unif = runif(1);

  quant(0) = pnorm(a,mean, sigma)[0] + unif[0]*(pnorm(b,mean,sigma)[0] - pnorm(a,mean,sigma)[0]);
  out = mean + sigma*qnorm(quant, 0.0, 1.0)[0];

  return(out);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to be in the interval
// (a,b) via rejection sampling.
// ======================================================================

// [[Rcpp::export]]

double norm_rs(double a, double b)
{
  double x;
  x = rnorm(1, 0.0, 1.0)[0];
  while( (x < a) || (x > b) ) x = rnorm(1, 0.0, 1.0)[0];
  return x;
}


// exp_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using exponential rejection sampling.
// ======================================================================

// [[Rcpp::export]]

double exp_rs(double a, double b)
{
  double  z, u, rate;

  //  Rprintf("in exp_rs");
  rate = 1/a;
  //1/a

  // Generate a proposal on (0, b-a)
  z = rexp(1, rate)[0];
  while(z > (b-a)) z = rexp(1, rate)[0];
  u = runif(1)[0];

  while( log(u) > (-0.5*z*z))
  {
    z = rexp(1, rate)[0];
    while(z > (b-a)) z = rexp(1, rate)[0];
    u = runif(1)[0];
  }
  return(z+a);
}




// rnorm_trunc( mu, sigma, lower, upper)
//
// generates one random normal RVs with mean 'mu' and standard
// deviation 'sigma', truncated to the interval (lower,upper), where
// lower can be -Inf and upper can be Inf.
//======================================================================

// [[Rcpp::export]]
double rnorm_trunc (double mu, double sigma, double lower, double upper)
{
  double a, b;
  double logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725;
  double z, tmp, lograt;

  a = (lower - mu)/sigma;
  b = (upper - mu)/sigma;

  // First scenario
  if( (a == R_NegInf) || (b == R_PosInf))
  {
    if(a == R_NegInf)
    {
      a = -b;
      b = R_PosInf;
    }

    // The two possibilities for this scenario
    if(a <= 0.45) z = norm_rs(a, b);
    else z = exp_rs(a, b);
    z = -z;
  }

  double output;
  output = sigma*z + mu;
  return (output);
}
