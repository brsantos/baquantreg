#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppGSL)]]

struct bessel_params {double beta; double lambda; double m;};

double funG(double y, void *params){
  struct bessel_params *p = (struct bessel_params *) params;

  double beta = (p->beta);
  double lambda = (p->lambda);
  double m = (p->m);

  return 0.5 * beta * pow(y,3) - pow(y,2) * (0.5 * beta * m + lambda + 1) + y * ((lambda - 1) * m - 0.5 * beta) + 0.5 * beta * m;
}

double rgigRcpp(double chi, double psi, double lambda){
  double alpha = sqrt(psi/chi);
  double beta = sqrt(psi * chi);
  double m = std::max((lambda - 1 + sqrt(pow(lambda - 1,2) + pow(beta,2)))/beta, 1e-8);

  const gsl_root_fsolver_type *T1;
  const gsl_root_fsolver_type *T2;
  gsl_root_fsolver *s1;
  gsl_root_fsolver *s2;

  gsl_function F1, F2;
  struct bessel_params params = {beta, lambda, m};

  double upper = m;
  while (funG(upper, &params) <= 0) upper = 2 * upper;

  // Encontrar soluções de funG (yP e yM) //
  F1.function = &funG; F2.function = &funG;
  F1.params = &params; F2.params = &params;

  T1 = gsl_root_fsolver_brent;
  T2 = gsl_root_fsolver_brent;
  s1 = gsl_root_fsolver_alloc (T1);
  s2 = gsl_root_fsolver_alloc (T2);

  int status1, status2;

  double x_lo1 = 0.0, x_hi1 = m;
  double x_lo2 = m, x_hi2 = upper;

  gsl_root_fsolver_set (s1, &F1, x_lo1, x_hi1);
  gsl_root_fsolver_set (s2, &F2, x_lo2, x_hi2);

  double yM, yP;

  // Obtendo primeira raiz.
  do
  {
    status1 = gsl_root_fsolver_iterate (s1);
    yM = gsl_root_fsolver_root (s1);
    x_lo1 = gsl_root_fsolver_x_lower (s1);
    x_hi1 = gsl_root_fsolver_x_upper (s1);
    status1 = gsl_root_test_interval (x_lo1, x_hi1,0, 0.000001);
  }
  while (status1 == GSL_CONTINUE);
  gsl_root_fsolver_free (s1);

  // Obtendo segunda raiz.
  do
  {
    status2 = gsl_root_fsolver_iterate (s2);
    yP = gsl_root_fsolver_root (s2);
    x_lo2 = gsl_root_fsolver_x_lower (s2);
    x_hi2 = gsl_root_fsolver_x_upper (s2);
    status2 = gsl_root_test_interval (x_lo2, x_hi2,0, 0.000001);
  }
  while (status2 == GSL_CONTINUE);
  gsl_root_fsolver_free (s2);

  double a = (yP - m) * pow(yP/m,0.5 * (lambda - 1)) * exp(-0.25 * beta * (yP + 1/yP - m - 1/m));
  double b = (yM - m) * pow(yM/m,0.5 * (lambda - 1)) * exp(-0.25 * beta * (yM + 1/yM - m - 1/m));
  double c = -0.25 * beta * (m + 1/m) + 0.5 * (lambda - 1) * log(m);

  double output;

  int needValue = 1;
  double Y;
  while (needValue==1) {
    double R1 = runif(1)[0];
    double R2 = runif(1)[0];
    Y = m + a * R2/R1 + b * (1 - R2)/R1;
    if (Y > 0) {
      if (-log(R1) >= -0.5 * (lambda - 1) * log(Y) + 0.25 * beta * (Y + 1/Y) + c){
        needValue = 0;
      }
    }
  }
  output = Y/alpha;
  return output;
}