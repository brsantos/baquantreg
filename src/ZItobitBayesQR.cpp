#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

using namespace Rcpp;   // inline does that for us already
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

struct bessel_params {double beta; double lambda; double m;};

double postDensMHstep(colvec parameters, String link, mat X, colvec varInd,
                      colvec priorBeta, mat priorSigma){

  int n = X.n_rows;

  double priorWeight, likelihoodInf, postDens;
  arma::colvec valLik(n), aux(n);
  NumericVector aux1(1);

  priorWeight = -0.5 * as_scalar((parameters.t() - priorBeta.t()) * priorSigma.i() * (parameters - priorBeta));

  aux = X * parameters;

  if (link=="logit"){
    for (int j=0; j < n; j++){
      if(varInd(j)==1) valLik[j] = log(exp(aux[j])/(1 + exp(aux[j])));
      else valLik[j] = -log(1 + exp(aux[j]));
    }
  }
  else{
    for (int j=0; j < n; j++){
      aux1[0] = as_scalar(X.row(j) * parameters);
      if(varInd(j)==1) valLik[j] = log(pnorm(aux1, 0.0, 1.0)[0]);
      else valLik[j] = log(1-pnorm(aux1, 0.0, 1.0)[0]);
    }
  }

  likelihoodInf = sum(valLik);
  postDens = likelihoodInf + priorWeight;
  return postDens;
}

double predicProb(colvec parameters, String link, rowvec x){
  double prob;
  NumericVector aux(1);
  if (link=="logit") prob = exp(as_scalar(x * parameters))/(1 + exp(as_scalar(x * parameters)));
  else{
    aux[0] = as_scalar(x * parameters);
    prob = pnorm(aux, 0.0, 1.0)[0];
  }
  return prob;
}

double funG(double y, void *params){
  struct bessel_params *p = (struct bessel_params *) params;

  double beta = (p->beta);
  double lambda = (p->lambda);
  double m = (p->m);

  return 0.5 * beta * pow(y,3) - pow(y,2) * (0.5 * beta * m + lambda + 1) + y * ((lambda - 1) * m - 0.5 * beta) + 0.5 * beta * m;
}

// [[Rcpp::export]]
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

// [[Rcpp::export]]
double rinvgammaRcpp(double shape, double scale){
  double out = pow(rgamma(1, shape, 1/scale)[0],-1);
  return out;
}

// [[Rcpp::export]]
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




////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
colvec mvrnormRcpp(colvec mu, mat Sigma){
  int p = mu.size();
  colvec eV(p), X(p);
  mat eVec(p,p);

  eig_sym(eV, eVec, Sigma);

  X = rnorm(p);

  mat produto = eVec * diagmat(sqrt(eV));

  colvec output;

  output =  mu + produto * X;
  return output;
}

// [[Rcpp::export]]
double Flaplace (double predictor, double tau, double sigma){
  double output;
  if (predictor > 0) output = tau * exp(-((1-tau)/sigma)*predictor);
  else output = 1 - (1-tau)*exp((tau/sigma)*predictor);
  return output;
}

// [[Rcpp::export]]
List ziTobitBayesQR(double tau, colvec y, mat X, int itNum, int thin,
                      colvec betaValue, double sigmaValue, colvec betaZeroValue,
                      double sigmaBetaZero, String link){

   RNGScope scope;

   int n, n1, n2, l, m, p, cont_n2, cont_l, cont_m;
   cont_n2 = 0;
   cont_l = 0;
   cont_m = 0;

   n = X.n_rows;
   p = X.n_cols;

   colvec varInd(n, fill::zeros);
   for (int dd = 0; dd < n; dd++){
     if (y[dd]==0){
       varInd[dd] = 1.0;
     }
   }

   n1 = sum(varInd);
   n2 = n - n1;

   double fObsLap, probCens;
   colvec censInd(n, fill::zeros);

   double theta, psi2, s0, n0, gama2, lambda, meanModel, sdModel, nTilde, sTilde, zValue,
            postNew, postAnt, probAcceptance, acceptRate, accept, termsSum;

   acceptRate = 0;

   mat betaSample(itNum/thin, p, fill::zeros), betaZeroSample(itNum/thin, p, fill::zeros);
   colvec sigmaSample(itNum/thin, fill::zeros);
   NumericVector InfPar(1), LimSup(1);

   colvec b0(p), mu(p), newBetaZeroValue(p), vetorUm(p, fill::ones), delta2(n);

   mat B0(p,p), sigmaMinusOne(p,p), diagU, sigmaMat(p,p), matIdsigmaBetaZero(p,p);

   // Hiperparâmetros da priori normal
   B0 = 100 * B0.eye(p,p);
   b0.fill(0);

   /* Constantes que dependem somente do tau estimado */
   theta = (1-2*tau)/(tau*(1-tau));
   psi2 = 2/(tau*(1-tau));

   /* Inicializando valores. */
   n0 = 3;
   s0 = 0.1;

   lambda = 0.5;
   InfPar[0] = R_NegInf; LimSup[0] = 0;
   /* *********************** */

   /* Inteiro responsável pela contagem das aceitações do modelo.*/
   accept = 0.0;

   matIdsigmaBetaZero = sigmaBetaZero * diagmat(vetorUm);

   IntegerVector seqRefresh = seq(1, itNum/100)*100;

   colvec zSample(n, fill::ones), aux(n), aux2;

   sigmaSample[0] = sigmaValue;

   mat matrizIndCens(itNum/thin, n, fill::zeros), matrizIndNCens(itNum/thin, n, fill::zeros);

   uvec ind_n2(n2, fill::zeros), ind_l(n1, fill::zeros), ind_m(n1, fill::zeros);

   double probZero;

   // Laço que vai fazer a atualização dos estados das cadeias.
   for(int k = 1; k < itNum; k++){
      for(int j = 0; j < thin; j++){

        if(is_true(any(k+1 == seqRefresh))){
          Rcout <<  "k = " << k + 1 << std::endl;
        }

        // Apenas para verificar
        // Rcout <<  "k = " << k << std::endl;

        //  Atualizando a variável Presença de censura

        ind_m.fill(0); ind_l.fill(0);
        censInd.fill(0);
        cont_m = 0; cont_l = 0;

        for (int ee = 0; ee < n; ee++){
          if (y[ee]==0){
            fObsLap = Flaplace(as_scalar(X.row(ee) * betaValue), tau, sigmaValue);
            probZero = 1-predicProb(betaZeroValue, link, X.row(ee));
            probCens = (probZero * fObsLap)/(1-probZero + probZero*fObsLap);
            if (runif(1)[0] < probCens){
              censInd[ee] = 1.0;
              matrizIndCens(k, ee) = 1.0;
              ind_m[cont_m] = ee;
              cont_m++;
            }
            else{
              ind_l[cont_l] = ee;
              matrizIndNCens(k, ee) = 1.0;
              cont_l++;
            }
          }
          else if (k==1){
            ind_n2[cont_n2] = ee;
            cont_n2++;
          }
        }

        m = sum(censInd);
        l = n1 - m;

        // if (l == 0) stop("All zero observations are estimated to be censored.");

        uvec ind_v(n2 + m, fill::zeros);

        if (m > 0) {
          uvec ind_m2 = ind_m.subvec(0, cont_m - 1);
          ind_v.subvec(0, m - 1) = ind_m2;
          ind_v.subvec(m, n2+m-1) = ind_n2;
        }
        else {
          ind_v = ind_n2;
        }

        colvec yS = y;
        mat X2 = X.rows(ind_v);
        colvec zSample2(n2+m, fill::ones);

        aux = X * betaValue;
        for (int aa = 0; aa < n; aa++){
           if (censInd(aa)==1.0){
              meanModel = aux[aa] + theta*zSample(aa);
              sdModel = sqrt(psi2*sigmaValue*zSample(aa));
              // yS[aa] = rtnormRcpp(InfPar, LimSup, meanModel, sdModel);
              yS[aa] = rnorm_trunc(meanModel, sdModel, InfPar[0], LimSup[0]);
//               Rcout << "aa = "<< aa << "; meanModel = " << meanModel;
//               Rcout << "; sdModel = " << sdModel;
//               Rcout << "; zSample = " << zSample(aa);
//               Rcout << "; yS = " << yS[aa] << std::endl;
           }
        }

        colvec y2 = yS(ind_v);

        delta2 = diagvec((1/(psi2*sigmaValue)) * diagmat(yS - aux) * diagmat(yS - aux));

        gama2 = 2/sigmaValue + (theta*theta)/(psi2*sigmaValue);

//         if (k>71){
//           Rcout << "sigmaValue = " << sigmaValue << std::endl;
//           Rcout << "gama2 = " << gama2  << std::endl;
//         }

        for(int o = 0; o < n; o++){
          if (varInd(o)==censInd(o)){
            delta2[o] = std::max(delta2[o], 1e-10);
            gama2 = std::max(gama2, 1e-10);
//             if (k>71){
//               Rcout << "o = " << o;
//               Rcout << "; zSample[o] = " << zSample[o];
//               Rcout << "; yS = " << yS[o];
//               Rcout << "; delta2 = " << delta2[o]  << std::endl;
//             }
            zSample[o] = rgigRcpp(delta2[o], gama2, lambda);
          }
        }

        // Atualizando a cadeia de sigma.

        zSample2 = zSample(ind_v);
        aux2 = aux.rows(ind_v);
        termsSum = as_scalar((y2 - aux2 - theta*zSample2).t() * diagmat(zSample2).i() * (y2 - aux2 - theta*zSample2));

        nTilde = n0 + 3*(n2+m);
        sTilde =  s0 + 2*sum(zSample2) + termsSum/psi2;

        sigmaValue = rinvgammaRcpp(nTilde/2,sTilde/2);

        // Atualizando beta do modelo quantílico
        diagU = diagmat(zSample2)*(psi2*sigmaValue);

        sigmaMinusOne = (X2.t() * diagU.i() * X2) + B0.i();

        sigmaMat = sigmaMinusOne.i();

        mu = sigmaMat * (B0.i()*b0 +
          (1/(psi2 * sigmaValue))*(X2.t() * diagmat(zSample2).i() * (y2 - theta*zSample2)));

        betaValue = mvrnormRcpp(mu, sigmaMat);

        // ********************************************* //
        // Estimação da parte da massa pontual do modelo //

        newBetaZeroValue = mvrnormRcpp(betaZeroValue, matIdsigmaBetaZero);

        postAnt = postDensMHstep(betaZeroValue, link, X, varInd - censInd, b0, B0);
        postNew = postDensMHstep(newBetaZeroValue, link, X, varInd - censInd, b0, B0);

        probAcceptance = std::min(0.0,postNew-postAnt);

        if(log(runif(1)[0])<probAcceptance){
          betaZeroValue = newBetaZeroValue;
          accept++;
        }

        // ********************************************* //

      }


      // Armazenando valores finais considerando parâmetro thin
      for(int jj = 0; jj < p; jj++){
        betaSample(k,jj) = betaValue[jj];
        betaZeroSample(k,jj) = betaZeroValue[jj];
      }

      sigmaSample[k] = sigmaValue;
   }

   acceptRate = accept/itNum;

   return List::create(
        Named("BetaSample") = betaSample,
        Named("BetaZeroSample") = betaZeroSample,
        Named("acceptRate") = acceptRate,
        Named("tau") = tau,
        Named("SigmaSample") = sigmaSample,
        Named("matrizIndCens") = matrizIndCens);
}
