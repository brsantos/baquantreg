// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

// BayesQR
List BayesQR(double tau, arma::colvec y, arma::mat X, int itNum, int thin, arma::colvec betaValue, double sigmaValue, arma::vec vSampleInit, double priorVar, int refresh, bool quiet, bool tobit);
RcppExport SEXP baquantreg_BayesQR(SEXP tauSEXP, SEXP ySEXP, SEXP XSEXP, SEXP itNumSEXP, SEXP thinSEXP, SEXP betaValueSEXP, SEXP sigmaValueSEXP, SEXP vSampleInitSEXP, SEXP priorVarSEXP, SEXP refreshSEXP, SEXP quietSEXP, SEXP tobitSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type itNum(itNumSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type betaValue(betaValueSEXP);
    Rcpp::traits::input_parameter< double >::type sigmaValue(sigmaValueSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type vSampleInit(vSampleInitSEXP);
    Rcpp::traits::input_parameter< double >::type priorVar(priorVarSEXP);
    Rcpp::traits::input_parameter< int >::type refresh(refreshSEXP);
    Rcpp::traits::input_parameter< bool >::type quiet(quietSEXP);
    Rcpp::traits::input_parameter< bool >::type tobit(tobitSEXP);
    __result = Rcpp::wrap(BayesQR(tau, y, X, itNum, thin, betaValue, sigmaValue, vSampleInit, priorVar, refresh, quiet, tobit));
    return __result;
END_RCPP
}
// tpBayesQR
List tpBayesQR(double tau, arma::colvec y, arma::mat X, int itNum, int thin, arma::vec betaValue, double sigmaValue, arma::vec vSampleInit, arma::vec gammaValue, double sigmaGamma, int link, double priorVar, int refresh, bool quiet);
RcppExport SEXP baquantreg_tpBayesQR(SEXP tauSEXP, SEXP ySEXP, SEXP XSEXP, SEXP itNumSEXP, SEXP thinSEXP, SEXP betaValueSEXP, SEXP sigmaValueSEXP, SEXP vSampleInitSEXP, SEXP gammaValueSEXP, SEXP sigmaGammaSEXP, SEXP linkSEXP, SEXP priorVarSEXP, SEXP refreshSEXP, SEXP quietSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type itNum(itNumSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type betaValue(betaValueSEXP);
    Rcpp::traits::input_parameter< double >::type sigmaValue(sigmaValueSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type vSampleInit(vSampleInitSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gammaValue(gammaValueSEXP);
    Rcpp::traits::input_parameter< double >::type sigmaGamma(sigmaGammaSEXP);
    Rcpp::traits::input_parameter< int >::type link(linkSEXP);
    Rcpp::traits::input_parameter< double >::type priorVar(priorVarSEXP);
    Rcpp::traits::input_parameter< int >::type refresh(refreshSEXP);
    Rcpp::traits::input_parameter< bool >::type quiet(quietSEXP);
    __result = Rcpp::wrap(tpBayesQR(tau, y, X, itNum, thin, betaValue, sigmaValue, vSampleInit, gammaValue, sigmaGamma, link, priorVar, refresh, quiet));
    return __result;
END_RCPP
}
// ziTobitBayesQR
List ziTobitBayesQR(double tau, arma::colvec y, arma::mat X, int itNum, int thin, arma::colvec betaValue, double sigmaValue, arma::colvec gammaValue, double sigmaGamma, int link, double priorVar, int refresh, bool quiet, int burnin);
RcppExport SEXP baquantreg_ziTobitBayesQR(SEXP tauSEXP, SEXP ySEXP, SEXP XSEXP, SEXP itNumSEXP, SEXP thinSEXP, SEXP betaValueSEXP, SEXP sigmaValueSEXP, SEXP gammaValueSEXP, SEXP sigmaGammaSEXP, SEXP linkSEXP, SEXP priorVarSEXP, SEXP refreshSEXP, SEXP quietSEXP, SEXP burninSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type itNum(itNumSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type betaValue(betaValueSEXP);
    Rcpp::traits::input_parameter< double >::type sigmaValue(sigmaValueSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type gammaValue(gammaValueSEXP);
    Rcpp::traits::input_parameter< double >::type sigmaGamma(sigmaGammaSEXP);
    Rcpp::traits::input_parameter< int >::type link(linkSEXP);
    Rcpp::traits::input_parameter< double >::type priorVar(priorVarSEXP);
    Rcpp::traits::input_parameter< int >::type refresh(refreshSEXP);
    Rcpp::traits::input_parameter< bool >::type quiet(quietSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    __result = Rcpp::wrap(ziTobitBayesQR(tau, y, X, itNum, thin, betaValue, sigmaValue, gammaValue, sigmaGamma, link, priorVar, refresh, quiet, burnin));
    return __result;
END_RCPP
}
// logLikelihoodAlpha
double logLikelihoodAlpha(double alpha, arma::mat aux, arma::mat diagU, arma::mat covMat, double jitter);
RcppExport SEXP baquantreg_logLikelihoodAlpha(SEXP alphaSEXP, SEXP auxSEXP, SEXP diagUSEXP, SEXP covMatSEXP, SEXP jitterSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aux(auxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type diagU(diagUSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type covMat(covMatSEXP);
    Rcpp::traits::input_parameter< double >::type jitter(jitterSEXP);
    __result = Rcpp::wrap(logLikelihoodAlpha(alpha, aux, diagU, covMat, jitter));
    return __result;
END_RCPP
}
// logLikelihoodAlpha2
double logLikelihoodAlpha2(double alpha, arma::mat aux, arma::mat diagU, arma::mat covMat, arma::mat covMat2, arma::mat covMatAux, double jitter, arma::uvec indices, int m);
RcppExport SEXP baquantreg_logLikelihoodAlpha2(SEXP alphaSEXP, SEXP auxSEXP, SEXP diagUSEXP, SEXP covMatSEXP, SEXP covMat2SEXP, SEXP covMatAuxSEXP, SEXP jitterSEXP, SEXP indicesSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aux(auxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type diagU(diagUSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type covMat(covMatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type covMat2(covMat2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type covMatAux(covMatAuxSEXP);
    Rcpp::traits::input_parameter< double >::type jitter(jitterSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    __result = Rcpp::wrap(logLikelihoodAlpha2(alpha, aux, diagU, covMat, covMat2, covMatAux, jitter, indices, m));
    return __result;
END_RCPP
}
// mhAlpha
double mhAlpha(double alpha, arma::mat aux, arma::mat diagU, arma::mat covMat, double tuneA, double jitter);
RcppExport SEXP baquantreg_mhAlpha(SEXP alphaSEXP, SEXP auxSEXP, SEXP diagUSEXP, SEXP covMatSEXP, SEXP tuneASEXP, SEXP jitterSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aux(auxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type diagU(diagUSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type covMat(covMatSEXP);
    Rcpp::traits::input_parameter< double >::type tuneA(tuneASEXP);
    Rcpp::traits::input_parameter< double >::type jitter(jitterSEXP);
    __result = Rcpp::wrap(mhAlpha(alpha, aux, diagU, covMat, tuneA, jitter));
    return __result;
END_RCPP
}
// mhAlpha2
double mhAlpha2(double alpha, arma::mat aux, arma::mat diagU, arma::mat covMat, arma::mat covMat2, arma::mat covMatAux, double tuneA, double jitter, arma::uvec indices, int m);
RcppExport SEXP baquantreg_mhAlpha2(SEXP alphaSEXP, SEXP auxSEXP, SEXP diagUSEXP, SEXP covMatSEXP, SEXP covMat2SEXP, SEXP covMatAuxSEXP, SEXP tuneASEXP, SEXP jitterSEXP, SEXP indicesSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aux(auxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type diagU(diagUSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type covMat(covMatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type covMat2(covMat2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type covMatAux(covMatAuxSEXP);
    Rcpp::traits::input_parameter< double >::type tuneA(tuneASEXP);
    Rcpp::traits::input_parameter< double >::type jitter(jitterSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    __result = Rcpp::wrap(mhAlpha2(alpha, aux, diagU, covMat, covMat2, covMatAux, tuneA, jitter, indices, m));
    return __result;
END_RCPP
}
// rgigRcpp
double rgigRcpp(double chi, double psi, double lambda);
RcppExport SEXP baquantreg_rgigRcpp(SEXP chiSEXP, SEXP psiSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type chi(chiSEXP);
    Rcpp::traits::input_parameter< double >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    __result = Rcpp::wrap(rgigRcpp(chi, psi, lambda));
    return __result;
END_RCPP
}
// logLikelihoodKappa
double logLikelihoodKappa(double kappa, arma::mat aux, arma::mat diagU, arma::mat covMat, arma::mat covMatInv, arma::vec spCoord1, arma::vec spCoord2, double alpha, double jitter, bool newkappa);
RcppExport SEXP baquantreg_logLikelihoodKappa(SEXP kappaSEXP, SEXP auxSEXP, SEXP diagUSEXP, SEXP covMatSEXP, SEXP covMatInvSEXP, SEXP spCoord1SEXP, SEXP spCoord2SEXP, SEXP alphaSEXP, SEXP jitterSEXP, SEXP newkappaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aux(auxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type diagU(diagUSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type covMat(covMatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type covMatInv(covMatInvSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type spCoord1(spCoord1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type spCoord2(spCoord2SEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type jitter(jitterSEXP);
    Rcpp::traits::input_parameter< bool >::type newkappa(newkappaSEXP);
    __result = Rcpp::wrap(logLikelihoodKappa(kappa, aux, diagU, covMat, covMatInv, spCoord1, spCoord2, alpha, jitter, newkappa));
    return __result;
END_RCPP
}
// mhKappa
double mhKappa(double kappa, arma::vec spCoord1, arma::vec spCoord2, arma::mat aux, arma::mat diagU, arma::mat covMat, arma::mat covMatInv, double tuneParam, double alpha, double jitter);
RcppExport SEXP baquantreg_mhKappa(SEXP kappaSEXP, SEXP spCoord1SEXP, SEXP spCoord2SEXP, SEXP auxSEXP, SEXP diagUSEXP, SEXP covMatSEXP, SEXP covMatInvSEXP, SEXP tuneParamSEXP, SEXP alphaSEXP, SEXP jitterSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type spCoord1(spCoord1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type spCoord2(spCoord2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type aux(auxSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type diagU(diagUSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type covMat(covMatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type covMatInv(covMatInvSEXP);
    Rcpp::traits::input_parameter< double >::type tuneParam(tuneParamSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type jitter(jitterSEXP);
    __result = Rcpp::wrap(mhKappa(kappa, spCoord1, spCoord2, aux, diagU, covMat, covMatInv, tuneParam, alpha, jitter));
    return __result;
END_RCPP
}
// spBayesQR
List spBayesQR(double tau, arma::colvec y, arma::mat X, int itNum, int thin, arma::colvec betaValue, double sigmaValue, arma::vec spCoord1, arma::vec spCoord2, double lambda, double tuneP, double alphaValue, double tuneA, double priorVar, int refresh, bool quiet, double jitter, bool includeAlpha, double tuneV, int kMT);
RcppExport SEXP baquantreg_spBayesQR(SEXP tauSEXP, SEXP ySEXP, SEXP XSEXP, SEXP itNumSEXP, SEXP thinSEXP, SEXP betaValueSEXP, SEXP sigmaValueSEXP, SEXP spCoord1SEXP, SEXP spCoord2SEXP, SEXP lambdaSEXP, SEXP tunePSEXP, SEXP alphaValueSEXP, SEXP tuneASEXP, SEXP priorVarSEXP, SEXP refreshSEXP, SEXP quietSEXP, SEXP jitterSEXP, SEXP includeAlphaSEXP, SEXP tuneVSEXP, SEXP kMTSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type itNum(itNumSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type betaValue(betaValueSEXP);
    Rcpp::traits::input_parameter< double >::type sigmaValue(sigmaValueSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type spCoord1(spCoord1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type spCoord2(spCoord2SEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type tuneP(tunePSEXP);
    Rcpp::traits::input_parameter< double >::type alphaValue(alphaValueSEXP);
    Rcpp::traits::input_parameter< double >::type tuneA(tuneASEXP);
    Rcpp::traits::input_parameter< double >::type priorVar(priorVarSEXP);
    Rcpp::traits::input_parameter< int >::type refresh(refreshSEXP);
    Rcpp::traits::input_parameter< bool >::type quiet(quietSEXP);
    Rcpp::traits::input_parameter< double >::type jitter(jitterSEXP);
    Rcpp::traits::input_parameter< bool >::type includeAlpha(includeAlphaSEXP);
    Rcpp::traits::input_parameter< double >::type tuneV(tuneVSEXP);
    Rcpp::traits::input_parameter< int >::type kMT(kMTSEXP);
    __result = Rcpp::wrap(spBayesQR(tau, y, X, itNum, thin, betaValue, sigmaValue, spCoord1, spCoord2, lambda, tuneP, alphaValue, tuneA, priorVar, refresh, quiet, jitter, includeAlpha, tuneV, kMT));
    return __result;
END_RCPP
}
// sppBayesQR
List sppBayesQR(double tau, arma::colvec y, arma::mat X, int itNum, int thin, arma::colvec betaValue, double sigmaValue, arma::vec spCoord1, arma::vec spCoord2, double lambda, double tuneP, arma::uvec indices, int m, double alphaValue, double tuneA, double priorVar, bool quiet, int refresh, double jitter, bool includeAlpha, double tuneV, int kMT);
RcppExport SEXP baquantreg_sppBayesQR(SEXP tauSEXP, SEXP ySEXP, SEXP XSEXP, SEXP itNumSEXP, SEXP thinSEXP, SEXP betaValueSEXP, SEXP sigmaValueSEXP, SEXP spCoord1SEXP, SEXP spCoord2SEXP, SEXP lambdaSEXP, SEXP tunePSEXP, SEXP indicesSEXP, SEXP mSEXP, SEXP alphaValueSEXP, SEXP tuneASEXP, SEXP priorVarSEXP, SEXP quietSEXP, SEXP refreshSEXP, SEXP jitterSEXP, SEXP includeAlphaSEXP, SEXP tuneVSEXP, SEXP kMTSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type itNum(itNumSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type betaValue(betaValueSEXP);
    Rcpp::traits::input_parameter< double >::type sigmaValue(sigmaValueSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type spCoord1(spCoord1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type spCoord2(spCoord2SEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type tuneP(tunePSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type indices(indicesSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type alphaValue(alphaValueSEXP);
    Rcpp::traits::input_parameter< double >::type tuneA(tuneASEXP);
    Rcpp::traits::input_parameter< double >::type priorVar(priorVarSEXP);
    Rcpp::traits::input_parameter< bool >::type quiet(quietSEXP);
    Rcpp::traits::input_parameter< int >::type refresh(refreshSEXP);
    Rcpp::traits::input_parameter< double >::type jitter(jitterSEXP);
    Rcpp::traits::input_parameter< bool >::type includeAlpha(includeAlphaSEXP);
    Rcpp::traits::input_parameter< double >::type tuneV(tuneVSEXP);
    Rcpp::traits::input_parameter< int >::type kMT(kMTSEXP);
    __result = Rcpp::wrap(sppBayesQR(tau, y, X, itNum, thin, betaValue, sigmaValue, spCoord1, spCoord2, lambda, tuneP, indices, m, alphaValue, tuneA, priorVar, quiet, refresh, jitter, includeAlpha, tuneV, kMT));
    return __result;
END_RCPP
}
