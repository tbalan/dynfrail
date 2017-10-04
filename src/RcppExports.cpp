// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Estep_id
NumericVector Estep_id(NumericVector events, NumericVector cvec, double aalpha, double ggamma, int dist, double pvfm, NumericVector times, double llambda);
RcppExport SEXP _dynfrail_Estep_id(SEXP eventsSEXP, SEXP cvecSEXP, SEXP aalphaSEXP, SEXP ggammaSEXP, SEXP distSEXP, SEXP pvfmSEXP, SEXP timesSEXP, SEXP llambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cvec(cvecSEXP);
    Rcpp::traits::input_parameter< double >::type aalpha(aalphaSEXP);
    Rcpp::traits::input_parameter< double >::type ggamma(ggammaSEXP);
    Rcpp::traits::input_parameter< int >::type dist(distSEXP);
    Rcpp::traits::input_parameter< double >::type pvfm(pvfmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type times(timesSEXP);
    Rcpp::traits::input_parameter< double >::type llambda(llambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(Estep_id(events, cvec, aalpha, ggamma, dist, pvfm, times, llambda));
    return rcpp_result_gen;
END_RCPP
}
// Vcov_adj_id2
int Vcov_adj_id2(NumericVector events, NumericVector cvec, double aalpha, double ggamma, int dist, double pvfm, NumericVector times, double llambda, NumericVector elp, NumericMatrix xelph, NumericMatrix tau, Rcpp::List interval_rows, NumericVector ez);
RcppExport SEXP _dynfrail_Vcov_adj_id2(SEXP eventsSEXP, SEXP cvecSEXP, SEXP aalphaSEXP, SEXP ggammaSEXP, SEXP distSEXP, SEXP pvfmSEXP, SEXP timesSEXP, SEXP llambdaSEXP, SEXP elpSEXP, SEXP xelphSEXP, SEXP tauSEXP, SEXP interval_rowsSEXP, SEXP ezSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cvec(cvecSEXP);
    Rcpp::traits::input_parameter< double >::type aalpha(aalphaSEXP);
    Rcpp::traits::input_parameter< double >::type ggamma(ggammaSEXP);
    Rcpp::traits::input_parameter< int >::type dist(distSEXP);
    Rcpp::traits::input_parameter< double >::type pvfm(pvfmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type times(timesSEXP);
    Rcpp::traits::input_parameter< double >::type llambda(llambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type elp(elpSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xelph(xelphSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type interval_rows(interval_rowsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ez(ezSEXP);
    rcpp_result_gen = Rcpp::wrap(Vcov_adj_id2(events, cvec, aalpha, ggamma, dist, pvfm, times, llambda, elp, xelph, tau, interval_rows, ez));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dynfrail_Estep_id", (DL_FUNC) &_dynfrail_Estep_id, 8},
    {"_dynfrail_Vcov_adj_id2", (DL_FUNC) &_dynfrail_Vcov_adj_id2, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_dynfrail(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
