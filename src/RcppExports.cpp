// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sce_lik
double sce_lik(arma::cx_cube R, arma::cx_cube X, arma::uvec des, arma::uvec des_pos, arma::uvec des_n, arma::uvec prune, arma::vec T, double S);
RcppExport SEXP _sce_sce_lik(SEXP RSEXP, SEXP XSEXP, SEXP desSEXP, SEXP des_posSEXP, SEXP des_nSEXP, SEXP pruneSEXP, SEXP TSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cx_cube >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::cx_cube >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type des(desSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type des_pos(des_posSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type des_n(des_nSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type prune(pruneSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type T(TSEXP);
    Rcpp::traits::input_parameter< double >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(sce_lik(R, X, des, des_pos, des_n, prune, T, S));
    return rcpp_result_gen;
END_RCPP
}
// sce_rec
arma::cube sce_rec(arma::cx_cube R, arma::cx_cube bw_R, arma::cx_cube X, arma::uvec des, arma::uvec des_pos, arma::uvec des_n, arma::uvec prune, arma::vec T, arma::uvec tip);
RcppExport SEXP _sce_sce_rec(SEXP RSEXP, SEXP bw_RSEXP, SEXP XSEXP, SEXP desSEXP, SEXP des_posSEXP, SEXP des_nSEXP, SEXP pruneSEXP, SEXP TSEXP, SEXP tipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cx_cube >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::cx_cube >::type bw_R(bw_RSEXP);
    Rcpp::traits::input_parameter< arma::cx_cube >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type des(desSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type des_pos(des_posSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type des_n(des_nSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type prune(pruneSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type T(TSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type tip(tipSEXP);
    rcpp_result_gen = Rcpp::wrap(sce_rec(R, bw_R, X, des, des_pos, des_n, prune, T, tip));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sce_sce_lik", (DL_FUNC) &_sce_sce_lik, 8},
    {"_sce_sce_rec", (DL_FUNC) &_sce_sce_rec, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_sce(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
