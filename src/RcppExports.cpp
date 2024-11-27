// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// S1_wls
arma::vec S1_wls(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S1_wls(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S1_wls(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S2_wls
arma::vec S2_wls(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S2_wls(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S2_wls(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S3_wls
arma::vec S3_wls(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S3_wls(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S3_wls(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S1_yang
arma::vec S1_yang(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S1_yang(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S1_yang(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S2_yang
arma::vec S2_yang(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S2_yang(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S2_yang(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S3_yang
arma::vec S3_yang(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S3_yang(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S3_yang(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// tertiles_yang
arma::vec tertiles_yang(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_tertiles_yang(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(tertiles_yang(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// quintiles_yang
arma::vec quintiles_yang(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_quintiles_yang(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(quintiles_yang(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// deciles_yang
arma::vec deciles_yang(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_deciles_yang(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(deciles_yang(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S1_yang_logit
arma::vec S1_yang_logit(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S1_yang_logit(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S1_yang_logit(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S2_yang_logit
arma::vec S2_yang_logit(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S2_yang_logit(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S2_yang_logit(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S3_yang_logit
arma::vec S3_yang_logit(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S3_yang_logit(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S3_yang_logit(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// tertiles_yang_logit
arma::vec tertiles_yang_logit(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_tertiles_yang_logit(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(tertiles_yang_logit(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// quintiles_yang_logit
arma::vec quintiles_yang_logit(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_quintiles_yang_logit(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(quintiles_yang_logit(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// deciles_yang_logit
arma::vec deciles_yang_logit(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_deciles_yang_logit(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(deciles_yang_logit(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S1_yang_lap
arma::vec S1_yang_lap(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S1_yang_lap(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S1_yang_lap(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S2_yang_lap
arma::vec S2_yang_lap(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S2_yang_lap(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S2_yang_lap(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S3_yang_lap
arma::vec S3_yang_lap(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S3_yang_lap(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S3_yang_lap(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// tertiles_yang_lap
arma::vec tertiles_yang_lap(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_tertiles_yang_lap(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(tertiles_yang_lap(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// quintiles_yang_lap
arma::vec quintiles_yang_lap(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_quintiles_yang_lap(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(quintiles_yang_lap(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// deciles_yang_lap
arma::vec deciles_yang_lap(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_deciles_yang_lap(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(deciles_yang_lap(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S1_bala
arma::vec S1_bala(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S1_bala(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S1_bala(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S2_bala
arma::vec S2_bala(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S2_bala(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S2_bala(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S3_bala
arma::vec S3_bala(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S3_bala(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S3_bala(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// tertiles_bala
arma::vec tertiles_bala(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_tertiles_bala(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(tertiles_bala(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// quintiles_bala
arma::vec quintiles_bala(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_quintiles_bala(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(quintiles_bala(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// deciles_bala
arma::vec deciles_bala(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_deciles_bala(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(deciles_bala(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S1_lap
arma::vec S1_lap(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S1_lap(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S1_lap(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S2_lap
arma::vec S2_lap(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S2_lap(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S2_lap(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S3_lap
arma::vec S3_lap(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S3_lap(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S3_lap(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S1_logit
arma::vec S1_logit(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S1_logit(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S1_logit(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S2_logit
arma::vec S2_logit(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S2_logit(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S2_logit(summary, n));
    return rcpp_result_gen;
END_RCPP
}
// S3_logit
arma::vec S3_logit(arma::vec summary, double n);
RcppExport SEXP _metaBLUEX_S3_logit(SEXP summarySEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type summary(summarySEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(S3_logit(summary, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_metaBLUEX_S1_wls", (DL_FUNC) &_metaBLUEX_S1_wls, 2},
    {"_metaBLUEX_S2_wls", (DL_FUNC) &_metaBLUEX_S2_wls, 2},
    {"_metaBLUEX_S3_wls", (DL_FUNC) &_metaBLUEX_S3_wls, 2},
    {"_metaBLUEX_S1_yang", (DL_FUNC) &_metaBLUEX_S1_yang, 2},
    {"_metaBLUEX_S2_yang", (DL_FUNC) &_metaBLUEX_S2_yang, 2},
    {"_metaBLUEX_S3_yang", (DL_FUNC) &_metaBLUEX_S3_yang, 2},
    {"_metaBLUEX_tertiles_yang", (DL_FUNC) &_metaBLUEX_tertiles_yang, 2},
    {"_metaBLUEX_quintiles_yang", (DL_FUNC) &_metaBLUEX_quintiles_yang, 2},
    {"_metaBLUEX_deciles_yang", (DL_FUNC) &_metaBLUEX_deciles_yang, 2},
    {"_metaBLUEX_S1_yang_logit", (DL_FUNC) &_metaBLUEX_S1_yang_logit, 2},
    {"_metaBLUEX_S2_yang_logit", (DL_FUNC) &_metaBLUEX_S2_yang_logit, 2},
    {"_metaBLUEX_S3_yang_logit", (DL_FUNC) &_metaBLUEX_S3_yang_logit, 2},
    {"_metaBLUEX_tertiles_yang_logit", (DL_FUNC) &_metaBLUEX_tertiles_yang_logit, 2},
    {"_metaBLUEX_quintiles_yang_logit", (DL_FUNC) &_metaBLUEX_quintiles_yang_logit, 2},
    {"_metaBLUEX_deciles_yang_logit", (DL_FUNC) &_metaBLUEX_deciles_yang_logit, 2},
    {"_metaBLUEX_S1_yang_lap", (DL_FUNC) &_metaBLUEX_S1_yang_lap, 2},
    {"_metaBLUEX_S2_yang_lap", (DL_FUNC) &_metaBLUEX_S2_yang_lap, 2},
    {"_metaBLUEX_S3_yang_lap", (DL_FUNC) &_metaBLUEX_S3_yang_lap, 2},
    {"_metaBLUEX_tertiles_yang_lap", (DL_FUNC) &_metaBLUEX_tertiles_yang_lap, 2},
    {"_metaBLUEX_quintiles_yang_lap", (DL_FUNC) &_metaBLUEX_quintiles_yang_lap, 2},
    {"_metaBLUEX_deciles_yang_lap", (DL_FUNC) &_metaBLUEX_deciles_yang_lap, 2},
    {"_metaBLUEX_S1_bala", (DL_FUNC) &_metaBLUEX_S1_bala, 2},
    {"_metaBLUEX_S2_bala", (DL_FUNC) &_metaBLUEX_S2_bala, 2},
    {"_metaBLUEX_S3_bala", (DL_FUNC) &_metaBLUEX_S3_bala, 2},
    {"_metaBLUEX_tertiles_bala", (DL_FUNC) &_metaBLUEX_tertiles_bala, 2},
    {"_metaBLUEX_quintiles_bala", (DL_FUNC) &_metaBLUEX_quintiles_bala, 2},
    {"_metaBLUEX_deciles_bala", (DL_FUNC) &_metaBLUEX_deciles_bala, 2},
    {"_metaBLUEX_S1_lap", (DL_FUNC) &_metaBLUEX_S1_lap, 2},
    {"_metaBLUEX_S2_lap", (DL_FUNC) &_metaBLUEX_S2_lap, 2},
    {"_metaBLUEX_S3_lap", (DL_FUNC) &_metaBLUEX_S3_lap, 2},
    {"_metaBLUEX_S1_logit", (DL_FUNC) &_metaBLUEX_S1_logit, 2},
    {"_metaBLUEX_S2_logit", (DL_FUNC) &_metaBLUEX_S2_logit, 2},
    {"_metaBLUEX_S3_logit", (DL_FUNC) &_metaBLUEX_S3_logit, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_metaBLUEX(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
