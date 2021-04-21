// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// LogSumExp
double LogSumExp(NumericVector logV);
RcppExport SEXP _classNB_LogSumExp(SEXP logVSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type logV(logVSEXP);
    rcpp_result_gen = Rcpp::wrap(LogSumExp(logV));
    return rcpp_result_gen;
END_RCPP
}
// C_sum2Mat
NumericMatrix C_sum2Mat(NumericMatrix A, NumericMatrix B);
RcppExport SEXP _classNB_C_sum2Mat(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type A(ASEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(C_sum2Mat(A, B));
    return rcpp_result_gen;
END_RCPP
}
// ParamVec
SEXP ParamVec(SEXP A_, SEXP B_, SEXP C_, SEXP D_);
RcppExport SEXP _classNB_ParamVec(SEXP A_SEXP, SEXP B_SEXP, SEXP C_SEXP, SEXP D_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type A_(A_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type B_(B_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type C_(C_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type D_(D_SEXP);
    rcpp_result_gen = Rcpp::wrap(ParamVec(A_, B_, C_, D_));
    return rcpp_result_gen;
END_RCPP
}
// OrdParam
List OrdParam(NumericVector initPr, NumericMatrix transPr, NumericVector esAvg, NumericVector eSize, NumericVector eProb);
RcppExport SEXP _classNB_OrdParam(SEXP initPrSEXP, SEXP transPrSEXP, SEXP esAvgSEXP, SEXP eSizeSEXP, SEXP eProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type initPr(initPrSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type transPr(transPrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type esAvg(esAvgSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eSize(eSizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eProb(eProbSEXP);
    rcpp_result_gen = Rcpp::wrap(OrdParam(initPr, transPr, esAvg, eSize, eProb));
    return rcpp_result_gen;
END_RCPP
}
// ForwardBackward_HMM
List ForwardBackward_HMM(NumericVector initPr, NumericMatrix transPr, NumericMatrix emisPr);
RcppExport SEXP _classNB_ForwardBackward_HMM(SEXP initPrSEXP, SEXP transPrSEXP, SEXP emisPrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type initPr(initPrSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type transPr(transPrSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type emisPr(emisPrSEXP);
    rcpp_result_gen = Rcpp::wrap(ForwardBackward_HMM(initPr, transPr, emisPr));
    return rcpp_result_gen;
END_RCPP
}
// EmisPr_NB
NumericMatrix EmisPr_NB(NumericVector X, NumericVector eSize, NumericVector eProb);
RcppExport SEXP _classNB_EmisPr_NB(SEXP XSEXP, SEXP eSizeSEXP, SEXP eProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eSize(eSizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eProb(eProbSEXP);
    rcpp_result_gen = Rcpp::wrap(EmisPr_NB(X, eSize, eProb));
    return rcpp_result_gen;
END_RCPP
}
// EmisPr_mixNB
NumericMatrix EmisPr_mixNB(NumericVector X, IntegerVector component, NumericVector eWeight, NumericVector eSize, NumericVector eProb);
RcppExport SEXP _classNB_EmisPr_mixNB(SEXP XSEXP, SEXP componentSEXP, SEXP eWeightSEXP, SEXP eSizeSEXP, SEXP eProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type component(componentSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eWeight(eWeightSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eSize(eSizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eProb(eProbSEXP);
    rcpp_result_gen = Rcpp::wrap(EmisPr_mixNB(X, component, eWeight, eSize, eProb));
    return rcpp_result_gen;
END_RCPP
}
// EM_NB
List EM_NB(NumericVector X, NumericVector initPr, NumericMatrix transPr, NumericVector eSize, NumericVector eProb, int iterMax, double EPSILON);
RcppExport SEXP _classNB_EM_NB(SEXP XSEXP, SEXP initPrSEXP, SEXP transPrSEXP, SEXP eSizeSEXP, SEXP eProbSEXP, SEXP iterMaxSEXP, SEXP EPSILONSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type initPr(initPrSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type transPr(transPrSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eSize(eSizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eProb(eProbSEXP);
    Rcpp::traits::input_parameter< int >::type iterMax(iterMaxSEXP);
    Rcpp::traits::input_parameter< double >::type EPSILON(EPSILONSEXP);
    rcpp_result_gen = Rcpp::wrap(EM_NB(X, initPr, transPr, eSize, eProb, iterMax, EPSILON));
    return rcpp_result_gen;
END_RCPP
}
// Viterbi_HMM
SEXP Viterbi_HMM(SEXP initPr_, SEXP transPr_, SEXP emisPr_);
RcppExport SEXP _classNB_Viterbi_HMM(SEXP initPr_SEXP, SEXP transPr_SEXP, SEXP emisPr_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type initPr_(initPr_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type transPr_(transPr_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type emisPr_(emisPr_SEXP);
    rcpp_result_gen = Rcpp::wrap(Viterbi_HMM(initPr_, transPr_, emisPr_));
    return rcpp_result_gen;
END_RCPP
}
// InitEM
List InitEM(NumericVector X, int nbS);
RcppExport SEXP _classNB_InitEM(SEXP XSEXP, SEXP nbSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nbS(nbSSEXP);
    rcpp_result_gen = Rcpp::wrap(InitEM(X, nbS));
    return rcpp_result_gen;
END_RCPP
}
// HMM_NB_CF
List HMM_NB_CF(NumericVector X, int nbS, bool viterbi, Nullable<NumericVector> initPr, Nullable<NumericMatrix> transPr, Nullable<NumericVector> eSize, Nullable<NumericVector> eProb, int iterMax, double EPSILON);
RcppExport SEXP _classNB_HMM_NB_CF(SEXP XSEXP, SEXP nbSSEXP, SEXP viterbiSEXP, SEXP initPrSEXP, SEXP transPrSEXP, SEXP eSizeSEXP, SEXP eProbSEXP, SEXP iterMaxSEXP, SEXP EPSILONSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nbS(nbSSEXP);
    Rcpp::traits::input_parameter< bool >::type viterbi(viterbiSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type initPr(initPrSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericMatrix> >::type transPr(transPrSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type eSize(eSizeSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type eProb(eProbSEXP);
    Rcpp::traits::input_parameter< int >::type iterMax(iterMaxSEXP);
    Rcpp::traits::input_parameter< double >::type EPSILON(EPSILONSEXP);
    rcpp_result_gen = Rcpp::wrap(HMM_NB_CF(X, nbS, viterbi, initPr, transPr, eSize, eProb, iterMax, EPSILON));
    return rcpp_result_gen;
END_RCPP
}
// HMM_NB_NS
List HMM_NB_NS(NumericVector X, int nbS, bool viterbi, int iterMax, double EPSILON);
RcppExport SEXP _classNB_HMM_NB_NS(SEXP XSEXP, SEXP nbSSEXP, SEXP viterbiSEXP, SEXP iterMaxSEXP, SEXP EPSILONSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nbS(nbSSEXP);
    Rcpp::traits::input_parameter< bool >::type viterbi(viterbiSEXP);
    Rcpp::traits::input_parameter< int >::type iterMax(iterMaxSEXP);
    Rcpp::traits::input_parameter< double >::type EPSILON(EPSILONSEXP);
    rcpp_result_gen = Rcpp::wrap(HMM_NB_NS(X, nbS, viterbi, iterMax, EPSILON));
    return rcpp_result_gen;
END_RCPP
}
// mixNB_NB_CF
List mixNB_NB_CF(NumericVector X, int nbS, int iterMax, double EPSILON);
RcppExport SEXP _classNB_mixNB_NB_CF(SEXP XSEXP, SEXP nbSSEXP, SEXP iterMaxSEXP, SEXP EPSILONSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nbS(nbSSEXP);
    Rcpp::traits::input_parameter< int >::type iterMax(iterMaxSEXP);
    Rcpp::traits::input_parameter< double >::type EPSILON(EPSILONSEXP);
    rcpp_result_gen = Rcpp::wrap(mixNB_NB_CF(X, nbS, iterMax, EPSILON));
    return rcpp_result_gen;
END_RCPP
}
// mixNB_NB_NS
List mixNB_NB_NS(NumericVector X, int nbS, int iterMax, double EPSILON);
RcppExport SEXP _classNB_mixNB_NB_NS(SEXP XSEXP, SEXP nbSSEXP, SEXP iterMaxSEXP, SEXP EPSILONSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type nbS(nbSSEXP);
    Rcpp::traits::input_parameter< int >::type iterMax(iterMaxSEXP);
    Rcpp::traits::input_parameter< double >::type EPSILON(EPSILONSEXP);
    rcpp_result_gen = Rcpp::wrap(mixNB_NB_NS(X, nbS, iterMax, EPSILON));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_classNB_LogSumExp", (DL_FUNC) &_classNB_LogSumExp, 1},
    {"_classNB_C_sum2Mat", (DL_FUNC) &_classNB_C_sum2Mat, 2},
    {"_classNB_ParamVec", (DL_FUNC) &_classNB_ParamVec, 4},
    {"_classNB_OrdParam", (DL_FUNC) &_classNB_OrdParam, 5},
    {"_classNB_ForwardBackward_HMM", (DL_FUNC) &_classNB_ForwardBackward_HMM, 3},
    {"_classNB_EmisPr_NB", (DL_FUNC) &_classNB_EmisPr_NB, 3},
    {"_classNB_EmisPr_mixNB", (DL_FUNC) &_classNB_EmisPr_mixNB, 5},
    {"_classNB_EM_NB", (DL_FUNC) &_classNB_EM_NB, 7},
    {"_classNB_Viterbi_HMM", (DL_FUNC) &_classNB_Viterbi_HMM, 3},
    {"_classNB_InitEM", (DL_FUNC) &_classNB_InitEM, 2},
    {"_classNB_HMM_NB_CF", (DL_FUNC) &_classNB_HMM_NB_CF, 9},
    {"_classNB_HMM_NB_NS", (DL_FUNC) &_classNB_HMM_NB_NS, 5},
    {"_classNB_mixNB_NB_CF", (DL_FUNC) &_classNB_mixNB_NB_CF, 4},
    {"_classNB_mixNB_NB_NS", (DL_FUNC) &_classNB_mixNB_NB_NS, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_classNB(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}