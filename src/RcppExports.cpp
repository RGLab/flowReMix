// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// MH
void MH(NumericMatrix lastMean, NumericMatrix estimatedRandomEffects, const NumericVector y, const NumericVector N, const NumericMatrix randomEffectSamp, const int i, IntegerVector popInd, arma::mat invcov, NumericVector accept, double iter, double rate, NumericVector unifs, const NumericVector nullEta, const NumericVector altEta, double updateLag);
RcppExport SEXP flowReMix_MH(SEXP lastMeanSEXP, SEXP estimatedRandomEffectsSEXP, SEXP ySEXP, SEXP NSEXP, SEXP randomEffectSampSEXP, SEXP iSEXP, SEXP popIndSEXP, SEXP invcovSEXP, SEXP acceptSEXP, SEXP iterSEXP, SEXP rateSEXP, SEXP unifsSEXP, SEXP nullEtaSEXP, SEXP altEtaSEXP, SEXP updateLagSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type lastMean(lastMeanSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type estimatedRandomEffects(estimatedRandomEffectsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type N(NSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix >::type randomEffectSamp(randomEffectSampSEXP);
    Rcpp::traits::input_parameter< const int >::type i(iSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type popInd(popIndSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type invcov(invcovSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type accept(acceptSEXP);
    Rcpp::traits::input_parameter< double >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type unifs(unifsSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type nullEta(nullEtaSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type altEta(altEtaSEXP);
    Rcpp::traits::input_parameter< double >::type updateLag(updateLagSEXP);
    MH(lastMean, estimatedRandomEffects, y, N, randomEffectSamp, i, popInd, invcov, accept, iter, rate, unifs, nullEta, altEta, updateLag);
    return R_NilValue;
END_RCPP
}
// zero
NumericVector zero(NumericVector accept);
RcppExport SEXP flowReMix_zero(SEXP acceptSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type accept(acceptSEXP);
    rcpp_result_gen = Rcpp::wrap(zero(accept));
    return rcpp_result_gen;
END_RCPP
}
// vecBetaBinomDens
NumericVector vecBetaBinomDens(NumericVector count, NumericVector N, NumericVector prob, double M);
RcppExport SEXP flowReMix_vecBetaBinomDens(SEXP countSEXP, SEXP NSEXP, SEXP probSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type count(countSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    Rcpp::traits::input_parameter< double >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(vecBetaBinomDens(count, N, prob, M));
    return rcpp_result_gen;
END_RCPP
}
// setNumericVectorToZero
void setNumericVectorToZero(NumericVector x);
RcppExport SEXP flowReMix_setNumericVectorToZero(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    setNumericVectorToZero(x);
    return R_NilValue;
END_RCPP
}
// subsetAssignGibbs
NumericMatrix subsetAssignGibbs(NumericVector y, NumericVector prop, NumericVector N, NumericMatrix isingCoefs, NumericVector nullEta, NumericVector altEta, NumericMatrix covariance, int nsamp, int nSubsets, int keepEach, int intSampSize, NumericVector MHcoef, IntegerVector popInd, NumericVector unifVec, NumericVector normVec, NumericVector dispersion, bool betaDispersion, IntegerVector preAssignment, double randomAssignProb, NumericVector mprobs, double preAssignCoef);
RcppExport SEXP flowReMix_subsetAssignGibbs(SEXP ySEXP, SEXP propSEXP, SEXP NSEXP, SEXP isingCoefsSEXP, SEXP nullEtaSEXP, SEXP altEtaSEXP, SEXP covarianceSEXP, SEXP nsampSEXP, SEXP nSubsetsSEXP, SEXP keepEachSEXP, SEXP intSampSizeSEXP, SEXP MHcoefSEXP, SEXP popIndSEXP, SEXP unifVecSEXP, SEXP normVecSEXP, SEXP dispersionSEXP, SEXP betaDispersionSEXP, SEXP preAssignmentSEXP, SEXP randomAssignProbSEXP, SEXP mprobsSEXP, SEXP preAssignCoefSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prop(propSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type isingCoefs(isingCoefsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nullEta(nullEtaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type altEta(altEtaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covariance(covarianceSEXP);
    Rcpp::traits::input_parameter< int >::type nsamp(nsampSEXP);
    Rcpp::traits::input_parameter< int >::type nSubsets(nSubsetsSEXP);
    Rcpp::traits::input_parameter< int >::type keepEach(keepEachSEXP);
    Rcpp::traits::input_parameter< int >::type intSampSize(intSampSizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MHcoef(MHcoefSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type popInd(popIndSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type unifVec(unifVecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type normVec(normVecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dispersion(dispersionSEXP);
    Rcpp::traits::input_parameter< bool >::type betaDispersion(betaDispersionSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type preAssignment(preAssignmentSEXP);
    Rcpp::traits::input_parameter< double >::type randomAssignProb(randomAssignProbSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mprobs(mprobsSEXP);
    Rcpp::traits::input_parameter< double >::type preAssignCoef(preAssignCoefSEXP);
    rcpp_result_gen = Rcpp::wrap(subsetAssignGibbs(y, prop, N, isingCoefs, nullEta, altEta, covariance, nsamp, nSubsets, keepEach, intSampSize, MHcoef, popInd, unifVec, normVec, dispersion, betaDispersion, preAssignment, randomAssignProb, mprobs, preAssignCoef));
    return rcpp_result_gen;
END_RCPP
}
// randomEffectCoordinateMH
NumericMatrix randomEffectCoordinateMH(NumericVector y, NumericVector N, int i, int nsamp, int nSubsets, NumericVector MHcoef, IntegerVector assignment, IntegerVector popInd, NumericVector eta, NumericVector randomEst, NumericVector condvar, NumericMatrix covariance, NumericMatrix invcov, NumericVector MHattempts, NumericVector MHsuccess, NumericVector unifVec, NumericVector dispersion, bool betaDispersion, int keepEach);
RcppExport SEXP flowReMix_randomEffectCoordinateMH(SEXP ySEXP, SEXP NSEXP, SEXP iSEXP, SEXP nsampSEXP, SEXP nSubsetsSEXP, SEXP MHcoefSEXP, SEXP assignmentSEXP, SEXP popIndSEXP, SEXP etaSEXP, SEXP randomEstSEXP, SEXP condvarSEXP, SEXP covarianceSEXP, SEXP invcovSEXP, SEXP MHattemptsSEXP, SEXP MHsuccessSEXP, SEXP unifVecSEXP, SEXP dispersionSEXP, SEXP betaDispersionSEXP, SEXP keepEachSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type nsamp(nsampSEXP);
    Rcpp::traits::input_parameter< int >::type nSubsets(nSubsetsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MHcoef(MHcoefSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type assignment(assignmentSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type popInd(popIndSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type randomEst(randomEstSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type condvar(condvarSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covariance(covarianceSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type invcov(invcovSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MHattempts(MHattemptsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MHsuccess(MHsuccessSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type unifVec(unifVecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dispersion(dispersionSEXP);
    Rcpp::traits::input_parameter< bool >::type betaDispersion(betaDispersionSEXP);
    Rcpp::traits::input_parameter< int >::type keepEach(keepEachSEXP);
    rcpp_result_gen = Rcpp::wrap(randomEffectCoordinateMH(y, N, i, nsamp, nSubsets, MHcoef, assignment, popInd, eta, randomEst, condvar, covariance, invcov, MHattempts, MHsuccess, unifVec, dispersion, betaDispersion, keepEach));
    return rcpp_result_gen;
END_RCPP
}
