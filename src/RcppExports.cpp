// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// CppFlowSstepList
List CppFlowSstepList(const List& subjectDataList, int nsamp, const int nSubsets, const int intSampSize, const NumericMatrix& isingCoefs, const NumericMatrix& covariance, const int keepEach, const NumericVector& MHcoef, const bool betaDispersion, const double randomAssignProb, const NumericVector& modelprobs, const double iterAssignCoef, const double prior, const bool zeroPosteriorProbs, const NumericVector& M, const NumericMatrix& invcov, const bool mixed, const bool sampleRandom, const LogicalVector& doNotSample, const bool markovChainEM);
RcppExport SEXP _flowReMix_CppFlowSstepList(SEXP subjectDataListSEXP, SEXP nsampSEXP, SEXP nSubsetsSEXP, SEXP intSampSizeSEXP, SEXP isingCoefsSEXP, SEXP covarianceSEXP, SEXP keepEachSEXP, SEXP MHcoefSEXP, SEXP betaDispersionSEXP, SEXP randomAssignProbSEXP, SEXP modelprobsSEXP, SEXP iterAssignCoefSEXP, SEXP priorSEXP, SEXP zeroPosteriorProbsSEXP, SEXP MSEXP, SEXP invcovSEXP, SEXP mixedSEXP, SEXP sampleRandomSEXP, SEXP doNotSampleSEXP, SEXP markovChainEMSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type subjectDataList(subjectDataListSEXP);
    Rcpp::traits::input_parameter< int >::type nsamp(nsampSEXP);
    Rcpp::traits::input_parameter< const int >::type nSubsets(nSubsetsSEXP);
    Rcpp::traits::input_parameter< const int >::type intSampSize(intSampSizeSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type isingCoefs(isingCoefsSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type covariance(covarianceSEXP);
    Rcpp::traits::input_parameter< const int >::type keepEach(keepEachSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type MHcoef(MHcoefSEXP);
    Rcpp::traits::input_parameter< const bool >::type betaDispersion(betaDispersionSEXP);
    Rcpp::traits::input_parameter< const double >::type randomAssignProb(randomAssignProbSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type modelprobs(modelprobsSEXP);
    Rcpp::traits::input_parameter< const double >::type iterAssignCoef(iterAssignCoefSEXP);
    Rcpp::traits::input_parameter< const double >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const bool >::type zeroPosteriorProbs(zeroPosteriorProbsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type invcov(invcovSEXP);
    Rcpp::traits::input_parameter< const bool >::type mixed(mixedSEXP);
    Rcpp::traits::input_parameter< const bool >::type sampleRandom(sampleRandomSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type doNotSample(doNotSampleSEXP);
    Rcpp::traits::input_parameter< const bool >::type markovChainEM(markovChainEMSEXP);
    rcpp_result_gen = Rcpp::wrap(CppFlowSstepList(subjectDataList, nsamp, nSubsets, intSampSize, isingCoefs, covariance, keepEach, MHcoef, betaDispersion, randomAssignProb, modelprobs, iterAssignCoef, prior, zeroPosteriorProbs, M, invcov, mixed, sampleRandom, doNotSample, markovChainEM));
    return rcpp_result_gen;
END_RCPP
}
// CppFlowSstepList_mc_vec
List CppFlowSstepList_mc_vec(const int nsubjects, const arma::mat& Y, const arma::mat& N, const arma::mat& subpopInd, arma::mat& clusterassignments, const arma::mat& nullEta, const arma::mat altEta, const arma::mat& rand, arma::vec& index, const arma::mat& preassign, const int nsamp, const int nsubsets, const int intSampSize, const arma::mat& isingCoefs, const arma::mat& covariance, const int keepEach, const arma::vec& MHcoef, const bool betaDispersion, const double randomAssignProb, const double iterAssignCoef, const double prior, const bool zeroPosteriorProbs, const arma::vec& M, const arma::mat& invcov, const bool mixed, const bool sampleRandom, const arma::vec& doNotSample, const bool markovChainEM, int cpus, int seed);
RcppExport SEXP _flowReMix_CppFlowSstepList_mc_vec(SEXP nsubjectsSEXP, SEXP YSEXP, SEXP NSEXP, SEXP subpopIndSEXP, SEXP clusterassignmentsSEXP, SEXP nullEtaSEXP, SEXP altEtaSEXP, SEXP randSEXP, SEXP indexSEXP, SEXP preassignSEXP, SEXP nsampSEXP, SEXP nsubsetsSEXP, SEXP intSampSizeSEXP, SEXP isingCoefsSEXP, SEXP covarianceSEXP, SEXP keepEachSEXP, SEXP MHcoefSEXP, SEXP betaDispersionSEXP, SEXP randomAssignProbSEXP, SEXP iterAssignCoefSEXP, SEXP priorSEXP, SEXP zeroPosteriorProbsSEXP, SEXP MSEXP, SEXP invcovSEXP, SEXP mixedSEXP, SEXP sampleRandomSEXP, SEXP doNotSampleSEXP, SEXP markovChainEMSEXP, SEXP cpusSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type nsubjects(nsubjectsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type subpopInd(subpopIndSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type clusterassignments(clusterassignmentsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type nullEta(nullEtaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type altEta(altEtaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type rand(randSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type index(indexSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type preassign(preassignSEXP);
    Rcpp::traits::input_parameter< const int >::type nsamp(nsampSEXP);
    Rcpp::traits::input_parameter< const int >::type nsubsets(nsubsetsSEXP);
    Rcpp::traits::input_parameter< const int >::type intSampSize(intSampSizeSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type isingCoefs(isingCoefsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type covariance(covarianceSEXP);
    Rcpp::traits::input_parameter< const int >::type keepEach(keepEachSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type MHcoef(MHcoefSEXP);
    Rcpp::traits::input_parameter< const bool >::type betaDispersion(betaDispersionSEXP);
    Rcpp::traits::input_parameter< const double >::type randomAssignProb(randomAssignProbSEXP);
    Rcpp::traits::input_parameter< const double >::type iterAssignCoef(iterAssignCoefSEXP);
    Rcpp::traits::input_parameter< const double >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const bool >::type zeroPosteriorProbs(zeroPosteriorProbsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type invcov(invcovSEXP);
    Rcpp::traits::input_parameter< const bool >::type mixed(mixedSEXP);
    Rcpp::traits::input_parameter< const bool >::type sampleRandom(sampleRandomSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type doNotSample(doNotSampleSEXP);
    Rcpp::traits::input_parameter< const bool >::type markovChainEM(markovChainEMSEXP);
    Rcpp::traits::input_parameter< int >::type cpus(cpusSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(CppFlowSstepList_mc_vec(nsubjects, Y, N, subpopInd, clusterassignments, nullEta, altEta, rand, index, preassign, nsamp, nsubsets, intSampSize, isingCoefs, covariance, keepEach, MHcoef, betaDispersion, randomAssignProb, iterAssignCoef, prior, zeroPosteriorProbs, M, invcov, mixed, sampleRandom, doNotSample, markovChainEM, cpus, seed));
    return rcpp_result_gen;
END_RCPP
}
// vecBetaBinomDens
NumericVector vecBetaBinomDens(NumericVector count, NumericVector N, NumericVector prob, double M);
RcppExport SEXP _flowReMix_vecBetaBinomDens(SEXP countSEXP, SEXP NSEXP, SEXP probSEXP, SEXP MSEXP) {
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
RcppExport SEXP _flowReMix_setNumericVectorToZero(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    setNumericVectorToZero(x);
    return R_NilValue;
END_RCPP
}
// subsetAssignGibbs
NumericMatrix subsetAssignGibbs(const NumericVector& y, const NumericVector& prop, const NumericVector& N, const NumericMatrix& isingCoefs, const NumericVector& nullEta, const NumericVector& altEta, const NumericMatrix& covariance, int nsamp, const int nSubsets, const int keepEach, const int intSampSize, const NumericVector& MHcoef, const IntegerVector& popInd, const NumericVector& unifVec, const NumericVector& normVec, const NumericVector& dispersion, const bool betaDispersion, const IntegerVector& preAssignment, const double randomAssignProb, const NumericVector& mprobs, const double preAssignCoef, const double prior, const bool zeroPosteriorProbs, const LogicalVector& doNotSample, NumericVector assignment, const int msize);
RcppExport SEXP _flowReMix_subsetAssignGibbs(SEXP ySEXP, SEXP propSEXP, SEXP NSEXP, SEXP isingCoefsSEXP, SEXP nullEtaSEXP, SEXP altEtaSEXP, SEXP covarianceSEXP, SEXP nsampSEXP, SEXP nSubsetsSEXP, SEXP keepEachSEXP, SEXP intSampSizeSEXP, SEXP MHcoefSEXP, SEXP popIndSEXP, SEXP unifVecSEXP, SEXP normVecSEXP, SEXP dispersionSEXP, SEXP betaDispersionSEXP, SEXP preAssignmentSEXP, SEXP randomAssignProbSEXP, SEXP mprobsSEXP, SEXP preAssignCoefSEXP, SEXP priorSEXP, SEXP zeroPosteriorProbsSEXP, SEXP doNotSampleSEXP, SEXP assignmentSEXP, SEXP msizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type prop(propSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type isingCoefs(isingCoefsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type nullEta(nullEtaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type altEta(altEtaSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type covariance(covarianceSEXP);
    Rcpp::traits::input_parameter< int >::type nsamp(nsampSEXP);
    Rcpp::traits::input_parameter< const int >::type nSubsets(nSubsetsSEXP);
    Rcpp::traits::input_parameter< const int >::type keepEach(keepEachSEXP);
    Rcpp::traits::input_parameter< const int >::type intSampSize(intSampSizeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type MHcoef(MHcoefSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type popInd(popIndSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type unifVec(unifVecSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type normVec(normVecSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type dispersion(dispersionSEXP);
    Rcpp::traits::input_parameter< const bool >::type betaDispersion(betaDispersionSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type preAssignment(preAssignmentSEXP);
    Rcpp::traits::input_parameter< const double >::type randomAssignProb(randomAssignProbSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mprobs(mprobsSEXP);
    Rcpp::traits::input_parameter< const double >::type preAssignCoef(preAssignCoefSEXP);
    Rcpp::traits::input_parameter< const double >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< const bool >::type zeroPosteriorProbs(zeroPosteriorProbsSEXP);
    Rcpp::traits::input_parameter< const LogicalVector& >::type doNotSample(doNotSampleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type assignment(assignmentSEXP);
    Rcpp::traits::input_parameter< const int >::type msize(msizeSEXP);
    rcpp_result_gen = Rcpp::wrap(subsetAssignGibbs(y, prop, N, isingCoefs, nullEta, altEta, covariance, nsamp, nSubsets, keepEach, intSampSize, MHcoef, popInd, unifVec, normVec, dispersion, betaDispersion, preAssignment, randomAssignProb, mprobs, preAssignCoef, prior, zeroPosteriorProbs, doNotSample, assignment, msize));
    return rcpp_result_gen;
END_RCPP
}
// simRandomEffectCoordinateMH
NumericMatrix simRandomEffectCoordinateMH(const NumericVector& y, const NumericVector& N, const int i, int nsamp, const int nSubsets, const NumericVector& MHcoef, const IntegerVector& assignment, const IntegerVector& popInd, const NumericVector& eta, const NumericVector& randomEstt, const NumericVector& condvar, const NumericMatrix& covariance, const NumericMatrix& invcov, NumericVector MHattempts, NumericVector MHsuccess, const NumericVector& unifVec, const NumericVector& dispersion, const bool betaDispersion, const int keepEach, const int msize);
RcppExport SEXP _flowReMix_simRandomEffectCoordinateMH(SEXP ySEXP, SEXP NSEXP, SEXP iSEXP, SEXP nsampSEXP, SEXP nSubsetsSEXP, SEXP MHcoefSEXP, SEXP assignmentSEXP, SEXP popIndSEXP, SEXP etaSEXP, SEXP randomEsttSEXP, SEXP condvarSEXP, SEXP covarianceSEXP, SEXP invcovSEXP, SEXP MHattemptsSEXP, SEXP MHsuccessSEXP, SEXP unifVecSEXP, SEXP dispersionSEXP, SEXP betaDispersionSEXP, SEXP keepEachSEXP, SEXP msizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type nsamp(nsampSEXP);
    Rcpp::traits::input_parameter< const int >::type nSubsets(nSubsetsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type MHcoef(MHcoefSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type assignment(assignmentSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type popInd(popIndSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type randomEstt(randomEsttSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type condvar(condvarSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type covariance(covarianceSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type invcov(invcovSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MHattempts(MHattemptsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MHsuccess(MHsuccessSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type unifVec(unifVecSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type dispersion(dispersionSEXP);
    Rcpp::traits::input_parameter< const bool >::type betaDispersion(betaDispersionSEXP);
    Rcpp::traits::input_parameter< const int >::type keepEach(keepEachSEXP);
    Rcpp::traits::input_parameter< const int >::type msize(msizeSEXP);
    rcpp_result_gen = Rcpp::wrap(simRandomEffectCoordinateMH(y, N, i, nsamp, nSubsets, MHcoef, assignment, popInd, eta, randomEstt, condvar, covariance, invcov, MHattempts, MHsuccess, unifVec, dispersion, betaDispersion, keepEach, msize));
    return rcpp_result_gen;
END_RCPP
}
// newMHsampler
void newMHsampler(NumericMatrix assign, NumericMatrix random, NumericVector initAssign, NumericVector initRand, NumericVector y, NumericVector N, int keepEach, double prior, NumericMatrix isingCoefs, IntegerVector preAssignment, NumericMatrix invcov, NumericMatrix covariance, NumericVector condvar, NumericVector dispersion, NumericVector nullEta, NumericVector altEta, IntegerVector popInd, NumericVector MHattempts, NumericVector MHsuccess, NumericVector MHcoef);
RcppExport SEXP _flowReMix_newMHsampler(SEXP assignSEXP, SEXP randomSEXP, SEXP initAssignSEXP, SEXP initRandSEXP, SEXP ySEXP, SEXP NSEXP, SEXP keepEachSEXP, SEXP priorSEXP, SEXP isingCoefsSEXP, SEXP preAssignmentSEXP, SEXP invcovSEXP, SEXP covarianceSEXP, SEXP condvarSEXP, SEXP dispersionSEXP, SEXP nullEtaSEXP, SEXP altEtaSEXP, SEXP popIndSEXP, SEXP MHattemptsSEXP, SEXP MHsuccessSEXP, SEXP MHcoefSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type assign(assignSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type random(randomSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type initAssign(initAssignSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type initRand(initRandSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type keepEach(keepEachSEXP);
    Rcpp::traits::input_parameter< double >::type prior(priorSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type isingCoefs(isingCoefsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type preAssignment(preAssignmentSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type invcov(invcovSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covariance(covarianceSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type condvar(condvarSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dispersion(dispersionSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nullEta(nullEtaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type altEta(altEtaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type popInd(popIndSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MHattempts(MHattemptsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MHsuccess(MHsuccessSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type MHcoef(MHcoefSEXP);
    newMHsampler(assign, random, initAssign, initRand, y, N, keepEach, prior, isingCoefs, preAssignment, invcov, covariance, condvar, dispersion, nullEta, altEta, popInd, MHattempts, MHsuccess, MHcoef);
    return R_NilValue;
END_RCPP
}
// weightedMean
NumericVector weightedMean(NumericVector x, NumericVector weights, bool na_rm);
RcppExport SEXP _flowReMix_weightedMean(SEXP xSEXP, SEXP weightsSEXP, SEXP na_rmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    rcpp_result_gen = Rcpp::wrap(weightedMean(x, weights, na_rm));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_flowReMix_CppFlowSstepList", (DL_FUNC) &_flowReMix_CppFlowSstepList, 20},
    {"_flowReMix_CppFlowSstepList_mc_vec", (DL_FUNC) &_flowReMix_CppFlowSstepList_mc_vec, 30},
    {"_flowReMix_vecBetaBinomDens", (DL_FUNC) &_flowReMix_vecBetaBinomDens, 4},
    {"_flowReMix_setNumericVectorToZero", (DL_FUNC) &_flowReMix_setNumericVectorToZero, 1},
    {"_flowReMix_subsetAssignGibbs", (DL_FUNC) &_flowReMix_subsetAssignGibbs, 26},
    {"_flowReMix_simRandomEffectCoordinateMH", (DL_FUNC) &_flowReMix_simRandomEffectCoordinateMH, 20},
    {"_flowReMix_newMHsampler", (DL_FUNC) &_flowReMix_newMHsampler, 20},
    {"_flowReMix_weightedMean", (DL_FUNC) &_flowReMix_weightedMean, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_flowReMix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
