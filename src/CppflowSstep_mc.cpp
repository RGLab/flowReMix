#include  <RcppArmadillo.h>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]
#include <omp.h>
#include "flowReMix.h"




subjectResult CppFlowSstep_mc(const SubjectDat subjectData, const int nsamp, const int nSubsets, const int intSampSize,
                  const arma::mat isingCoefs, const arma::mat covariance, const int keepEach, const arma::vec MHcoef,
                  const bool betaDispersion, const double randomAssignProb, const arma::vec modelprobs,
                  const double iterAssignCoef, const double prior, const bool zeroPosteriorProbs,
                  const arma::vec M, const arma::mat invcov, const bool mixed, const bool sampleRandom,
                  const arma::vec doNotSample,const  bool markovChainEM, const int msize) {
  try{
    if(doNotSample.size() != nSubsets){
      throw std::domain_error("doNotSample argument has invalid length!");
    }
    arma::vec condvar = (1.0)/(invcov.diag());

    arma::vec unifVec = flowReMix::myrunif(nsamp*nSubsets);
    arma::vec normVec = flowReMix::myrnorm(intSampSize);
    arma::mat assignmentMat;
    if(mixed){
      arma::mat assignmentMat(1.0,nSubsets);
      std::fill(assignmentMat.begin(),assignmentMat.end(),1.0);
    }else{

      assignmentMat = subsetAssignGibbs_mc( subjectData.y,  subjectData.prop,  subjectData.N,
                                         isingCoefs,
                                         subjectData.nullEta,  subjectData.altEta,
                                         covariance,
                                         nsamp,  nSubsets,  keepEach,  intSampSize,
                                         MHcoef,
                                         subjectData.popInd,
                                         unifVec,  normVec,
                                         M,  betaDispersion,
                                         subjectData.preAssignment,
                                         randomAssignProb,
                                         modelprobs,  iterAssignCoef,
                                         prior,  zeroPosteriorProbs,
                                         doNotSample,
                                         subjectData.subjassign, msize);


    }
    unifVec = flowReMix::myrunif(nsamp * nSubsets);
    arma::vec newEta(subjectData.altEta.size());
    newEta = subjectData.nullEta;//should copy
    arma::rowvec assignment = assignmentMat.row(assignmentMat.n_rows-1);

    std::vector<int> responderSubset = flowReMix::match2(subjectData.popInd,flowReMix::which2(assignment));

    // newEta[responderSubset] = altEta[responderSubset];
    flowReMix::mapByVec(newEta,subjectData.altEta,responderSubset);
    arma::vec randomEst = subjectData.rand; // randomEst will  be modified?

    arma::vec MHattempts(nSubsets);
    arma::vec MHsuccess(nSubsets);
    arma::mat randomMat;

    if(sampleRandom) {
      int index = subjectData.index;
      randomMat = simRandomEffectCoordinateMH_mc(subjectData.y, subjectData.N, subjectData.index,
                                              nsamp, nSubsets, MHcoef,
                                              assignment, subjectData.popInd, newEta,
                                              randomEst, condvar, covariance, invcov,
                                              MHattempts, MHsuccess, unifVec,
                                              M, betaDispersion, keepEach, msize);
    } else {
      //TODO set randomMat to zero, resize etc.
    }
    arma::vec rate = MHsuccess/MHattempts;

    subjectResult retval{assignmentMat,randomMat,rate};
    return(retval);
  } catch (std::exception &ex){
    forward_exception_to_r(ex);
  } catch (...){
    Rf_error("C++ exception (unknown reason)");
  }
  subjectResult ret{};
  return ret;
}





/*** R
rv144 = readRDS("~/Dropbox/GoTeam/Projects/flowReMix/rv144_aggregate_models.rds")[[2]]
# debug(flowReMix::flowSstep)
stabilityGraph(rv144,sampleNew = TRUE,cpus = 6,reps = 10)
*/
