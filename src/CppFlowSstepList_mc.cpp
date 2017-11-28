#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include "flowReMix.h"


//[[Rcpp::export]]
List CppFlowSstepList_mc(const List subjectDataList,  int nsamp,const int nSubsets,const int intSampSize,
                         const arma::mat isingCoefs,const  arma::mat covariance, const int keepEach,const  arma::vec MHcoef,
                         const  bool betaDispersion,const double randomAssignProb,const arma::vec modelprobs,
                         const double iterAssignCoef, const  double prior,const  bool zeroPosteriorProbs,
                         const arma::vec M,const  arma::mat invcov,const  bool mixed,const  bool sampleRandom,
                         const arma::vec doNotSample,const  bool markovChainEM) {


  nsamp = floor(nsamp / keepEach) * keepEach ;
  int msize = int(nsamp/keepEach);

  std::vector<SubjectDat> subjects;

  for(int i=0;i<subjectDataList.length();++i){
    const List subjectData = (SEXP) subjectDataList[i];
    SubjectDat s(subjectData,nSubsets,markovChainEM);
    subjects.push_back(s);
  }

  Progress prog(subjectDataList.length(), true);
  std::vector<subjectResult> result;

  std::vector<subjectResult> currentResult;
  for(int i=0;i<subjects.size();++i){
    subjectResult r = CppFlowSstep_mc(subjects[i], nsamp, nSubsets, intSampSize, isingCoefs, covariance,
                                      keepEach, MHcoef, betaDispersion,  randomAssignProb,  modelprobs,
                                      iterAssignCoef,  prior,  zeroPosteriorProbs,
                                      M,  invcov,  mixed,  sampleRandom,
                                      doNotSample,  markovChainEM, msize);

    currentResult.push_back(r);
    if(Progress::check_abort()){
      break;
    }
    prog.increment();
  }

  List resultsList(currentResult.size());
  std::transform(currentResult.begin(),currentResult.end(),resultsList.begin(),[] (subjectResult a){
    List retval = List::create(
      Named("assign") = wrap(a.assignmentMat),
      Named("rand") = wrap(a.randomMat),
      Named("rate") = wrap(a.rate)
    );
    return retval;
  });
  return(resultsList);
}
