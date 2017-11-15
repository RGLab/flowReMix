#include  <RcppArmadillo.h>
#include <algorithm>
using namespace Rcpp;

extern NumericMatrix subsetAssignGibbs(const NumericVector& y, const NumericVector& prop, const NumericVector& N, const NumericMatrix& isingCoefs, const NumericVector& nullEta, const NumericVector& altEta, const NumericMatrix& covariance, int nsamp, const int nSubsets, const int keepEach, const int intSampSize, const NumericVector& MHcoef, const IntegerVector& popInd, const NumericVector& unifVec, const NumericVector& normVec, const NumericVector& dispersion, const bool betaDispersion, const IntegerVector& preAssignment, const double randomAssignProb, const NumericVector& mprobs, const double preAssignCoef, const double prior, const bool zeroPosteriorProbs, const LogicalVector& doNotSample, NumericVector assignment,const int msize);
extern NumericMatrix simRandomEffectCoordinateMH(const NumericVector& y, const NumericVector& N, const int i, int nsamp, const int nSubsets, const NumericVector& MHcoef, const IntegerVector& assignment, const IntegerVector& popInd, const NumericVector& eta, const NumericVector& randomEstt, const NumericVector& condvar, const NumericMatrix& covariance, const NumericMatrix& invcov, NumericVector MHattempts, NumericVector MHsuccess, const NumericVector& unifVec, const NumericVector& dispersion, const bool betaDispersion, const int keepEach, const int msize);


List CppFlowSstep(const List& subjectData, const int nsamp, const int nSubsets, const int intSampSize,
              const NumericMatrix& isingCoefs, const NumericMatrix& covariance, const int keepEach, const NumericVector& MHcoef,
              const bool betaDispersion, const double randomAssignProb, const NumericVector& modelprobs,
              const double iterAssignCoef, const double prior, const bool zeroPosteriorProbs,
              const NumericVector& M, const NumericMatrix& invcov, const bool mixed, const bool sampleRandom,
              const LogicalVector& doNotSample,const  bool markovChainEM, const int msize) {
  try{
    if(doNotSample.length() != nSubsets){
      throw std::domain_error("doNotSample argument has invalid length!");
    }
    NumericVector condvar(invcov.nrow());

    auto invdiag = [](const NumericMatrix &m, NumericVector &v){
      for(int i=0;i<m.nrow();++i){
        v(i) = 1.0/m(i,i);
      }
      return;
    };

    invdiag(invcov,condvar);

    DataFrame dat = (SEXP) subjectData["dat"];
    DataFrame pre = (SEXP) subjectData["pre"];
    NumericVector rand = (SEXP) subjectData["rand"];
    int index = as<int>(subjectData["index"]);
    NumericVector N = (SEXP) dat["N"];
    NumericVector y = (SEXP) dat["y"];
    NumericVector nullEta = (SEXP) dat["nullEta"];
    NumericVector altEta = (SEXP) dat["altEta"];
    IntegerVector popInd = (SEXP) dat["subpopInd"];
    NumericVector prop = y/N;
    NumericVector unifVec = Rcpp::runif(nsamp*nSubsets);
    NumericVector normVec = Rcpp::rnorm(intSampSize);
    IntegerVector preAssignment = (SEXP) pre["assign"];
    NumericMatrix assignmentMat;
    NumericVector subjassign;
    if(mixed){
      NumericMatrix assignmentMat(1.0,nSubsets);
      std::fill(assignmentMat.begin(),assignmentMat.end(),1.0);
    }else{
      // check if subjectData has an "assign" slot.
      if(markovChainEM || ((CharacterVector) subjectData.names()).cend() == std::find(((CharacterVector) subjectData.names()).cbegin(),((CharacterVector) subjectData.names()).cend(),"assign")){
          subjassign = NumericVector(nSubsets);
      }else{
         subjassign = (SEXP) subjectData["assign"];
      }
      assignmentMat = subsetAssignGibbs( y,  prop,  N,
                                                      isingCoefs,
                                                      nullEta,  altEta,
                                                      covariance,
                                                      nsamp,  nSubsets,  keepEach,  intSampSize,
                                                      MHcoef,
                                                      popInd,
                                                      unifVec,  normVec,
                                                      M,  betaDispersion,
                                                      preAssignment,
                                                      randomAssignProb,
                                                      modelprobs,  iterAssignCoef,
                                                      prior,  zeroPosteriorProbs,
                                                      doNotSample,
                                                      subjassign, msize);

    }
      unifVec = runif(nsamp * nSubsets);
      NumericVector newEta(altEta.length());
      newEta = nullEta;//should copy
      NumericVector assignmentRow = assignmentMat(assignmentMat.nrow()-1,_);
      IntegerVector assignment = as<IntegerVector>(assignmentRow);
      auto which2 = [](IntegerVector vec){
        std::vector<int> result;
        for(int i=0; i<vec.length();++i){
          if(vec(i)==1){
            result.push_back(i+1);
          }
        }
        IntegerVector v = wrap(result);
        return(v);
      };
      auto match2 = [] (IntegerVector v1, IntegerVector v2){
        IntegerVector vres = match(v1,v2);
        LogicalVector res(vres.length());
        for(int i=0;i<vres.length();++i){
          if(vres(i)>0){
            res(i)=true;
          }else{
            res(i)=false;
          }
        }
        return(res);
      };
      LogicalVector responderSubset = match2(popInd,which2(assignment));

      newEta[responderSubset] = altEta[responderSubset];

      NumericVector randomEst = rand; // randomEst will  be modified?

      NumericVector MHattempts(nSubsets);
      NumericVector MHsuccess(nSubsets);
      NumericMatrix randomMat;

      if(sampleRandom) {
        int index = (as<IntegerVector>(subjectData["index"])).at(0);
        randomMat = simRandomEffectCoordinateMH(y, N, index,
                                                nsamp, nSubsets, MHcoef,
                                                assignment, popInd, newEta,
                                                randomEst, condvar, covariance, invcov,
                                                MHattempts, MHsuccess, unifVec,
                                                M, betaDispersion, keepEach, msize);
      } else {
        randomMat = NILSXP;
      }
      NumericVector rate = MHsuccess/MHattempts;
      List retval = List::create(
        Named("assign") = assignmentMat,
        Named("rand") = randomMat,
        Named("rate") = rate
      );

      return(retval);
  } catch (std::exception &ex){
      forward_exception_to_r(ex);
  } catch (...){
    Rf_error("C++ exception (unknown reason)");
  }

  return(List());
}


//[[Rcpp::export]]
List CppFlowSstepList(const List& subjectDataList, int nsamp, const int nSubsets,const int intSampSize,
                      const NumericMatrix& isingCoefs, const NumericMatrix& covariance, const int keepEach, const NumericVector& MHcoef,
                      const bool betaDispersion,const double randomAssignProb,const NumericVector& modelprobs,
                      const double iterAssignCoef,  const double prior, const bool zeroPosteriorProbs,
                      const NumericVector& M, const NumericMatrix& invcov, const bool mixed, const bool sampleRandom,
                      const LogicalVector& doNotSample, const bool markovChainEM) {
  nsamp = floor(nsamp / keepEach) * keepEach ;
  int msize = int(nsamp/keepEach);

  List results(subjectDataList.length());
  for(int i=0;i<subjectDataList.length();++i){
    SEXP le = subjectDataList[i];
    List subjectData = le;
    List currentResult = CppFlowSstep(subjectData, nsamp, nSubsets,intSampSize,isingCoefs,covariance,
                                      keepEach,MHcoef,betaDispersion,  randomAssignProb,  modelprobs,
                                      iterAssignCoef,  prior,  zeroPosteriorProbs,
                                      M,  invcov,  mixed,  sampleRandom,
                                      doNotSample,  markovChainEM, msize);
    results[i] = currentResult;
  }
  return(results);
}
/*** R
rv144 = readRDS("~/Dropbox/GoTeam/Projects/flowReMix/rv144_aggregate_models.rds")[[2]]
# debug(flowReMix::flowSstep)
stabilityGraph(rv144,sampleNew = TRUE,cpus = 6,reps = 10)
*/
