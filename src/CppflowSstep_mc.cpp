#include  <RcppArmadillo.h>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
// [[Rcpp::plugins(openmp)]]
#include <omp.h>

extern arma::mat subsetAssignGibbs_mc(const arma::vec y, const arma::vec prop, const arma::vec N, const arma::mat isingCoefs, const arma::vec nullEta, const arma::vec altEta, const arma::mat covariance, int nsamp, const int nSubsets, const int keepEach, const int intSampSize, const arma::vec MHcoef, const arma::vec popInd, const arma::vec unifVec, const arma::vec normVec, const arma::vec dispersion, const bool betaDispersion, const arma::vec preAssignment, const double randomAssignProb, const arma::vec mprobs, const double preAssignCoef, const double prior, const bool zeroPosteriorProbs, const arma::vec doNotSample, arma::vec assignment,const int msize);
extern arma::mat simRandomEffectCoordinateMH_mc(const arma::vec y, const arma::vec N, const int i, int nsamp, const int nSubsets, const arma::vec MHcoef, const arma::vec assignment, const arma::vec popInd, const arma::vec eta, const arma::vec randomEstt, const arma::vec condvar, const arma::mat covariance, const arma::mat invcov, arma::vec MHattempts, arma::vec MHsuccess, const arma::vec unifVec, const arma::vec dispersion, const bool betaDispersion, const int keepEach, const int msize);
extern double expit(double);
extern double betaBinomDens(int count, int N, double prob, double M);


auto printDims(arma::vec a,std::string c){
  Rcout<<c<<" "<<a.n_rows<<" x "<<a.n_cols<<"\n";
};

auto printDims(arma::mat a,std::string c){
  Rcout<<c<<" "<<a.n_rows<<" x "<<a.n_cols<<"\n";
};

class SubjectDat {
public:
  arma::vec rand, N, y, nullEta, altEta, popInd, preAssignment, prop, subjassign;
  int index;
  SubjectDat(const List& l,const int &nSubsets,const bool &markovChainEM){
    DataFrame dat = (SEXP) l["dat"];
    DataFrame pre = (SEXP) l["pre"];

    NumericVector _rand = (SEXP) l["rand"];
    rand = arma::vec(_rand.begin(),_rand.length(),false);

    int index = as<int>(l["index"]);

    NumericVector _N = (SEXP) dat["N"];
    N = arma::vec(_N.begin(),_N.length(),false);

    NumericVector _y = (SEXP) dat["y"];
    y = arma::vec(_y.begin(),_y.length(),false);

    NumericVector _nullEta = (SEXP) dat["nullEta"];
    nullEta = arma::vec(_nullEta.begin(),_nullEta.length(),false);

    NumericVector _altEta = (SEXP) dat["altEta"];
    altEta = arma::vec(_altEta.begin(),_altEta.length(),false);

    IntegerVector _popInd = (SEXP) dat["subpopInd"];
    popInd = arma::vec(_popInd.length());
    std::copy(_popInd.begin(),_popInd.end(),popInd.begin());

    IntegerVector _preAssignment = (SEXP) pre["assign"];
    preAssignment = arma::vec(_preAssignment.length());
    std::copy(_preAssignment.begin(),_preAssignment.end(),preAssignment.begin());


    index = as<IntegerVector>(l["index"]).at(0);

    prop = y/N;

    if(markovChainEM || ((CharacterVector) l.names()).cend() == std::find(((CharacterVector) l.names()).cbegin(),((CharacterVector) l.names()).cend(),"assign")){
      subjassign = arma::vec(nSubsets);
    }else{
      subjassign = arma::vec(((NumericVector)(SEXP) l["assign"]));
    }
  };
};

struct subjectResult{
  arma::mat assignmentMat;
  arma::mat randomMat;
  arma::vec rate;
};

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

    auto myrnorm = [](int s){
      arma::vec r(s);
      std::transform(r.begin(),r.end(),r.begin(),[](auto a){return (double) R::rnorm(0,1);});
      return(r);
    };
    auto myrunif = [](int s){
      arma::vec r(s);
      std::transform(r.begin(),r.end(),r.begin(),[](auto a){return (double) R::runif(0,1);});
      return(r);
    };
    arma::vec unifVec = myrunif(nsamp*nSubsets);
    arma::vec normVec = myrnorm(intSampSize);
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
    unifVec = myrunif(nsamp * nSubsets);
    arma::vec newEta(subjectData.altEta.size());
    newEta = subjectData.nullEta;//should copy
    arma::rowvec assignment = assignmentMat.row(assignmentMat.n_rows-1);
    auto which2 = [](arma::rowvec vec){
      std::vector<double> result;
      for(int i=0; i<vec.size();++i){
        if(vec(i)==1.0){
          result.push_back(i+1.0);
        }
      }
      return(result);
    };
    auto match2 = [] (arma::vec v1, std::vector<double> v2){
      std::vector<int> res(v1.size());
      for(int i=0;i<v1.size();++i){
        double ati = v1.at(i);
        for(int j=0;j<v2.size();++j){
          if(v2.at(j)==ati){
            res.at(i) = 1;
          }
        }
      }
      return(res);
    };
    std::vector<int> responderSubset = match2(subjectData.popInd,which2(assignment));
    auto mapByVec = [] (auto &dest,auto source,auto inds){
      for(int i=0;i<dest.size();i++){
        if(inds.at(i)==1){
          dest.at(i)=source.at(i);
        }
      }
    };
    // newEta[responderSubset] = altEta[responderSubset];
    mapByVec(newEta,subjectData.altEta,responderSubset);
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




/*** R
rv144 = readRDS("~/Dropbox/GoTeam/Projects/flowReMix/rv144_aggregate_models.rds")[[2]]
# debug(flowReMix::flowSstep)
stabilityGraph(rv144,sampleNew = TRUE,cpus = 6,reps = 10)
*/
