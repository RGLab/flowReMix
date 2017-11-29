#include <RcppArmadillo.h>
#include <boost/throw_exception.hpp>
#include <random>
#include <omp.h>
using namespace Rcpp;

#ifndef FLOWREMIX_H
#define FLOWREMIX_H


namespace ParallelNormalGenerator {
  static std::vector<std::mt19937> generatorlist;
static void initialize(int t){
  for(int i=0;i<t;i++){
    std::mt19937 g(i);
    generatorlist.push_back(g);
  }
}
static double generate(double mean,double sigma){
  std::normal_distribution<double> d(mean,sigma);
  return d(generatorlist.at(omp_get_thread_num()));
}
}

namespace ParallelUnifGenerator {
static std::vector<std::mt19937> generatorlist;
static void initialize(int t){
  for(int i=0;i<t;i++){
    std::mt19937 g(i);
    generatorlist.push_back(g);
  }
}
static double generate(double lower,double upper){
  std::uniform_real_distribution<double> d(lower,upper);
  return d(generatorlist.at(omp_get_thread_num()));
}
}
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


arma::mat subsetAssignGibbs_mc(const arma::vec y,
                               const arma::vec prop,const arma::vec N,
                               const arma::mat isingCoefs,const arma::vec nullEta,
                               const arma::vec altEta, const arma::mat covariance,
                               int nsamp, const int nSubsets,
                               const int keepEach, const int intSampSize,
                               const arma::vec MHcoef, const arma::vec popInd,
                               const arma::vec unifVec, const arma::vec normVec,
                               const arma::vec dispersion, const bool betaDispersion,
                               const arma::vec preAssignment, const double randomAssignProb,
                               const arma::vec mprobs, const double preAssignCoef,
                               const double prior, const bool zeroPosteriorProbs,
                               const arma::vec doNotSample, arma::vec assignment,
                               const int msize);

arma::mat simRandomEffectCoordinateMH_mc(const arma::vec y, const arma::vec N,
                                         const int i, int nsamp, const int nSubsets,
                                         const arma::vec MHcoef, const arma::vec assignment,
                                         const arma::vec popInd, const arma::vec eta,
                                         const arma::vec randomEstt, const arma::vec condvar,
                                         const arma::mat covariance, const arma::mat invcov,
                                         arma::vec& MHattempts, arma::vec& MHsuccess,
                                         const arma::vec unifVec, const arma::vec dispersion,
                                         const bool betaDispersion, const int keepEach,
                                         const int msize);
double expit(double);
double betaBinomDens(int count, int N, double prob, double M);

subjectResult CppFlowSstep_mc(const SubjectDat subjectData, const int nsamp, const int nSubsets, const int intSampSize,
                              const arma::mat isingCoefs, const arma::mat covariance, const int keepEach, const arma::vec MHcoef,
                              const bool betaDispersion, const double randomAssignProb, const arma::vec modelprobs,
                              const double iterAssignCoef, const double prior, const bool zeroPosteriorProbs,
                              const arma::vec M, const arma::mat invcov, const bool mixed, const bool sampleRandom,
                              const arma::vec doNotSample,const  bool markovChainEM, const int msize);

arma::mat computeRandomEta_arma(arma::vec eta, arma::vec vsample);
arma::mat computeRandomEta_arma(arma::mat eta, arma::vec vsample);
arma::vec computeIntegratedDensities_arma(arma::mat logdens);
inline arma::vec computeBinomDensity_arma(const arma::vec subsetCount,
                                   const arma::vec subsetN,
                                   const arma::mat randomEta,
                                   const double M,
                                   const bool betaDispersion);
arma::vec logit_arma(arma::vec p);

void printDims(arma::vec a,std::string c);
void printDims(arma::mat a,std::string c);
void printDims(arma::uvec a,std::string c);

double binomDensityForMH_arma(arma::vec count, arma::vec N,
                              arma::vec eta, double proposal,
                              double M, bool betaDispersion);
double computeConditionalMean_arma(int subset,
                                   arma::vec condvar,
                                   arma::mat invcov,
                                   arma::vec randomEst);

List CppFlowSstepList_mc_vec(const int nsubjects,const arma::mat Y,
                             const arma::mat N,const arma::mat subpopInd,
                             const arma::mat nullEta, const arma::mat altEta,
                             const arma::mat rand, arma::vec index, arma::mat preassign,
                             const int nsamp,const int nSubsets,
                             int intSampSize , const arma::mat isingCoefs,
                             const arma::mat covariance,const int keepEach,
                             const arma::vec MHcoef, bool betaDispersion,
                             const double randomAssignProb,
                             const double iterAssignCoef, const double prior,
                             const bool zeroPosteriorProbs,
                             const arma::vec M, const arma::mat invcov, const bool mixed,
                             const bool sampleRandom,
                             const arma::vec doNotSample,const bool markovChainEM);


double lbeta_cpp(double a, double b);
double lchoose_cpp(double n, double k);

namespace flowReMix
{

// lambdas for the flowReMix namespace
auto myrnorm = [](int s){
  arma::vec r(s);
  std::transform(r.begin(),r.end(),r.begin(),[](auto a){return (double) R::rnorm(0,1);});
  return(r);
};
auto myrunif = [](int s){
  arma::vec r(s);
  std::transform(r.begin(),r.end(),r.begin(),[](auto a){
    return ParallelUnifGenerator::generate(0,1);
    // return (double) R::runif(0,1);
    });
  return(r);
};





auto myrnorm3 = [](int s, double mean, double sigma){
  arma::vec r(s);
  std::transform(r.begin(),r.end(),r.begin(),[&mean,&sigma](auto a){
    // return (double) R::rnorm(mean,sigma);
    return ParallelNormalGenerator::generate(mean,sigma);
    });
  return(r);
};


auto myrnorm2 = [](int rows, int cols){
  arma::mat r(rows,cols);
  std::transform(r.begin(),r.end(),r.begin(),[](auto a){return (double) R::rnorm(0,1);});
  return(r);
};
auto myrunif2 = [](int rows, int cols){
  arma::mat r(rows,cols);
  std::transform(r.begin(),r.end(),r.begin(),[](auto a){return (double) R::runif(0,1);});
  return(r);
};

auto which2 = [](arma::rowvec vec){
  std::vector<double> result;
  for(int i=0; i<vec.size();++i){
    if(vec(i)==1.0){
      result.push_back(i+1.0);
    }
  }
  return(result);
};

auto find = [] (arma::vec pi,int subset){
  std::vector<int> ind;
  for(int i = 0; i<pi.size();i++){
    if(pi(i)==subset){
      ind.push_back(i);;
    }
  }
  arma::uvec ind2(ind.size());
  std::copy(ind.begin(),ind.end(),ind2.begin());
  return(ind2);
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
auto mapByVec = [] (auto &dest,auto source,auto inds){
  for(int i=0;i<dest.size();i++){
    if(inds.at(i)==1){
      dest(i)=source(i);
    }
  }
};


auto pmax = [] (const arma::vec& a,const arma::vec& b) {
  auto max = [](double a, double b){
    if(a>b){
      return(a);
    }else{
      return b;
    }
  };

  arma::vec mx(a.size());
  for(int i=0;i<a.size();i++){
    mx(i) = max(a(i),b(i));
  }
  return(mx);
};

auto dnorm4 = [] (const arma::vec& a,const double& b, const double& c, const bool& islog) -> arma::vec {
  arma::vec result(a.size());
  for(int i=0;i<a.size();i++){
    result.at(i) = R::dnorm4(a.at(i),b,c,islog);
  }
  return(result);
};


}

#endif /* FLOWREMIX_H */

