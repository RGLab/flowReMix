#include <RcppArmadillo.h>
#include "flowReMix.h"
#include <exception>
#include <boost/exception/all.hpp>
// [[Rcpp::depends(RcppProgress)]]
//[[Rcpp::plugins(openmp)]]
#ifdef _OPENMP
#include <omp.h>
#endif
#include <progress.hpp>

using namespace Rcpp;

void printDims(arma::vec a,std::string c){
  Rcout<<c<<" "<<a.n_rows<<" x "<<a.n_cols<<"\n";
};
void printDims(arma::uvec a,std::string c){
  Rcout<<c<<" "<<a.n_rows<<" x "<<a.n_cols<<"\n";
};

void printDims(arma::mat a,std::string c){
  Rcout<<c<<" "<<a.n_rows<<" x "<<a.n_cols<<"\n";
};




//[[Rcpp::export]]
List CppFlowSstepList_mc_vec(const int nsubjects,const arma::mat& Y,
                        const arma::mat& N,const arma::mat& subpopInd, arma::mat& clusterassignments,
                        const arma::mat& nullEta, const arma::mat altEta,
                        const arma::mat& rand, arma::vec& index, const arma::mat& preassign,
                        const int nsamp,const int nsubsets,
                        const int intSampSize , const arma::mat& isingCoefs,
                        const arma::mat& covariance,const int keepEach,
                        const arma::vec& MHcoef, const bool betaDispersion,
                        const double randomAssignProb,
                        const double iterAssignCoef, const double prior,
                        const bool zeroPosteriorProbs,
                        const arma::vec& M, const arma::mat& invcov, const bool mixed,
                        const bool sampleRandom,
                        const arma::vec& doNotSample,const bool markovChainEM, int cpus, int seed){
  int nsamp_floor = floor(nsamp / keepEach) * keepEach ;
  int mat_size = int(nsamp_floor/keepEach);
  Progress prog(nsubjects,true);
  try{
    if(nsubjects != Y.n_cols){
      throw std::domain_error("nsubjects and dimensions of Y don't match!");
    }
    if(nsubsets != rand.n_rows){
      throw std::domain_error("nSubsets and dimensions of rand don't match!");
    }
    // allocate objects for return values;
    arma::cube assignmentMats(mat_size,nsubsets,nsubjects);
    arma::cube randomeffectMats(mat_size,nsubsets,nsubjects);
    arma::mat MHsuccessrates(nsubsets,nsubjects);
    arma::mat  proportions = Y/N; //should be a double since Y, N, and proportions are mat which is Mat<double>
    //compute the conditional variance
    arma::vec condvar = (1.0)/(invcov.diag());


    //preallocate and precompute random samples
    arma::mat unifVec = flowReMix::myrunif2(nsamp_floor*nsubsets,nsubjects);
    arma::mat normVec = flowReMix::myrnorm2(intSampSize,nsubjects);

    //some conditional operations
    if(mixed){
      std::fill(assignmentMats.begin(),assignmentMats.end(),1.0);
    }

    //Gibbs for each subject
    arma::cube clusterDensities(intSampSize,2,nsubjects) ;
    arma::mat iterPosteriors(nsubsets,nsubjects) ;
    ParallelNormalGenerator::initialize(cpus, seed);
    ParallelUnifGenerator::initialize(cpus, seed);

#pragma omp parallel shared(unifVec,nsamp_floor,normVec,clusterassignments,proportions, assignmentMats,prog, preassign, clusterDensities, flowReMix::dnorm4, flowReMix::pmax, flowReMix::myrnorm3) num_threads(cpus)
{
  int abort = 0;

#pragma omp for
  for(int subject=0;subject<nsubjects;++subject){
      int unifPosition = 0 ;
      double isingOffset = 0 ;
      int assignNum = 0;
      arma::uvec subject_indicator(1);
      subject_indicator(0) = subject;

      int sample = 0, subset=0;

      for( sample = 0; sample < nsamp_floor ; sample++){
        for( subset = 0; subset < nsubsets ; subset++){
          if(!abort){
          if(doNotSample(subset) == 1.0 || (preassign(subset,subject) != -1 &  iterAssignCoef < 10e-4 & !zeroPosteriorProbs)) {
            clusterassignments(subset,subject) = 0 ;
            continue ;
          } else if(preassign(subset,subject) != -1 & iterAssignCoef < 10e-4 & !zeroPosteriorProbs) {
            clusterassignments(subset,subject) = preassign(subset,subject) ;
            continue;
          }else if(preassign(subset,subject) == 0 ) {
            isingOffset = -prior ;
          } else if(preassign(subset,subject) == 1) {
            isingOffset = prior ;
          } else {
            isingOffset = 0 ;
          }
          // isingOffset = 0 ;
            arma::uvec subset_indicator = flowReMix::find(subpopInd.col(subject),subset+1);
          if(subset_indicator.size() == 0) {
            continue ;
          }
          arma::mat subsetNullEta, subsetAltEta, subsetN, subsetProp, subsetCount;
          subsetNullEta = nullEta(subset_indicator,subject_indicator);
          subsetAltEta = altEta(subset_indicator,subject_indicator);
          subsetN = N(subset_indicator,subject_indicator);
          subsetProp = proportions(subset_indicator,subject_indicator)  ;
          subsetCount = Y(subset_indicator,subject_indicator) ;

          // Some cell populations for a specific subject may be missing
          double sigmaHat = sqrt(covariance(subset, subset)) ;
          arma::vec empEta = logit_arma(flowReMix::pmax(subsetProp, 1.0 / subsetN)) ;

          // integrating densities

          double mcoef;

          mcoef = MHcoef(subset) ;

          mcoef = std::max(1.0, mcoef) ;
          double muHat = 0;
          arma::vec vsample;
          double prevMuHat;
          for(int cluster = 0; cluster < 2; cluster++) {
            arma::mat eta;
            if(cluster == 0) {
              eta = subsetNullEta ;
            } else {
              eta = subsetAltEta ;
            }
            // arma::vec etaResid = empEta - eta;

            //the next conditions are totally pointless.
            if(cluster == 0) {
              // muHat = mean(etaResid) ; //this is removed because we are just doing random walk.
              muHat = 0;

              vsample = flowReMix::myrnorm3(intSampSize, muHat, sigmaHat * mcoef) ; //vsample is fixed across clusters.

            }
            arma::vec importanceWeights;

            arma::vec sampNormDens = flowReMix::dnorm4(vsample, muHat, sigmaHat * mcoef, true) ;
            arma::vec normDens = flowReMix::dnorm4(vsample, 0, sigmaHat, true) ; // presumably we're tuning mcoef.
            importanceWeights = normDens - sampNormDens ;

            arma::mat randomEta;

            randomEta = computeRandomEta_arma(eta, vsample) ;

            double disp;

            disp = M(subset);

            // arma::vec binomDensity = computeBinomDensity_arma(subsetCount, subsetN, randomEta, disp, betaDispersion) ;

            int sampSize = randomEta.n_rows;
            int subsetSize = subsetCount.size() ;
            int count, n;
            double density, prob ;

            arma::vec binomDensity(sampSize) ;
            int i = 0, j = 0;

            for( i = 0; i < sampSize ; i++) {
              density = 0;
              for( j = 0; j < subsetSize; j++) {
                prob = randomEta(i, j) ;
                prob = expit(prob) ;
                count = subsetCount(j) ;
                n = subsetN(j) ;
                if(betaDispersion & (disp < 150000) ) {
                  density += betaBinomDens(count, n ,prob, disp) ;
                } else {
                  density += R::dbinom(count, n, prob, 1) / subsetSize ;
                }
              }
              binomDensity(i) = density ;
            }


            clusterDensities.slice(subject).col(cluster) = binomDensity + importanceWeights ;

          }

          arma::vec integratedDensities ;

          integratedDensities = computeIntegratedDensities_arma(clusterDensities.slice(subject)) ;


          double priorProb = 0.5;
          if(sample >= 0) {

            clusterassignments(subset,subject) = 1 ;


            priorProb = expit(sum(isingCoefs.row(subset) * clusterassignments.col(subject)) + isingOffset) ;

          } else {
            priorProb = 0.5 ;
          }

          double pResponder;
          {
            double densityRatio;

            densityRatio = integratedDensities(0) / integratedDensities(1) * (1.0 - priorProb) / priorProb ;
            pResponder = 1.0 / (1.0 + densityRatio) ;
          }

          if(preassign(subset,subject) == 1) {
            pResponder = 1 - (1 - pResponder) * iterAssignCoef ;
          } else if(preassign(subset,subject) == 0) {
            pResponder = pResponder * iterAssignCoef ;
          }

          if(unifVec(unifPosition++,subject) < pResponder) {
            clusterassignments(subset,subject) = 1 ;
          } else {
            clusterassignments(subset,subject) = 0 ;
          }

          if(Progress::check_abort()){abort=1;}
        }
        }
        if((sample % keepEach) == 0) {

          for(int i = 0 ; i < clusterassignments.col(subject).n_elem ; i ++ ){
            assignmentMats(assignNum, i,subject) = clusterassignments(i,subject) ;
          }
          assignNum++ ;

        }
      }


      arma::vec unifVec2;

      unifVec2 = flowReMix::myrunif(nsamp_floor * nsubsets);

      arma::vec newEta(altEta.n_rows);
      newEta = nullEta.col(subject);//should copy
      arma::rowvec assignment;

      assignment = assignmentMats.slice(subject).row(assignmentMats.n_rows-1);

      std::vector<int> responderSubset;

      responderSubset = flowReMix::match2(subpopInd.col(subject),flowReMix::which2(assignment));


      // newEta[responderSubset] = altEta[responderSubset];

      flowReMix::mapByVec(newEta,altEta.col(subject),responderSubset);


      arma::vec randomEst;

      randomEst = rand.col(subject); // randomEst will  be modified?


      arma::vec MHattempts(nsubsets,arma::fill::zeros);
      arma::vec MHsuccess(nsubsets,arma::fill::zeros);

      if(sampleRandom) {
        double index_subject = index(subject);

        randomeffectMats.slice(subject) = simRandomEffectCoordinateMH_mc(Y.col(subject), N.col(subject), index_subject,
                                                   nsamp_floor, nsubsets, MHcoef,
                                                   assignment.t(), subpopInd.col(subject), newEta,
                                                   randomEst, condvar, covariance, invcov,
                                                   MHattempts, MHsuccess, unifVec2,
                                                   M, betaDispersion, keepEach, mat_size);

      }
      // std::cout<<MHsuccess;
      // std::cout<<MHattempts;
      MHsuccessrates.col(subject) = MHsuccess/MHattempts;

      prog.increment();
    }

}

    //Iterate over columns and assign to lists or just return a list of the three matrices..
      List retval = List::create(
        Named("assign") = wrap(assignmentMats),
        Named("rand") = wrap(randomeffectMats),
        Named("rate") = wrap(MHsuccessrates)
      );
      std::cout<<"\n";
      return retval;
  }
  catch (std::exception &ex){
    forward_exception_to_r(ex);
  } catch (...){
    Rf_error("C++ exception (unknown reason)");
  }
  return(List());
}
