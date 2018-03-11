// Copyright (C) 2018 by Fred Hutchinson Cancer Research Center

#include <RcppArmadillo.h>
#include "./flowReMix.h"
// [[Rcpp::depends(RcppArmadillo)]]
using  namespace Rcpp;


arma::vec computeIntegratedDensities_arma(arma::mat logdens) {
  double maxdens = logdens.max();

  arma::vec result(2);
  result(0) = sum(exp(logdens.col(0) - maxdens));
  result(1) = sum(exp(logdens.col(1) - maxdens));

  return result ;
}


arma::mat computeRandomEta_arma(arma::vec eta, arma::vec vsample) {
  int m = vsample.size() ;

  arma::mat result(m, eta.size()) ;
  for(int i = 0 ; i < m ; i++) {
    result.row(i) = eta.t();
    result.row(i) = result.row(i)+vsample[i];
  }
  return result ;
}


//eta should be K x 1.
arma::mat computeRandomEta_arma(arma::mat eta, arma::vec vsample) {
  int m = vsample.size() ;
  arma::mat result(m, eta.size()) ;
  for(int i = 0 ; i < m ; i++) {
    result.row(i) = eta.t();
    result.row(i) = result.row(i)+vsample(i);
  }
  return result ;
}

inline arma::vec computeBinomDensity_arma(const arma::vec subsetCount,
                                          const arma::vec subsetN,
                                          const arma::mat randomEta,
                                          const double M,
                                          const bool betaDispersion) {
  int sampSize = randomEta.n_rows;
  int subsetSize = subsetCount.size() ;
  int count, N;
  double density, prob ;

  arma::vec binomDensity(sampSize) ;
  for(int i = 0; i < sampSize ; i++) {
    density = 0;
    for(int j = 0; j < subsetSize; j++) {
      prob = randomEta(i, j) ;
      prob = expit(prob) ;
      count = subsetCount(j) ;
      N = subsetN(j) ;
      if(betaDispersion && (M < 150000) ) {
        density += betaBinomDens(count, N ,prob, M) ;
      } else {
        density += R::dbinom(count, N, prob, 1) / subsetSize ;
      }
    }
    binomDensity(i) = density ;
  }

  return binomDensity ;
}


arma::mat subsetAssignGibbs_mc(const arma::vec y, const arma::vec prop,const  arma::vec N,
                               const arma::mat isingCoefs,
                               const arma::vec nullEta, const arma::vec altEta,
                               const arma::mat covariance,
                               const int nsamp, const int nSubsets, const int keepEach, const int intSampSize,
                               const arma::vec MHcoef,
                               const arma::vec popInd,
                               const arma::vec unifVec, const arma::vec normVec,
                               const arma::vec dispersion, const bool betaDispersion,
                               const arma::vec preAssignment,
                               const double randomAssignProb,
                               const arma::vec mprobs, const double preAssignCoef,
                               const double prior,const bool zeroPosteriorProbs,
                               const arma::vec doNotSample,
                               arma::vec assignment,const int msize) {

  arma::vec subsetNullEta, subsetAltEta, empEta, eta, etaResid ;
  arma::vec subsetProp, subsetCount, subsetN ;
  arma::vec vsample, sampNormDens, normDens, importanceWeights ;
  arma::vec binomDensity ;
  arma::mat randomEta ;
  arma::vec integratedDensities;
  double sigmaHat, muHat, prevMuHat ;
  double priorProb, densityRatio, pResponder ;
  int k, j, m;

  int assignNum = 0;
  arma::mat clusterDensities(intSampSize,2) ;
  arma::vec iterPosteriors(nSubsets) ;
  // nsamp = floor(nsamp / keepEach) * keepEach ;
  arma::mat assignmentMatrix(msize, nSubsets) ;
  //no need for this, just call the passed in init variable as assignment and reuse it.
  // NumericVector assignment(nSubsets) ;
  // for(int i = 0; i < nSubsets ; i ++) {
  //   assignment[i] = init[i] ;
  // }

  int unifPosition = 0 ;
  double isingOffset = 0 ;
  for(m = 0; m < nsamp ; m++) {
    for(j = 0; j < nSubsets ; j++) {
      if(doNotSample[j]==1.0) {
        assignment[j] = 0 ;
        continue ;
      } else if(preAssignment[j] != -1 & preAssignCoef < 10e-4 & !zeroPosteriorProbs) {
        assignment[j] = preAssignment[j] ;
        continue;
      }else if(preAssignment[j] == 0 ) {
        isingOffset = -prior ;
      } else if(preAssignment[j] == 1) {
        isingOffset = prior ;
      } else {
        isingOffset = 0 ;
      }
      arma::uvec ipop = find(popInd == (j+1));

      subsetNullEta = nullEta.elem(ipop) ;
      subsetAltEta = altEta.elem(ipop) ;

      // Some cell populations for a specific subject may be missing
      if(subsetNullEta.size() == 0) {
        continue ;
      }



      subsetN = N.elem(ipop);
      sigmaHat = sqrt(covariance(j, j)) ;
      subsetProp = prop.elem(ipop)  ;
      empEta = logit_arma(flowReMix::pmax(subsetProp, 1.0 / subsetN)) ;
      subsetCount = y.elem(ipop) ;

      // integrating densities
      double mcoef = MHcoef[j] ;
      mcoef = std::max(1.0, MHcoef[j]) ;
      for(k = 0; k < 2; k++) {
        if(k == 0) {
          eta = subsetNullEta ;
        } else {
          eta = subsetAltEta ;
        }
        etaResid = empEta - eta ;
        if(k == 0) {
          muHat = mean(etaResid) ;
          muHat = 0;
          // if(j == 0) {
          //   Rcpp::Rcout<<muHat<<" " ;
          // }
          vsample = rnorm(intSampSize, muHat, sigmaHat * mcoef) ;
        } else {
          prevMuHat = muHat ;
          muHat = mean(etaResid) ;
          muHat = 0;
          // if(j == 0) {
          //   Rcpp::Rcout<<muHat<<" " ;
          // }
          vsample = vsample - prevMuHat + muHat ;
        }




        // sampNormDens = dnorm(vsample, muHat, sigmaHat * mcoef, TRUE) ;
        sampNormDens = flowReMix::dnorm4(vsample, muHat, sigmaHat * mcoef, true) ;

        normDens = flowReMix::dnorm4(vsample, 0, sigmaHat, true) ;


        importanceWeights = normDens - sampNormDens ;


        randomEta = computeRandomEta_arma(eta, vsample) ;



        binomDensity = computeBinomDensity_arma(subsetCount, subsetN, randomEta,
                                                dispersion(j), betaDispersion) ;





        clusterDensities.col(k) = binomDensity + importanceWeights ;
      }


      integratedDensities = computeIntegratedDensities_arma(clusterDensities) ;


      if(m >= 0) {
        assignment[j] = 1 ;
        // int nRespond = sum(assignment) ;
        // double multiadjust = std::log(mprobs[nRespond]) - std::log(mprobs[nRespond - 1]) ;
        priorProb = expit(sum(isingCoefs.row(j) * assignment) + isingOffset) ;
        // if(j == 0) {
        //   Rcpp::Rcout<<priorProb<<" " ;
        // }
      } else {
        priorProb = 0.5 ;
      }


      densityRatio = integratedDensities[0] / integratedDensities[1] * (1.0 - priorProb) / priorProb ;
      pResponder = 1.0 / (1.0 + densityRatio) ;
      // if(j == 0) {
      //   Rcpp::Rcout<<integratedDensities[0]<<" "<<integratedDensities[1]<<" " ;
      //   Rcpp::Rcout<<pResponder<<"\n " ;
      // }

      if(preAssignment[j] == 1) {
        pResponder = 1 - (1 - pResponder) * preAssignCoef ;
      } else if(preAssignment[j] == 0) {
        pResponder = pResponder * preAssignCoef ;
      }


      if(unifVec[unifPosition++] < pResponder) {
        assignment[j] = 1 ;
      } else {
        assignment[j] = 0 ;
      }
    }
    if((m % keepEach) == 0) {
      for(int i = 0 ; i < assignment.size() ; i ++ ){
        assignmentMatrix(assignNum, i) = assignment[i] ;
      }
      assignNum++ ;
    }
  }

  return assignmentMatrix;
}
