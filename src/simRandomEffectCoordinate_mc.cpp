// Copyright (C) 2018 by Fred Hutchinson Cancer Research Center
#include <RcppArmadillo.h>
#include "./flowReMix.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

double binomDensityForMH_arma(arma::vec count, arma::vec N,
                              arma::vec eta, double proposal,
                              double M, bool betaDispersion) {
  double prob;
  double density = 0;

  for (int i = 0; i < eta.size(); i++) {
    prob = expit(eta(i) + proposal);
    if (betaDispersion & (M < 200000)) {
      density += betaBinomDens(count(i), N(i), prob, M);
    } else {
      density += R::dbinom(count(i), N(i), prob, true);
    }
  }

  return density;
}

double computeConditionalMean_arma(int subset,
                                   arma::vec condvar,
                                   arma::mat invcov,
                                   arma::vec randomEst) {
  double conditionalMean = 0;

  for (int i = 0; i < randomEst.size() ; i++) {
    if (i != subset) {
      conditionalMean += invcov(subset, i) * randomEst[i];
    }
  }
  conditionalMean = conditionalMean * condvar[subset];
  return conditionalMean;
}



void simRandomEffectCoordinateMH_mc(arma::mat& sampleMatrix,
                                         const arma::vec y,
                                         const arma::vec N,
                                         const int i,
                                         int nsamp,
                                         const int nSubsets,
                                         const arma::vec MHcoef,
                                         const arma::vec assignment,
                                         const arma::vec popInd,
                                         const arma::vec eta,
                                         const arma::vec randomEst_in,
                                         const arma::vec condvar,
                                         const arma::mat covariance,
                                         const arma::mat invcov,
                                         arma::vec& MHattempts,
                                         arma::vec& MHsuccess,
                                         const arma::vec dispersion,
                                         const bool betaDispersion,
                                         const int keepEach,
                                         const int msize) {
  int m, j;
  arma::vec subsetEta;
  arma::vec subsetCount, subsetN;
  double current, proposal, sqrtsig;
  double condmean, newdens, olddens;
  arma::vec randomEst = randomEst_in;  // should copy

  int unifIndex = 0;

  // int nsamp = floor(nsamp / keepEach) * keepEach ;
  // arma::mat sampleMatrix(msize, nSubsets);
  int assignNum = 0;


  for (m = 0; m < nsamp ; m++) {
    for (j = 0; j < nSubsets; j++) {
      arma::uvec ipop = find(popInd == (j+1));
      if (ipop.size() == 0) {
        continue;
      }
      MHattempts(j) +=  1;
      subsetEta = eta.elem(ipop);
      // subsetEta = eta[popInd == (j + 1)] ;
      subsetCount = y.elem(ipop);
      // subsetCount = y[popInd == (j + 1)] ;
      // subsetN = N[popInd == (j + 1)] ;
      subsetN = N.elem(ipop);
      sqrtsig = sqrt(covariance(j, j));
      // Some cell populations for a specific subject may be missing


      current = randomEst(j);
      condmean = computeConditionalMean_arma(j, condvar, invcov, randomEst);
      proposal = flowReMix::myrnorm(1)(0) * sqrtsig * MHcoef(j) + current;
      newdens = binomDensityForMH_arma(subsetCount, subsetN,
                                       subsetEta, proposal,
                                       dispersion(j), betaDispersion);
      olddens = binomDensityForMH_arma(subsetCount, subsetN,
                                       subsetEta, current,
                                       dispersion(j), betaDispersion);

      newdens = newdens + R::dnorm4(proposal, condmean, sqrt(condvar(j)), true);
      olddens = olddens + R::dnorm4(current, condmean, sqrt(condvar(j)), true);

      // if(unifVec(unifIndex++) < std::exp(newdens - olddens))  {
      if (R::runif(0, 1) < std::exp(newdens - olddens)) {
        randomEst(j) = proposal;
        MHsuccess(j) += 1;
      }
    }

    if ((m % keepEach) == 0) {
      for (int i = 0 ; i < randomEst.size(); i++) {
        sampleMatrix(assignNum, i) = randomEst(i);
      }
      assignNum++;
    }
  }

  return;
}
