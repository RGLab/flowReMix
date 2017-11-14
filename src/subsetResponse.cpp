#include  <RcppArmadillo.h>
using namespace Rcpp;

// headers

NumericVector expit(NumericVector x) ;

void print(NumericVector x);

// helper functions

double expit(double x) {
  return 1.0 / (1.0 + std::exp(-x)) ;
}

double unifZeroOne() {
  return rand() / (RAND_MAX + 1.) ;
}

void copyVec(NumericVector from, NumericVector to) {
  for(int i = 0 ; i < from.length() ; i++) {
    to[i] = from[i] ;
  }
}

void copyVec(IntegerVector from, IntegerVector to) {
  for(int i = 0 ; i < from.length() ; i++) {
    to[i] = from[i] ;
  }
}


NumericVector logit(NumericVector p) {
  NumericVector result = log(p / (1.0 - p)) ;
  return result ;
}

NumericMatrix computeRandomEta(NumericVector eta, NumericVector vsample) {
  int m = vsample.length() ;
  NumericMatrix result(m, eta.length()) ;

  for(int i = 0 ; i < m ; i++) {
    result(i, _) = vsample[i] + eta ;
  }

  return result ;
}

double betaBinomDens(int count, int N, double prob, double M) {
  double a = M * prob ;
  double b = M * (1 - prob) ;

  double logdens = R::lchoose(N, count) ;
  logdens += R::lbeta(count + a, N - count + b) ;
  logdens -= R::lbeta(a, b) ;
  return logdens ;
}

// [[Rcpp::export]]
NumericVector vecBetaBinomDens(NumericVector count, NumericVector N,
                               NumericVector prob, double M) {
  NumericVector logdens(count.length()) ;
  for(int i = 0 ; i < count.length() ; i++) {
    logdens[i] = betaBinomDens(count[i], N[i], prob[i], M) ;
  }

  return logdens ;
}


NumericVector computeBinomDensity(NumericVector subsetCount,
                                  NumericVector subsetN,
                                  NumericMatrix randomEta,
                                  double M,
                                  bool betaDispersion) {
  int sampSize = randomEta.nrow() ;
  int subsetSize = subsetCount.length() ;
  int count, N;
  double density, prob ;

  NumericVector binomDensity(sampSize) ;

  for(int i = 0; i < sampSize ; i++) {
    density = 0;
    for(int j = 0; j < subsetSize; j++) {
      prob = randomEta(i, j) ;
      prob = expit(prob) ;
      count = subsetCount[j] ;
      N = subsetN[j] ;
      if(betaDispersion & (M < 150000) ) {
        density += betaBinomDens(count, N ,prob, M) ;
      } else {
        density += R::dbinom(count, N, prob, 1) / subsetSize ;
      }
    }
    binomDensity[i] = density ;
  }

  return binomDensity ;
}

NumericVector computeIntegratedDensities(NumericMatrix logdens) {
  double maxdens = max(logdens) ;

  NumericVector result(2) ;
  result[0] = sum(exp(logdens(0, _) - maxdens)) ;
  result[1] = sum(exp(logdens(1, _) - maxdens)) ;

  return result ;
}

// [[Rcpp::export]]
void setNumericVectorToZero(NumericVector x) {
  int length = x.length() ;
  for(int i = 0; i < length ; i++) {
    x[i] = 0 ;
  }
}

// [[Rcpp::export]]
NumericMatrix subsetAssignGibbs(NumericVector y, NumericVector prop, NumericVector N,
                                NumericMatrix isingCoefs,
                                NumericVector nullEta, NumericVector altEta,
                                NumericMatrix covariance,
                                int nsamp, int nSubsets, int keepEach, int intSampSize,
                                NumericVector MHcoef,
                                IntegerVector popInd,
                                NumericVector unifVec, NumericVector normVec,
                                NumericVector dispersion, bool betaDispersion,
                                IntegerVector preAssignment,
                                double randomAssignProb,
                                NumericVector mprobs, double preAssignCoef,
                                double prior, bool zeroPosteriorProbs,
                                LogicalVector doNotSample,
                                NumericVector assignment) {
  NumericVector subsetNullEta, subsetAltEta, empEta, eta, etaResid ;
  NumericVector subsetProp, subsetCount, subsetN ;
  NumericVector vsample, sampNormDens, normDens, importanceWeights ;
  NumericVector binomDensity ;
  NumericMatrix randomEta ;
  NumericVector integratedDensities;
  double sigmaHat, muHat, prevMuHat ;
  double priorProb, densityRatio, pResponder ;
  int k, j, m;

  int assignNum = 0;
  NumericMatrix clusterDensities(2, intSampSize) ;
  NumericVector iterPosteriors(nSubsets) ;
  nsamp = floor(nsamp / keepEach) * keepEach ;
  NumericMatrix assignmentMatrix(int(nsamp / keepEach), nSubsets) ;
  //no need for this, just call the passed in init as assignment and reuse it.
  // NumericVector assignment(nSubsets) ;
  // for(int i = 0; i < nSubsets ; i ++) {
  //   assignment[i] = init[i] ;
  // }

  int unifPosition = 0 ;
  double isingOffset = 0 ;
  for(m = 0; m < nsamp ; m++) {
    for(j = 0; j < nSubsets ; j++) {
      if(doNotSample[j]) {
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

      subsetNullEta = nullEta[popInd == (j + 1)] ;
      subsetAltEta = altEta[popInd == (j + 1)] ;

      // Some cell populations for a specific subject may be missing
      if(subsetNullEta.length() == 0) {
        continue ;
      }

      subsetN = N[popInd == (j + 1)];
      sigmaHat = sqrt(covariance(j, j)) ;
      subsetProp = prop[popInd == (j + 1)]  ;
      empEta = logit(pmax(subsetProp, 1 / subsetN)) ;
      subsetCount = y[popInd == (j + 1)] ;

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

        sampNormDens = dnorm(vsample, muHat, sigmaHat * mcoef, TRUE) ;
        normDens = dnorm(vsample, 0, sigmaHat, TRUE) ;
        importanceWeights = normDens - sampNormDens ;
        randomEta = computeRandomEta(eta, vsample) ;
        binomDensity = computeBinomDensity(subsetCount, subsetN, randomEta,
                                           dispersion[j], betaDispersion) ;
        clusterDensities(k, _) = binomDensity + importanceWeights ;
      }

      integratedDensities = computeIntegratedDensities(clusterDensities) ;
      if(m >= 0) {
        assignment[j] = 1 ;
        // int nRespond = sum(assignment) ;
        // double multiadjust = std::log(mprobs[nRespond]) - std::log(mprobs[nRespond - 1]) ;
        priorProb = expit(sum(isingCoefs(j, _) * assignment) + isingOffset) ;
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
      for(int i = 0 ; i < assignment.length() ; i ++ ){
        assignmentMatrix(assignNum, i) = assignment[i] ;
      }
      assignNum++ ;
    }
  }


  return assignmentMatrix;
}



double computeConditionalMean(int subset,
                              NumericVector condvar,
                              NumericMatrix invcov,
                              NumericVector randomEst) {
  double conditionalMean = 0;

  for(int i = 0; i < randomEst.length() ; i++) {
    if(i != subset) {
      conditionalMean += invcov(subset, i) * randomEst[i] ;
    }
  }

  conditionalMean = conditionalMean * condvar[subset] ;
  return conditionalMean;
}

double binomDensityForMH(NumericVector count, NumericVector N,
                         NumericVector eta, double proposal,
                         double M, bool betaDispersion) {
  double prob ;
  double density = 0;

  for(int i = 0; i < eta.length(); i++) {
    prob = expit(eta[i] + proposal) ;
    if(betaDispersion & (M < 200000)) {
      density += betaBinomDens(count[i], N[i], prob, M) ;
    } else {
      density += R::dbinom(count[i], N[i], prob, true) ;
    }
  }

  return density ;
}


// [[Rcpp::export]]
NumericMatrix simRandomEffectCoordinateMH(NumericVector y, NumericVector N,
                              int i, int nsamp, int nSubsets,
                              NumericVector MHcoef,
                              IntegerVector assignment,
                              IntegerVector popInd,
                              NumericVector eta,
                              NumericVector randomEstt,
                              NumericVector condvar,
                              NumericMatrix covariance,
                              NumericMatrix invcov,
                              NumericVector MHattempts, NumericVector MHsuccess,
                              NumericVector unifVec,
                              NumericVector dispersion, bool betaDispersion,
                              int keepEach) {
  int m, j ;
  NumericVector subsetEta ;
  NumericVector subsetCount, subsetN ;
  double current, proposal, sqrtsig ;
  double condmean, newdens, olddens ;
  NumericVector randomEst = clone(randomEstt) ;

  int unifIndex = 0;

  nsamp = floor(nsamp / keepEach) * keepEach ;
  NumericMatrix sampleMatrix(int(nsamp / keepEach), nSubsets) ;
  int assignNum = 0 ;

  for(m = 0; m < nsamp ; m++) {
    for(j = 0; j < nSubsets; j++) {
      MHattempts[j] +=  1 ;
      subsetEta = eta[popInd == (j + 1)] ;
      subsetCount = y[popInd == (j + 1)] ;
      subsetN = N[popInd == (j + 1)] ;
      sqrtsig = sqrt(covariance(j, j)) ;

      // Some cell populations for a specific subject may be missing
      if(subsetEta.length() == 0) {
        continue ;
      }

      current = randomEst[j] ;
      condmean = computeConditionalMean(j, condvar, invcov, randomEst) ;
      proposal = rnorm(1)[0] * sqrtsig * MHcoef[j] + current ;
      newdens = binomDensityForMH(subsetCount, subsetN, subsetEta, proposal,
                                  dispersion[j], betaDispersion) ;
      olddens = binomDensityForMH(subsetCount, subsetN, subsetEta, current,
                                  dispersion[j], betaDispersion) ;

      newdens = newdens + R::dnorm(proposal, condmean, sqrt(condvar[j]), TRUE) ;
      olddens = olddens + R::dnorm(current, condmean, sqrt(condvar[j]), TRUE) ;

      if(unifVec[unifIndex++] < std::exp(newdens - olddens))  {
        randomEst[j] = proposal ;
        MHsuccess[j] += 1 ;
      }
    }

    if((m % keepEach) == 0) {
      for(int i = 0 ; i < randomEst.length() ; i ++ ){
        sampleMatrix(assignNum, i) = randomEst[i] ;
      }
      assignNum++ ;
    }
  }

  return sampleMatrix ;
}

// [[Rcpp::export]]
void newMHsampler(NumericMatrix assign, NumericMatrix random,
                  NumericVector initAssign, NumericVector initRand,
                  NumericVector y, NumericVector N,
                  int keepEach, double prior,
                  NumericMatrix isingCoefs,
                  IntegerVector preAssignment,
                  NumericMatrix invcov, NumericMatrix covariance,
                  NumericVector condvar,
                  NumericVector dispersion,
                  NumericVector nullEta, NumericVector altEta,
                  IntegerVector popInd,
                  NumericVector MHattempts, NumericVector MHsuccess,
                  NumericVector MHcoef) {
  NumericVector cAssign = clone(initAssign) ;
  NumericVector cRand = clone(initRand) ;
  NumericVector eta ;
  NumericVector subsetN, subsetCount ;
  double sconst;

  int nSubsets = assign.ncol() ;
  int row = 0 ;
  int iteration = 0 ;
  int zProposal ;
  int prevAssign ;
  double isingProb, isingOffset ;
  double randProposal ;
  double proposalDensity, prevDensity, MHratio ;
  double condMean ;

  while(row < assign.nrow()) {
    iteration++ ;
    for(int i = 0 ; i < nSubsets ; i++) {
      if((preAssignment[i] != -1) & (prior > 1000)) {
        cAssign[i] = preAssignment[i] ;
        continue ;
      }else if(preAssignment[i] == 0 ) {
        isingOffset = -prior ;
      } else if(preAssignment[i] == 1) {
        isingOffset = prior ;
      } else {
        isingOffset = 0 ;
      }

      // Getting relevant quantities
      subsetN = N[popInd == (i + 1)];
      subsetCount = y[popInd == (i + 1)] ;

      // Sampling assignment
      prevAssign = cAssign[i] ;
      cAssign[i] = 1 ;
      isingProb = expit(sum(isingCoefs(i, _) * cAssign) + isingOffset) ;
      cAssign[i] = prevAssign ;
      zProposal = R::rbinom(1, isingProb) ;

      // Sampling random effect
      condMean = computeConditionalMean(i, condvar, invcov, cRand) ;
      randProposal = R::rnorm(condMean, std::sqrt(condvar[i])) ;

      // Computing MH ratio
      if(zProposal < 0.5) {
        eta = nullEta[popInd == (i + 1)] ;
      } else {
        eta = altEta[popInd == (i + 1)] ;
      }
      proposalDensity = binomDensityForMH(subsetCount, subsetN, eta, randProposal, dispersion[i], true) ;

      if(cAssign[i] <  0.5) {
        eta = nullEta[popInd == (i + 1)] ;
      } else {
        eta = altEta[popInd == (i + 1)] ;
      }
      prevDensity = binomDensityForMH(subsetCount, subsetN, eta, cRand[i], dispersion[i], true) ;

      MHratio = std::exp(proposalDensity - prevDensity) ;
      // Rejecting/accepting
      if(runif(1)[0] < MHratio) {
        cAssign[i] = zProposal;
        cRand[i] = randProposal ;
      }
    }

    if(iteration % keepEach == 0) {
      //Rcpp::Rcout<<cAssign[0]<<" "<<cRand[0]<<"\n" ;
      for(int i = 0 ; i < cAssign.length() ; i ++ ){
        assign(row, i) = cAssign[i] ;
      }
      row++ ;
    }
  }

  row = 0;
  iteration = 0;
  while(row < random.nrow()) {
    for(int j = 0; j < cRand.length() ; j++) {
      int s = j;
      sconst = MHcoef[s] ;
      subsetN = N[popInd == (s + 1)];
      subsetCount = y[popInd == (s + 1)] ;
      if(cAssign[s] == 0) {
        eta = nullEta[popInd == (s + 1)] ;
      } else {
        eta = altEta[popInd == (s + 1)] ;
      }

      condMean = computeConditionalMean(s, condvar, invcov, cRand) ;
      randProposal = R::rnorm(cRand[s], sconst * std::sqrt(covariance(s, s))) ;
      proposalDensity = binomDensityForMH(subsetCount, subsetN, eta, randProposal, dispersion[s], true) ;
      proposalDensity += R::dnorm(randProposal, condMean, std::sqrt(condvar[s]), 1) ;
      prevDensity = binomDensityForMH(subsetCount, subsetN, eta, cRand[s], dispersion[s], true) ;
      prevDensity += R::dnorm(cRand[s], condMean, std::sqrt(condvar[s]), 1) ;
      MHratio = std::exp(proposalDensity - prevDensity) ;

      // Rejecting/accepting
      MHattempts[s] += 1 ;
      if(runif(1)[0] < MHratio) {
        MHsuccess[s] += 1 ;
        cRand[s] = randProposal ;
      }
    }

    if(iteration % keepEach == 0) {
      //Rcpp::Rcout<<cAssign[0]<<" "<<cRand[0]<<"\n" ;
      for(int i = 0 ; i < cAssign.length() ; i ++ ){
        random(row, i) = cRand[i] ;
      }
      row++ ;
    }
    iteration++ ;
  }
}



// void newMHsampler(NumericMatrix assign, NumericMatrix random,
//                   NumericVector initAssign, NumericVector initRand,
//                   NumericVector y, NumericVector N,
//                   int keepEach, double prior,
//                   NumericMatrix isingCoefs,
//                   IntegerVector preAssignment,
//                   NumericMatrix invcov, NumericMatrix covariance,
//                   NumericVector condvar,
//                   NumericVector dispersion,
//                   NumericVector nullEta, NumericVector altEta,
//                   IntegerVector popInd,
//                   NumericVector MHattempts, NumericVector MHsuccess,
//                   NumericVector MHcoef) {
//   NumericVector cAssign = clone(initAssign) ;
//   NumericVector cRand = clone(initRand) ;
//   NumericVector eta ;
//   NumericVector subsetN, subsetCount ;
//   double sconst;
//
//   int nSubsets = assign.ncol() ;
//   int row = 0 ;
//   int iteration = 0 ;
//   int zProposal ;
//   int prevAssign ;
//   double isingProb, isingOffset ;
//   double randProposal ;
//   double proposalDensity, prevDensity, MHratio ;
//   double condMean ;
//   int start = 0;
//
//   while(row < assign.nrow()) {
//     iteration++ ;
//     for(int i = 0 ; i < nSubsets ; i++) {
//       // Getting relevant quantities
//       subsetN = N[popInd == (i + 1)];
//       subsetCount = y[popInd == (i + 1)] ;
//
//       // Sampling assignment
//       if((preAssignment[i] != -1) & (prior > 1000)) {
//         cAssign[i] = preAssignment[i] ;
//         continue ;
//       }else if(preAssignment[i] == 0 ) {
//         isingOffset = -prior ;
//       } else if(preAssignment[i] == 1) {
//         isingOffset = prior ;
//       } else {
//         isingOffset = 0 ;
//       }
//
//       prevAssign = cAssign[i] ;
//       cAssign[i] = 1 ;
//       isingProb = expit(sum(isingCoefs(i, _) * cAssign) + isingOffset) ;
//       double zeroDens = std::log(1 - isingProb) ;
//       double oneDens = std::log(isingProb) ;
//       eta = altEta[popInd == (i + 1)] ;
//       oneDens += binomDensityForMH(subsetCount, subsetN, eta, cRand[i], dispersion[i], true) ;
//       cAssign[i] = 0;
//       eta = nullEta[popInd == (i + 1)] ;
//       zeroDens += binomDensityForMH(subsetCount, subsetN, eta, cRand[i], dispersion[i], true) ;
//       MHratio = 1 / (1 + std::exp(zeroDens - oneDens)) ;
//       if(runif(1)[0] < MHratio) {
//         cAssign[i] = 1 ;
//       } else {
//         cAssign[i] = 0 ;
//       }
//     }
//
//     // Sampling random effect
//     if(++start == cRand.length()) {
//       start = 0 ;
//     }
//     int s = start;
//     for(int j = 0; j < cRand.length() * 2; j++) {
//       if(s == cRand.length()) {
//         s = 0 ;
//       }
//       sconst = MHcoef[s] ;
//       subsetN = N[popInd == (s + 1)];
//       subsetCount = y[popInd == (s + 1)] ;
//       condMean = computeConditionalMean(s, condvar, invcov, cRand) ;
//       if(cAssign[s] == 0) {
//         eta = nullEta[popInd == (s + 1)] ;
//       } else {
//         eta = altEta[popInd == (s + 1)] ;
//       }
//
//       randProposal = R::rnorm(cRand[s], sconst * std::sqrt(covariance(s, s))) ;
//       proposalDensity = binomDensityForMH(subsetCount, subsetN, eta, randProposal, dispersion[s], true) ;
//       proposalDensity += R::dnorm(randProposal, condMean, std::sqrt(condvar[s]), 1) ;
//       prevDensity = binomDensityForMH(subsetCount, subsetN, eta, cRand[s], dispersion[s], true) ;
//       prevDensity += R::dnorm(cRand[s], condMean, std::sqrt(condvar[s]), 1) ;
//       MHratio = std::exp(proposalDensity - prevDensity) ;
//
//       // Rejecting/accepting
//       MHattempts[s] += 1 ;
//       if(runif(1)[0] < MHratio) {
//         MHsuccess[s] += 1 ;
//         cRand[s] = randProposal ;
//       }
//       s++ ;
//     }
//
//     if(iteration % keepEach == 0) {
//       //Rcpp::Rcout<<cAssign[0]<<" "<<cRand[0]<<"\n" ;
//       for(int i = 0 ; i < cAssign.length() ; i ++ ){
//         assign(row, i) = cAssign[i] ;
//         random(row, i) = cRand[i] ;
//       }
//       row++ ;
//     }
//   }
// }





