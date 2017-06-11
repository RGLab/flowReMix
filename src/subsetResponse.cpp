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
  double logdens = std::lgamma(N  + 1) ;
  logdens += std::lgamma(count + a) ;
  logdens += std::lgamma(N - count + b) ;
  logdens += std::lgamma(a + b) ;
  logdens -= std::lgamma(count + 1) ;
  logdens -= std::lgamma(N - count + 1) ;
  logdens -= std:: lgamma(N + a + b) ;
  logdens -= std::lgamma(a) ;
  logdens -= std::lgamma(b) ;
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
      if(betaDispersion) {
        density += betaBinomDens(count, N ,prob, M) ;
      } else {
        density += R::dbinom(count, N, prob, TRUE) / subsetSize ;
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
  for(int i; i < length ; i++) {
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
                                NumericVector mprobs) {
  NumericVector subsetNullEta, subsetAltEta, empEta, eta, etaResid ;
  NumericVector subsetProp, subsetCount, subsetN ;
  NumericVector vsample, sampNormDens, normDens, importanceWeights ;
  NumericVector binomDensity ;
  NumericMatrix randomEta ;
  NumericVector integratedDensities;
  double sigmaHat, muHat ;
  double priorProb, densityRatio, pResponder ;
  int k, j, m;

  int assignNum = 0;
  NumericMatrix clusterDensities(2, intSampSize) ;
  NumericVector iterPosteriors(nSubsets) ;
  nsamp = floor(nsamp / keepEach) * keepEach ;
  NumericMatrix assignmentMatrix(int(nsamp / keepEach), nSubsets) ;
  NumericVector assignment(nSubsets) ;

  int unifPosition = 0 ;
  for(m = 0; m < nsamp ; m++) {
    for(j = 0; j < nSubsets ; j++) {
      if(preAssignment[j] != -1) {
        assignment[j] = preAssignment[j] ;
        continue ;
      }

      subsetNullEta = nullEta[popInd == (j + 1)] ;
      subsetAltEta = altEta[popInd == (j + 1)] ;

      // Some cell populations for a specific subject may be missing
      if(subsetNullEta.length() == 0) {
        continue ;
      }

      sigmaHat = sqrt(covariance(j, j)) ;
      subsetProp = prop[popInd == (j + 1)]  ;
      empEta = logit(subsetProp + 10e-5) ;
      subsetCount = y[popInd == (j + 1)] ;
      subsetN = N[popInd == (j + 1)];

      // integrating densities
      for(k = 0; k < 2; k++) {
        if(k == 0) {
          eta = subsetNullEta ;
        } else {
          eta = subsetAltEta ;
        }
        etaResid = empEta - eta ;
        muHat = mean(etaResid) ;
        vsample = rnorm(intSampSize, muHat, sigmaHat * MHcoef[j]) ;
        sampNormDens = dnorm(vsample, muHat, sigmaHat * MHcoef[j], TRUE) ;
        normDens = dnorm(vsample, 0, sigmaHat, TRUE) ;
        importanceWeights = normDens - sampNormDens ;
        randomEta = computeRandomEta(eta, vsample) ;
        binomDensity = computeBinomDensity(subsetCount, subsetN, randomEta,
                                           dispersion[j], betaDispersion) ;
        clusterDensities(k, _) = binomDensity + importanceWeights ;
      }

      integratedDensities = computeIntegratedDensities(clusterDensities) ;
      if(m >= 1) {
        assignment[j] = 1 ;
        int nRespond = sum(assignment) ;
        double multiadjust = std::log(mprobs[nRespond]) - std::log(mprobs[nRespond - 1]) ;
        priorProb = expit(sum(isingCoefs(j, _) * assignment) + multiadjust) ;
      } else {
        priorProb = 0.5 ;
      }

      densityRatio = integratedDensities[0] / integratedDensities[1] * (1.0 - priorProb) / priorProb ;
      pResponder = 1.0 / (1.0 + densityRatio) ;
      pResponder = std::max(pResponder, randomAssignProb) ;
      pResponder = std::min(pResponder, 1 - randomAssignProb) ;
      if(unifVec[unifPosition++] < pResponder) {
        assignment[j] = 1 ;
      } else {
        assignment[j] = 0 ;
      }
    }

    if((m % keepEach) == 0) {
      assignmentMatrix(assignNum++, _) = assignment ;
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
    if(betaDispersion) {
      density += betaBinomDens(count[i], N[i], prob, M) ;
    } else {
      density += R::dbinom(count[i], N[i], prob, TRUE) ;
    }
  }

  return density ;
}


// [[Rcpp::export]]
NumericMatrix randomEffectCoordinateMH(NumericVector y, NumericVector N,
                              int i, int nsamp, int nSubsets,
                              NumericVector MHcoef,
                              IntegerVector assignment,
                              IntegerVector popInd,
                              NumericVector eta,
                              NumericVector randomEst,
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
      sampleMatrix(assignNum++, _) = randomEst ;
    }
  }

  return sampleMatrix ;
}

