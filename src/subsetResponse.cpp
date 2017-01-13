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

NumericVector computeBinomDensity(NumericVector subsetCount, NumericVector subsetN,
                                  NumericMatrix randomEta) {
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
      density += R::dbinom(count, N, prob, TRUE) / subsetSize;
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
NumericMatrix subsetAssignGibbs(NumericVector y, NumericVector prop, NumericVector N,
                                NumericMatrix isingCoefs,
                                NumericVector nullEta, NumericVector altEta,
                                NumericMatrix covariance,
                                int nsamp, int nSubsets, int keepEach, double MHcoef,
                                IntegerVector popInd) {

  NumericVector subsetNullEta, subsetAltEta, empEta, eta, etaResid ;
  NumericVector subsetProp, subsetCount, subsetN ;
  NumericVector vsample, sampNormDens, normDens, importanceWeights ;
  NumericVector binomDensity ;
  NumericMatrix randomEta ;
  NumericVector integratedDensities;
  double sigmaHat, muHat ;
  double priorProb, densityRatio, pResponder ;
  int k, j, m;

  int intSampSize = 100 ;
  int assignNum = 0;
  NumericMatrix clusterDensities(2, intSampSize) ;
  NumericVector iterPosteriors(nSubsets) ;
  NumericMatrix assignmentMatrix(int(floor(nsamp / keepEach)), nSubsets) ;
  NumericVector assignment(nSubsets) ;

  for(m = 0; m < nsamp ; m++) {
    for(j = 0; j < nSubsets ; j++) {
      subsetNullEta = nullEta[popInd == (j + 1)] ;
      subsetAltEta = altEta[popInd == (j + 1)] ;

      sigmaHat = sqrt(covariance(j, j)) ;
      subsetProp = prop[popInd == (j + 1)]  ;
      empEta = logit(subsetProp + 0.00001) ;
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
        vsample = rnorm(intSampSize, muHat, sigmaHat * MHcoef) ;
        sampNormDens = dnorm(vsample, muHat, sigmaHat * MHcoef, TRUE) ;
        normDens = dnorm(vsample, muHat, sigmaHat, TRUE) ;
        importanceWeights = normDens - sampNormDens ;
        randomEta = computeRandomEta(eta, vsample) ;
        binomDensity = computeBinomDensity(subsetCount, subsetN, randomEta) ;
        clusterDensities(k, _) = binomDensity + importanceWeights ;
      }

      integratedDensities = computeIntegratedDensities(clusterDensities) ;
      if(m >= 1) {
        priorProb = expit(sum(isingCoefs(j, _) * assignment)) ;
      } else {
        priorProb = 0.5 ;
      }

      densityRatio = integratedDensities[0] / integratedDensities[1] * (1.0 - priorProb) / priorProb ;
      pResponder = 1.0 / (1.0 + densityRatio) ;
      if(unifZeroOne() < pResponder) {
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



