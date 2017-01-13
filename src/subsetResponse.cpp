#include  <RcppArmadillo.h>
using namespace Rcpp;

NumericVector expit(NumericVector x) ;

void print(NumericVector x);

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

NumericVector computeResponderProb(NumericMatrix coefficients, NumericVector assignment, int subsetIndex, int nSubsets) {
  NumericVector eta = coefficients(subsetIndex, 0) ;
  for(int i = 0; i < nSubsets - 1; i++) {
    if(i != subsetIndex) {
      eta[0] += coefficients(subsetIndex, i + 1) * assignment[i] ;
    }
  }

  eta = expit(eta) ;
  return eta;
}

NumericVector computeBinomDensity(NumericVector subsetCount, NumericVector subsetN,
                                  NumericMatrix randomEta) {
  int sampSize = randomEta.nrow() ;
  int subsetSize = subsetCount.length() ;
  double density ;
  NumericVector prob ;

  NumericVector binomDensity(sampSize) ;

  for(int i = 0; i < sampSize ; i++) {
    density = 0;
    for(int j = 0; j < subsetSize; j++) {
      // print(randomEta(i, _)) ;
      prob = randomEta(i, j) ;
      prob = expit(prob) ;
      double probability = prob[0] ;
      // Rcpp::Rcout<<"inner "<<i<<" "<< j <<"\n";
      int count = subsetCount[j] ;
      int N = subsetN[j] ;
      density += R::dbinom(count, N, probability, TRUE) ;
      // Rcpp::Rcout<<"inner "<<i<<" "<< j <<"\n";
    }
    binomDensity[i] = density ;
  }

  return binomDensity ;
}

NumericVector computeIntegratedDensities(NumericMatrix logdens) {
  double maxdens = max(logdens) ;
  //Rcpp::Rcout<< maxdens<<" ";
  NumericVector result(2) ;

  for(int i = 0; i < 2 ; i++) {
    for(int j = 0; j < logdens.ncol() ; j++) {
      result(i) += exp(logdens(i, j) - maxdens) ;
      // Rcpp::Rcout<< logdens(i, j)<<" ";
    }
  }

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
  NumericVector pResponder, densityRatio, integratedDensities;
  double sigmaHat, muHat;
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
        // print(importanceWeights) ;
        binomDensity = computeBinomDensity(subsetCount, subsetN, randomEta) ;
        // print(binomDensity) ;
        // Rcpp::Rcout<<"M J= "<<m << " " << j << "\n";
        clusterDensities(k, _) = binomDensity + importanceWeights ;
      }

      if(m >= 1) {
        pResponder = computeResponderProb(isingCoefs, assignment, j, nSubsets) ;
      } else {
        pResponder = 0.5 ;
      }

      // Rcpp::Rcout<<j << "\n";
      integratedDensities = computeIntegratedDensities(clusterDensities) ;
      // print(integratedDensities) ;
      densityRatio = integratedDensities[0] / integratedDensities[1] * (1.0 - pResponder) / pResponder ;
      // print(densityRatio) ;
      pResponder = 1.0 / (1.0 + densityRatio) ;
      // print(densityRatio) ;
      assignment[j] = R::rbinom(1, pResponder[0]) ;
    }

    if((m % keepEach) == 0) {
      assignmentMatrix(assignNum++, _) = assignment ;
      Rcpp::Rcout<<assignment[0]<<assignment[1]<<"\n";
    }
  }

  return assignmentMatrix;
}



