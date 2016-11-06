#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

NumericVector expit(NumericVector x) {
  NumericVector r;
  r = 1.0 / (1.0 + exp(-x));
  return(r);
}

void print(NumericVector x){
  for(int j=0;j<x.length();++j){
    Rcpp::Rcout<<x(j)<<" ";
  }
  Rcpp::Rcout<<std::endl;
}

NumericVector dbinom_vec(NumericVector y, NumericVector N, NumericVector mu){
  NumericVector v(y.length());
  for(int i=0;i<mu.length();++i){
    v(i) = Rf_dbinom(y(i),N(i),mu(i),true);
  }
  return(v);
}

// [[Rcpp::export]]
void MH(const List randomSampList,  NumericMatrix lastMean,  NumericMatrix estimatedRandomEffects, const NumericVector y,const  NumericVector N, const NumericMatrix randomEffectSamp, const NumericVector eta, const int i, IntegerVector popInd,  arma::mat invcov,  NumericVector accept,double iter,double rate) {
  int k = 0;
  NumericVector diff;
  int j = 0;
  NumericVector randomSamp;
  NumericVector mu;
  NumericVector currentlogLik, newlogLik;
  for(k = 0; k < 2; ++k){
    randomSamp = randomSampList[k];
    NumericVector currentSamp = lastMean(2*i - 2 + k,_);
    NumericVector randEst = estimatedRandomEffects(2*i - 2 + k,_);
    mu = currentSamp[popInd-1];
    mu = expit(eta+mu);
    NumericVector devr = currentSamp - randEst;
    arma::vec dev = arma::vec(devr.begin(),devr.length(),false);
    currentlogLik = sum(dbinom_vec(y, N, mu)) - 0.5 * dev.t() * invcov * dev;
    for(j=0;j<randomEffectSamp.ncol();++j){
      NumericVector newSamp = randomEffectSamp(_,j);
      newSamp = newSamp + currentSamp;
      mu = newSamp[popInd-1];
      mu = expit(eta + mu);
      devr = newSamp - randEst;
      dev = arma::vec(devr.begin(),devr.length(),false);
      newlogLik = sum(dbinom_vec(y, N, mu)) - 0.5 * dev.t() * invcov * dev;
      diff = exp(newlogLik - currentlogLik);
      if((runif(1))(0) < diff(0)) {
        currentSamp = newSamp;
        currentlogLik = newlogLik;
        accept = accept + 1;
      }
    }
    lastMean(2*i - 2 + k,_) = currentSamp;
    NumericVector currentEst = estimatedRandomEffects(2*i - 2 + k,_);
    estimatedRandomEffects(2*i - 2 + k,_) = currentEst + (currentSamp - currentEst) / pow((iter + 1.0),rate);
  }
}


