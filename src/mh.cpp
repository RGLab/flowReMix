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
    v(i) = Rf_dbinom(y(i),N(i),mu(i),1);
  }
  return(v);
}

// [[Rcpp::export]]
void MH(NumericMatrix lastMean,  NumericMatrix estimatedRandomEffects,
        const NumericVector y,const  NumericVector N,
        const NumericMatrix randomEffectSamp,
        const int i, IntegerVector popInd,  arma::mat invcov,
        NumericVector accept,double iter,double rate, NumericVector unifs,
        const NumericVector nullEta, const NumericVector altEta) {
  int k = 0;
  NumericVector diff;
  int j = 0;
  int I = 0;
  I = i-1;
  NumericVector eta;
  NumericVector mu;
  NumericVector currentlogLik, newlogLik;
  for(k = 0; k < 2; k++){
    if(k == 0) {
      eta = nullEta ;
    } else {
      eta = altEta ;
    }
    NumericVector currentSamp = lastMean(2*I + k,_);
    mu = currentSamp[popInd-1];
    mu = expit(eta+mu);

    arma::vec dev = arma::vec(currentSamp.begin(),currentSamp.length(),true);
    currentlogLik = sum(dbinom_vec(y, N, mu)) - 0.5 * dev.t() * invcov * dev;
    for(j=0;j<randomEffectSamp.ncol();j++){
      NumericVector newSamp = randomEffectSamp(_,j);
      newSamp = newSamp + currentSamp;
      mu = newSamp[popInd-1];
      mu = expit(eta + mu);
      dev = arma::vec(newSamp.begin(),newSamp.length(),true);
      newlogLik = sum(dbinom_vec(y, N, mu)) - 0.5 * dev.t() * invcov * dev;
      diff = exp(newlogLik - currentlogLik);
      if(unifs(j) < diff(0)) {
        currentSamp = newSamp;
        currentlogLik = newlogLik;
        accept(0) = accept(0)+1;
      }
      lastMean(2*I + k,_) = currentSamp;
      NumericVector currentEst = estimatedRandomEffects(2*I + k,_);
      estimatedRandomEffects(2*I + k,_) = currentEst + ((currentSamp - currentEst) / pow(std::max(iter - 5.0, 1.0), rate));
    }
    // Rprintf("Accepted %d on subject %d\n",accept,I);
   }
}

// [[Rcpp::export]]
NumericVector zero(NumericVector accept){
  accept(0)=0;
  return(accept);
}


