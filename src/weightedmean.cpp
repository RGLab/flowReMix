#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
// [[Rcpp::plugins(cpp17)]]
//[[Rcpp::export]]
NumericVector weightedMean(NumericVector x, NumericVector weights, bool na_rm)
{

  if (weights.length() != x.length())
    throw std::domain_error("'x' and 'w' must have the same length");
    if(na_rm){
      LogicalVector nas = !is_na(x);
      weights = weights[nas];
      x = x[nas];
    }
    std::vector<double> w(weights.begin(),weights.end());
    std::vector<double> X(x.begin(),x.end());
    auto weighedmean = [] (std::vector<double> &X,std::vector<double> &w){
      double sum = 0;
      double sumw = 0;
      for(int i=0;i<X.size();++i){
        double wi =  w.at(i);
        double xi = X.at(i);
        sum += xi*wi;
        sumw += wi;
      }
      return(sum/sumw);
    };
    double res;
    res=weighedmean(X,w);
    return(wrap(res));
}

/*** R
# library(microbenchmark)
# microbenchmark(weighted.mean(x = x,w = weights,na.rm = FALSE),
# weightedMean(x = x,w = weights,na_rm = FALSE),100)
weights = runif(10)
x = rnorm(10)
all.equal(weightedMean(x = x,w = weights,na_rm = FALSE), weighted.mean(x,weights,na.rm=FALSE))
*/
