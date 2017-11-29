#include <math.h>
#include <boost/math/special_functions/beta.hpp>
#include <tgmath.h>
#define ISINT(x) !(fabs((x) - std::nearbyint(x)) > 1e-7)
// double  lfastchoose(double n, double k)
// {
//   double tmp  = log(boost::math::beta(n - k + 1., k + 1.));
//   return -log(n + 1.) - tmp;
// }
/* mathematically the same:
less stable typically, but useful if n-k+1 < 0 : */
static
  double lfastchoose2(double n, double k)
  {
    double r;
    r = lgamma(n - k + 1.);
    return std::lgamma(n + 1.) - std::lgamma(k + 1.) - r;
  }

double lbeta_cpp(double a, double b)
{
  double r;
  r = std::lgamma(a+b);
  return std::lgamma(a) + std::lgamma(b) - r;
}

double lchoose_cpp(double n, double k)
{
  double k0 = k;
  k = nearbyint(k);
  /* NaNs propagated correctly */
static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required");
  if(isnan(n) || isnan(k)) return n + k;
  if (fabs(k - k0) > 1e-7)
    // std::cout<<"'k' "<<k0<<" must be integer, rounded to "<<k;
  if (k < 2) {
    if (k <	 0) return -std::numeric_limits<double>::infinity();
    if (k == 0) return 0.;
    /* else: k == 1 */
    return log(fabs(n));
  }
  /* else: k >= 2 */
  if (n < 0) {
    return lchoose_cpp(-n+ k-1, k);
  }
  else if (ISINT(n)) {
    n = std::nearbyint(n);
    if(n < k) return -std::numeric_limits<double>::infinity();
    /* k <= n :*/
    if(n - k < 2) return lchoose_cpp(n, n-k); /* <- Symmetry */
    /* else: n >= k+2 */
    return lfastchoose2(n, k);
  }
  /* else non-integer n >= 0 : */
  if (n < k-1) {
    return lfastchoose2(n, k);
  }
  return lfastchoose2(n, k);
 }


