#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List comvar2_cpp(NumericVector x, NumericVector y, NumericVector nxny, int B){
  // Compare the variances of two independent groups.
  double estx = var(x);
  double esty = var(y);
  double diff = estx - esty;
  NumericVector v(2);
  // v(1) = x.length();
  // v(2) = y.length();
  int nmin = min(nxny);
  
// NumericVector bootdiff(B);  // Preallocate storage for bootstrap results
// 
// // Perform bootstrap 
// for(int i = 0; i < B; i++) {
//   
//   // Bootstrap: sample with replacement
//   NumericVector bootx = sample(x, nmin, true);
//   NumericVector booty = sample(y, nmin, true);
//   bootdiff(i) = var(bootx) - var(booty);
//   
// }
// 
// bootdiff = bootdiff.sort();
// int ilow = 15;
// int ihi = 584;
// if(nmin < 250) {
//   ilow = 13; 
//   ihi = 586; 
// }
// if(nmin < 180) {
//   ilow = 10;
//   ihi = 589;
// }
// if(nmin < 80) {
//   ilow = 7;
//   ihi = 592;
// }
// if(nmin < 40) {
//   ilow = 6;
//   ihi = 593;
// }
// 
// int cilow = round((ilow/599)*B) + 1;
// int cihi = round((ihi/599)*B);
// NumericVector ci(2);
// ci(1) = bootdiff(cilow);
// ci(2) = bootdiff(cihi);

// Return results
List res;
res["diff"] = diff;
res["est_x"] = estx;
res["est_y"] = esty;
// res["ci"] = ci;
res["v"] = v;
res["nmin"] = nmin;
return res;
}