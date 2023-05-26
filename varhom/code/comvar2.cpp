#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ci_cpp(NumericVector x,int d,int B){
  return round(x/599*B,d);
  }

// [[Rcpp::export]]
// List comvar2_cpp(NumericVector x, NumericVector y, int n, int B){
NumericVector comvar2_cpp(NumericVector x, NumericVector y, int n, int B){
  // Compare the variances of two independent groups.
  // double estx = var(x);
  // double esty = var(y);
  // double diff = estx - esty;
  // NumericVector v(2);
  // v(1) = x.length();
  // v(2) = y.length();
  // int nmin = min(nxny);
  
  NumericVector bootdiff(B);  // Preallocate storage for bootstrap results
  
  // Perform bootstrap 
  for(int i = 0; i < B; i++) {
    
    // Bootstrap: sample with replacement
    NumericVector bootx = sample(x, n, true);
    NumericVector booty = sample(y, n, true);
    bootdiff(i) = var(bootx) - var(booty);
    
  }

  // int ilow = 15;
  // int ihi = 584;
  // if(n < 250) {
  //   ilow = 13; 
  //   ihi = 586; 
  // }
  // if(n < 180) {
  //   ilow = 10;
  //   ihi = 589;
  // }
  // if(n < 80) {
  //   ilow = 7;
  //   ihi = 592;
  // }
  // if(n < 40) {
  //   ilow = 6;
  //   ihi = 593;
  // }

  bootdiff = bootdiff.sort();
  return bootdiff;
  // NumericVector ci(2);
  // NumericVector ci_lo = ci_cpp(ilow,0,B);
  // ci(0) = bootdiff[ci_lo[0]];
  // NumericVector ci_hi = ci_cpp(ihi,0,B) -1;
  // ci(1) = bootdiff[ci_hi[0]];
  
  // Return results
  // List res;
  // res["diff"] = diff;
  // res["est_x"] = estx;
  // res["est_y"] = esty;
  // res["ilo"] = ilow;
  // res["ihi"] = ihi;
  // res["lo"] = ci_lo;
  // res["hi"] = ci_hi; 
  // res["ci"] = ci;
  // res["B"] = B;
  // return res;
}

// Rcpp::cppFunction("NumericVector mr(NumericVector x,int d) {return round(x/599*1000,d);}")

