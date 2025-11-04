#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
#include <string>
#include <cmath>
#include <limits> // for std::numeric_limits

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::export]]
int calculateDTH(NumericVector temperature, NumericVector daylength, double G, double Th, double Lc, double A, double B, double DVIstar) {
  double DVI = 0.0;
  double DVR=0.0;
  int DTH = 0;
  // DVIstar=0; Not used for now and reset to 0 (fixed)
  
  for (int i = 0; i < temperature.length(); ++i) {
    double T = temperature(i);
    double L = daylength(i);
    // double DVR1 = (1 / G) * (1 / (1 + exp(-A * (T - Th)))); not used when DVIstar=0
    double DVR1=0; //dummy
    double DVR2 = (1 / G) * (1 - exp(B * (L - Lc))) / (1 + exp(-A * (T - Th)));
    if (DVI<DVIstar) {DVR=DVR1;}
    else if (DVI>=DVIstar && L<Lc) {DVR=DVR2;}
    else {DVR=0;} // DVI>=DVIstar && L>=Lc
    
    DVI += DVR;
    if (DVI > 1) {
      DTH = i;
      break;
    }
  }
  return DTH;
}

// [[Rcpp::export]]
int calculateDTH_optimized(const NumericVector& temperature, const NumericVector& daylength, 
                           double G, double Th, double Lc, double A, double B, double DVIstar) {
  double DVI = 0.0;
  double DVR = 0.0;
  int DTH = 0;
  
  // Precompute constant
  double invG = 1.0 / G;
  
  // Use iterators for efficiency
  auto temp_it = temperature.begin();
  auto day_it = daylength.begin();
  
  while (temp_it != temperature.end() && day_it != daylength.end()) {
    double T = *temp_it;
    double L = *day_it;
    
    // Precompute common terms to avoid redundant calculations
    double exp_A = 1.0 / (1.0 + exp(-A * (T - Th)));
    double DVR2 = invG * (1 - exp(B * (L - Lc))) * exp_A;
    
    // Simplify conditionals,  MUST BE fixed when DVIstar is used
    if (DVI >= DVIstar && L < Lc) {
      DVR = DVR2;
    } else {
      DVR = 0.0; // Covers both DVI < DVIstar and DVI >= DVIstar && L >= Lc
    }
    
    DVI += DVR;
    
    if (DVI > 1.0) {
      return DTH;
    }
    
    ++temp_it;
    ++day_it;
    ++DTH;
  }
  
  return DTH; // If DVI never exceeds 1, return the full length
}


// [[Rcpp::export]]
double MSE_cpp(const Rcpp::NumericVector& dth,
               const Rcpp::NumericMatrix& temperatures,
               const Rcpp::NumericMatrix& daylengths,
               double G, double Th, double Lc,
               double A, double B, double DVIstar) {
  
  size_t number_of_sites = temperatures.ncol();
  Rcpp::NumericVector model(number_of_sites);
  
  // Parallelized loop using OpenMP
#pragma omp parallel for
  for (size_t i = 0; i < number_of_sites; ++i) {
    model[i] = calculateDTH_optimized(temperatures(Rcpp::_, i), 
                                daylengths(Rcpp::_, i), 
                                G, Th, Lc, A, B, DVIstar);
  }
  
  // Compute mean squared error
  double mse = 0.0;
  for (size_t i = 0; i < number_of_sites; ++i) {
    mse += std::pow(model[i] - dth[i], 2);
  }
  mse /= number_of_sites;
  
  return mse;
}
