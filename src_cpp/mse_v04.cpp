#include <Rcpp.h>
#include <RcppParallel.h>
#include <cmath>

using namespace Rcpp;
using namespace RcppParallel;

// [[Rcpp::export]]
double calculateDTH(NumericVector temperature,
                    NumericVector daylength,
                    double G, double Th, double Lc,
                    double A, double B, double DVIstar) {
  
  double DVI = 0.0;
  
  const int n = temperature.size();
  const double invG = 1.0 / G;
  
  for (int i = 0; i < n; ++i) {
    const double T = temperature[i];
    const double L = daylength[i];
    
    // Temperature response (0..1)
    const double fT = 1.0 / (1.0 + std::exp(-A * (T - Th)));
    
    // Pre-DVI* (juvenile): temperature-driven only
    const double DVR1 = invG * fT;
    
    // Post-DVI*: photoperiod effect when L < Lc, else 0
    // NOTE: (1 - exp(B*(L-Lc))) is >0 when L < Lc and B>0
    double DVR2 = invG * (1.0 - std::exp(B * (L - Lc))) * fT;
    
    double DVR = 0.0;
    if (DVI < DVIstar) {
      DVR = DVR1;
    } else {
      if (L < Lc) DVR = DVR2;
      else       DVR = 0.0;
    }
    
    // Safety: never allow negative accumulation
    if (DVR < 0.0) DVR = 0.0;
    
    DVI += DVR;
    
    if (DVI >= 1.0) {
      // Return day count (1..n)
      return static_cast<double>(i + 1);
    }
  }
  
  // If never reaches heading, return max window length
  return static_cast<double>(n);
}


// [[Rcpp::export]]
double calculateDTH_optimized(const NumericVector& temperature,
                              const NumericVector& daylength,
                              double G, double Th, double Lc,
                              double A, double B, double DVIstar) {
  
  double DVI = 0.0;
  double DTH = 0.0;
  
  const double invG = 1.0 / G;
  
  auto t_it = temperature.begin();
  auto l_it = daylength.begin();
  
  while (t_it != temperature.end() && l_it != daylength.end()) {
    const double T = *t_it;
    const double L = *l_it;
    
    const double fT  = 1.0 / (1.0 + std::exp(-A * (T - Th)));
    const double DVR1 = invG * fT;
    double DVR2 = invG * (1.0 - std::exp(B * (L - Lc))) * fT;
    
    double DVR = 0.0;
    
    if (DVI < DVIstar) {
      DVR = DVR1;
    } else {  
      if (L < Lc) DVR = DVR2;
      else       DVR = 0.0;
    }
    
    if (DVR < 0.0) DVR = 0.0;
    
    const double prev_DVI = DVI;
    DVI += DVR;
    
    // Count days as we step forward (each loop = 1 day)
    DTH += 1.0;
    
    if (DVI >= 1.0) {
      // Optional: fractional day interpolation within the last day
      // Only valid if DVR > 0
      if (DVR > 0.0) {
        const double need = 1.0 - prev_DVI;   // remaining DVI to reach 1
        const double frac = need / DVR;       // fraction of the day used
        // We already added 1 full day above, so subtract unused portion
        // If frac in (0,1], then unused = (1 - frac)
        return DTH - (1.0 - frac);
      } else {
        // Shouldn't happen if DVI increased, but safe fallback
        return DTH;
      }
    }
    
    ++t_it;
    ++l_it;
  }
  
  return DTH;  // If never exceeds 1, equals number of simulated days
}


// [[Rcpp::export]]
double MSE_cpp(const NumericVector& dth,
               const NumericMatrix& temperatures,
               const NumericMatrix& daylengths,
               double G, double Th, double Lc,
               double A, double B, double DVIstar) {
  
  const size_t number_of_sites = temperatures.ncol();
  NumericVector model(number_of_sites);
  
#pragma omp parallel for
  for (size_t i = 0; i < number_of_sites; ++i) {
    model[i] = calculateDTH_optimized(temperatures(_, i),
                                      daylengths(_, i),
                                      G, Th, Lc, A, B, DVIstar);
  }
  
  double mse = 0.0;
  for (size_t i = 0; i < number_of_sites; ++i) {
    const double err = model[i] - dth[i];
    mse += err * err;
  }
  mse /= static_cast<double>(number_of_sites);
  
  return mse;
}
