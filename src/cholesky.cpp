#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
void cholupL_Rcpp(arma::mat& L,
                  arma::vec& x){
  
  // It would be best to use this function with a wrapper, as it modifies BOTH
  // L and x in-place.
  
  int p = x.n_elem;
  int pp1 = p + 1;
  double r;
  double c;
  double s;
  double y;
  
  for(auto [j, a, b] = std::tuple{0, L.begin(), x.begin()}; j < p; j++, a += pp1, ++b){
    
    r = sqrt(*a * *a + *b * *b);
    c = *a / r;
    s = *b / r;
    
    for(auto [u, v] = std::tuple{b, a}; u != x.end(); ++u, ++v){
      
      // The variables iterate over...
      // u: iterates over x
      // v: iterates over L's column
      y = *v;
      *v = c * *v + s * *u;
      *u = s * y - c * *u;
      
    }
    
  }
  
}

// [[Rcpp::export]]
void cholupU_Rcpp(arma::mat& U,
                  arma::vec& x){
  
  // It would be best to use this function with a wrapper, as it modifies BOTH
  // U and x in-place.
  
  int p = x.n_elem;
  int pp1 = p + 1;
  double r;
  double c;
  double s;
  double y;
  
  for(auto [j, a, b] = std::tuple{0, U.begin(), x.begin()}; j < p; j++, a += pp1, ++b){
    
    r = sqrt(*a * *a + *b * *b);
    c = *a / r;
    s = *b / r;
    
    for(auto [u, v] = std::tuple{b, a}; u != x.end(); ++u, v += p){
      
      // The variables iterate over...
      // u: iterates over x
      // v: iterates over L's column
      y = *v;
      *v = c * *v + s * *u;
      *u = s * y - c * *u;
      
    }
    
  }
  
}

// [[Rcpp::export]]
void choldownL_Rcpp(arma::mat& L,
                    arma::vec& x){
  
  // It would be best to use this function with a wrapper, as it modifies BOTH
  // L and x in-place.
  
  int p = x.n_elem;
  int pp1 = p + 1;
  double r;
  double c;
  double s;
  double y;
  
  for(auto [j, a, b] = std::tuple{0, L.begin(), x.begin()}; j < p; j++, a += pp1, ++b){
    
    r = sqrt(*a * *a - *b * *b);
    c = *a / r;
    s = *b / r;
    
    for(auto [u, v] = std::tuple{b, a}; u != x.end(); ++u, ++v){
      
      // The variables iterate over...
      // v: iterates over x
      // u: iterates over L's column
      y = *v;
      *v = c * *v - s * *u;
      *u = c * *u - s * y;
      
    }
    
  }
  
}

// [[Rcpp::export]]
void choldownU_Rcpp(arma::mat& U,
                    arma::vec& x){
  
  // It would be best to use this function with a wrapper, as it modifies BOTH
  // U and x in-place.
  
  int p = x.n_elem;
  int pp1 = p + 1;
  double r;
  double c;
  double s;
  double y;
  
  for(auto [j, a, b] = std::tuple{0, U.begin(), x.begin()}; j < p; j++, a += pp1, ++b){
    
    r = sqrt(*a * *a - *b * *b);
    c = *a / r;
    s = *b / r;
    
    for(auto [u, v] = std::tuple{b, a}; u != x.end(); ++u, v += p){
      
      // The variables iterate over...
      // u: iterates over x
      // v: iterates over L's column
      y = *v;
      *v = c * *v - s * *u;
      *u = c * *u - s * y;
      
    }
    
  }
  
}
