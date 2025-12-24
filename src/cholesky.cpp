#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
void LtoU_Rcpp(arma::mat& U,
               const arma::mat& L){
  
  // This function copies the lower-triangular part of L (including main
  // diagonal) into the upper-triangular part of U. It modifies U in-place, so
  // best to use with a wrapper function which sets the dimension of U.
  
  int p = L.n_cols;
  int pp1 = p + 1;
  
  for(auto [j, a, b] = std::tuple{0, L.begin(), U.begin()}; j < p; j++, a += pp1, b += pp1){
    
    for(auto [ell, u] = std::tuple{a, b}; ell != L.end_col(j); ++ell, u += p){
      
      *u = *ell;
      
    }
    
  }
  
}

// [[Rcpp::export]]
void UtoL_Rcpp(arma::mat& L,
               const arma::mat& U){
  
  // This function copies the upper-triangular part of U (including main
  // diagonal) into the lower-triangular part of L. It modifies L in-place, so
  // best to use with a wrapper function which sets the dimension of L.
  
  int p = L.n_cols;
  int pp1 = p + 1;
  
  for(auto [j, a, b] = std::tuple{0, L.begin(), U.begin()}; j < p; j++, a += pp1, b += pp1){
    
    for(auto [ell, u] = std::tuple{a, b}; ell != L.end_col(j); ++ell, u += p){
      
      *ell = *u;
      
    }
    
  }
  
}

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

// [[Rcpp::export]]
void choldropL_Rcpp(arma::mat& L0,
                    const arma::mat& L,
                    const int& k){
  
  // The matrix L0 should have size = (nrow(L) - 1) * (ncol(L) - 1); it is
  // modified in-place.
  
  int km1 = k - 1;
  int p = L.n_cols;
  int pp1 = p + 1;
  int pm1 = p - 1;
  
  if(km1 == pm1){
    
    // If k is the last row/column, we copy up until that:
    for(int j = 0; j < L0.n_cols; j++){
      
      for(auto [u, v] = std::tuple{L.begin_col(j) + j, L0.begin_col(j) + j}; v != L0.end_col(j); ++u, ++v){
        
        *v = *u;
        
      }
      
    }
    
  } else if(km1 == 0){
    
    // Perform rank-1 update on lower-right corner
    int pk = p - k;
    arma::colvec x(pk);
    for(auto [i, j] = std::tuple{L.begin_col(km1) + k, x.begin()}; j != x.end(); ++i, ++j){
      
      *j = *i;
      
    }
    double r;
    double c;
    double s;
    int jm1;
    
    for(auto [j, a, b] = std::tuple{k, L.begin_col(k) + k, x.begin()}; j < p; j++, a += pp1, ++b){
      
      r = sqrt(*a * *a + *b * *b);
      c = *a / r;
      s = *b / r;
      jm1 = j - 1;
      
      for(auto [u, v, w] = std::tuple{b, a, L0.begin_col(jm1) + jm1}; u != x.end(); ++u, ++v, ++w){
        
        // The variables iterate over...
        // u: iterates over x
        // v: iterates over L's column
        // w: iterates over L0's column
        *w = c * *v + s * *u;
        *u = c * *u - s * *v;
        
      }
      
    }
    
  } else {
    
    // Copy the left side...
    for(int j = 0; j < km1; j++){
      
      // Copy the top-left (above k) corner
      for(auto [u, v] = std::tuple{L.begin_col(j) + j, L0.begin_col(j) + j}; u != L.begin_col(j) + km1; ++u, ++v){
        
        *v = *u;
        
      }
      
      // Copy the bottom-left rows after kth row
      for(auto [u, v] = std::tuple{L.begin_col(j) + k, L0.begin_col(j) + km1}; u != L.end_col(j); ++u, ++v){
        
        *v = *u;
        
      }
      
    }
    
    // Perform rank-1 update on lower-right corner
    int pk = p - k;
    arma::colvec x(pk);
    for(auto [i, j] = std::tuple{L.begin_col(km1) + k, x.begin()}; j != x.end(); ++i, ++j){
      
      *j = *i;
      
    }
    double r;
    double c;
    double s;
    int jm1;
    
    for(auto [j, a, b] = std::tuple{k, L.begin_col(k) + k, x.begin()}; j < p; j++, a += pp1, ++b){
      
      r = sqrt(*a * *a + *b * *b);
      c = *a / r;
      s = *b / r;
      jm1 = j - 1;
      
      for(auto [u, v, w] = std::tuple{b, a, L0.begin_col(jm1) + jm1}; u != x.end(); ++u, ++v, ++w){
        
        // The variables iterate over...
        // u: iterates over x
        // v: iterates over L's column
        // w: iterates over L0's column
        *w = c * *v + s * *u;
        *u = c * *u - s * *v;
        
      }
      
    }
    
  }
  
}

// [[Rcpp::export]]
void choldropU_Rcpp(arma::mat& U0,
                    const arma::mat& U,
                    const int& k){
  
  // The matrix U0 should have size = (nrow(U) - 1) * (ncol(U) - 1); it is
  // modified in-place.
  
  int km1 = k - 1;
  int p = U.n_cols;
  int pm1 = p - 1;
  int pp1 = p + 1;
  int jm1;
  
  if(km1 == pm1){
    
    // If k is the last row/column, we copy up until that:
    for(int j = 0; j < U0.n_cols; j++){
      
      for(auto [u, v, w] = std::tuple{U.begin_col(j), U0.begin_col(j), 0}; w <= j; ++u, ++v, w++){
        
        *v = *u;
        
      }
      
    }
    
  } else if(km1 == 0){
    
    // Perform rank-1 update on lower-right corner
    int pk = p - k;
    arma::colvec x(pk);
    for(auto [i, j] = std::tuple{U.begin_col(k) + km1, x.begin()}; j != x.end(); i += p, ++j){
      
      *j = *i;
      
    }
    double r;
    double c;
    double s;
    
    for(auto [j, a, b] = std::tuple{k, U.begin(), x.begin()}; j < p; j++, a += pp1, ++b){
      
      jm1 = j - 1;
      r = sqrt(*a * *a + *b * *b);
      c = *a / r;
      s = *b / r;
      
      for(auto [u, v, w] = std::tuple{b, a, U0.begin_col(jm1) + jm1}; u != x.end(); ++u, v += p, w += p){
        
        // The variables iterate over...
        // u: iterates over x
        // v: iterates over U's row
        // w: iterates over U0's row
        *w = c * *v + s * *u;
        *u = c * *u - s * *v;
        
      }
      
    }
    
  } else {
    
    // Copy the top side...
    int pk = p - k;
    for(int j = 0; j < km1; j++){
      
      // Copy the top-left (above k) corner
      for(auto [u, v, w] = std::tuple{U.begin_col(j), U0.begin_col(j), 0}; w <= j; ++u, ++v, w++){
        
        *v = *u;
        
      }
      
      // Copy the top-right columns after kth column
      for(auto [u, v, w] = std::tuple{U.begin_col(k) + j, U0.begin_col(km1) + j, 0}; w < pk; u += p, v += pm1, w++){
        
        *v = *u;
        
      }
      
    }
    
    // Perform rank-1 update on lower-right corner
    arma::colvec x(pk);
    for(auto [i, j] = std::tuple{U.begin_col(k) + km1, x.begin()}; j != x.end(); i += p, ++j){
      
      *j = *i;
      
    }
    double r;
    double c;
    double s;
    
    for(auto [j, a, b] = std::tuple{k, U.begin_col(k) + k, x.begin()}; j < p; j++, a += pp1, ++b){
      
      jm1 = j - 1;
      r = sqrt(*a * *a + *b * *b);
      c = *a / r;
      s = *b / r;
      
      for(auto [u, v, w] = std::tuple{b, a, U0.begin_col(jm1) + jm1}; u != x.end(); ++u, v += p, w += pm1){
        
        // The variables iterate over...
        // u: iterates over x
        // v: iterates over U's row
        // w: iterates over U0's row
        *w = c * *v + s * *u;
        *u = c * *u - s * *v;
        
      }
      
    }
    
  }
  
}

// [[Rcpp::export]]
void choladdL_Rcpp(arma::mat& Lout,
                   const arma::mat& Lin,
                   const arma::vec& z,
                   const int& k){
  
  int pin = Lin.n_cols;
  int pinp1 = pin + 1;
  int pout = z.n_elem;
  int poutp1 = pout + 1;
  int km1 = k - 1;
  
  // There are three cases, corresponding to which new row/column k will be:
  // 1. km1 = 0
  // 2. 0 < km1 < p
  // 3. km1 == p
  
  // No matter the case, we set w = z(km1)
  double w = z(km1);
  
  // Case 1.
  if(km1 == 0){
    
    arma::vec y = z.tail(pout - k);
    
    // Fix the kth diagonal entry
    double omega = sqrt(w);
    Lout(km1, km1) = omega;
    
    // Copy L21, and calculate/fix y <-- L21 * x
    for(auto [b, k1] = std::tuple{y.begin(), Lout.begin_col(km1) + k}; b != y.end(); ++b, ++k1){
      
      *b /= omega;
      *k1 = *b;
      
    }
    
    // choldown of lower-right corner
    double r;
    double c;
    double s;
    for(auto [b, k0, k1] = std::tuple{y.begin(), Lin.begin_col(km1) + km1, Lout.begin_col(k) + k}; b != y.end(); ++b, k0 += pinp1, k1 += poutp1){
      
      // j iterates over the number of entries;
      r = sqrt(*k0 * *k0 - *b * *b);
      c = *k0 / r;
      s = *b / r;
      
      for(auto [u, v0, v1] = std::tuple{b, k0, k1}; u != y.end(); ++u, ++v0, ++v1){
        
        // The variables iterate over...
        // u: iterates over y
        // v0: iterates over Lin's column
        // v1: iterates over Lout's column
        *v1 = c * *v0 - s * *u;
        *u = c * *u - s * *v0;
        
      }
      
    }
    
  } else if(km1 < pin){
    
    arma::vec x = z.head(km1);
    arma::vec y = z.tail(pout - k);
    
    // Copy L11 and calculate x <-- L11^{-1}x
    for(auto [j, a, k0, k1, kk] = std::tuple{0, x.begin(), Lin.begin(), Lout.begin(), Lout.begin() + km1}; a != x.end(); j++, ++a, k0 += pinp1, k1 += poutp1, kk += pout){
      
      // Fix current entry of x:
      *a /= *k0;
      *k1 = *k0;
      *kk = *a;
      
      // Loop through remaining entries in the column
      if(j < km1){
        
        for(auto [u, v0, v1] = std::tuple{a + 1, k0 + 1, k1 + 1}; u != x.end(); ++u, ++v0, ++v1){
          
          // u iterates through x
          // v0 iterates through Lin's column
          // v1 iterates through Lout's column
          *u -= *v0 * *a;
          *v1 = *v0;
          
        }
        
      }
      
    }
    
    // Fix the kth diagonal entry
    double omega = sqrt(w - arma::dot(x, x));
    Lout(km1, km1) = omega;
    
    // Copy L21, and calculate/fix y <-- L21 * x
    
    // Outer loop will be over y (rows of Lin, Lout)
    for(auto [b, j] = std::tuple{y.begin(), km1}; b != y.end(); ++b, j++){
      
      // Inner loop over x (columns of Lin, Lout)
      for(auto [a, k0, k1] = std::tuple{x.begin(), Lin.begin() + j, Lout.begin() + j + 1}; a != x.end(); ++a, k0 += pin, k1 += pout){
        
        *b -= *k0 * *a;
        *k1 = *k0;
        
      }
      
    }
    for(auto [b, k1] = std::tuple{y.begin(), Lout.begin_col(km1) + k}; b != y.end(); ++b, ++k1){
      
      *b /= omega;
      *k1 = *b;
      
    }
    
    // choldown of lower-right corner
    double r;
    double c;
    double s;
    for(auto [b, k0, k1] = std::tuple{y.begin(), Lin.begin_col(km1) + km1, Lout.begin_col(k) + k}; b != y.end(); ++b, k0 += pinp1, k1 += poutp1){
      
      // j iterates over the number of entries; 
      r = sqrt(*k0 * *k0 - *b * *b);
      c = *k0 / r;
      s = *b / r;
      
      for(auto [u, v0, v1] = std::tuple{b, k0, k1}; u != y.end(); ++u, ++v0, ++v1){
        
        // The variables iterate over...
        // u: iterates over y
        // v0: iterates over Lin's column
        // v1: iterates over Lout's column
        *v1 = c * *v0 - s * *u;
        *u = c * *u - s * *v0;
        
      }
      
    }
    
  } else {
    
    arma::vec x = z.head(km1);
    
    // Copy L11 and calculate x <-- L11^{-1}x
    for(auto [j, a, k0, k1, kk] = std::tuple{0, x.begin(), Lin.begin(), Lout.begin(), Lout.begin() + km1}; a != x.end(); j++, ++a, k0 += pinp1, k1 += poutp1, kk += pout){
      
      // Fix current entry of x:
      *a /= *k0;
      *k1 = *k0;
      *kk = *a;
      
      // Loop through remaining entries in the column
      if(j < km1){
        
        for(auto [u, v0, v1] = std::tuple{a + 1, k0 + 1, k1 + 1}; u != x.end(); ++u, ++v0, ++v1){
          
          // u iterates through x
          // v0 iterates through Lin's column
          // v1 iterates through Lout's column
          *u -= *v0 * *a;
          *v1 = *v0;
          
        }
        
      }
      
    }
    
    // Fix the kth diagonal entry
    double omega = sqrt(w - arma::dot(x, x));
    Lout(km1, km1) = omega;
    
  }
  
}