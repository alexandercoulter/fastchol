#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
void cholRcpp_Rcpp(arma::mat& R,
                   const arma::mat& M,
                   const bool& lower = true){
  
  // This function is a wrapper for Rcpp library Cholesky decomposition, which
  // permits specification of the lower/upper Cholesky factor. Calculating the
  // lower triangular factor tends to be twice as fast as calculating the upper,
  // so it is the default mode.
  
  // This function modifies R in-place, so it is advised to use this function
  // with its own wrapper.
  
  if(lower){
    
    chol(R, M, "lower");
    
  } else {
    
    chol(R, M, "upper");
    
  }
  
}

// [[Rcpp::export]]
void LtoU_Rcpp(arma::mat& U,
               const arma::mat& L){
  
  // This function copies the lower-triangular part of L (including main
  // diagonal) into the upper-triangular part of U. It modifies U in-place, so
  // best to use with a wrapper function which sets the dimension of U.
  
  unsigned int p = L.n_cols;
  unsigned int pp1 = p + 1;
  
  for(auto [j, a, b] = std::tuple{0u, L.begin(), U.begin()}; j < p; j++, a += pp1, b += pp1){
    
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
  
  unsigned int p = L.n_cols;
  unsigned int pp1 = p + 1;
  
  for(auto [j, a, b] = std::tuple{0u, L.begin(), U.begin()}; j < p; j++, a += pp1, b += pp1){
    
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
  
  auto l = L.begin();
  auto xptr = x.memptr();
  unsigned int p = x.n_elem;
  unsigned int j, k;
  double r, c, s, y, z;
  
  for(j = 0; j < p; j++){
    
    l += j;
    y = *l;
    z = xptr[j];
    
    r = sqrt(y * y + z * z);
    c = y / r;
    s = z / r;
    
    // Set first dimension
    *l = r;
    ++l;
    
    // Calculate for remaining dimensions
    for(k = j + 1; k < p; k++, ++l){
      
      // The variables iterate over...
      // z: iterates over x
      // y: iterates over L's column
      y = *l;
      z = xptr[k];
      
      *l = c * y + s * z;
      xptr[k] = s * y - c * z;
      
    }
    
  }
  
}

// [[Rcpp::export]]
void cholupKL_Rcpp(arma::mat& L,
                   arma::mat& X){
  
  // It would be best to use this function with a wrapper, as it modifies BOTH
  // L and x in-place.
  
  unsigned int p = X.n_rows;
  unsigned int K = X.n_cols;
  unsigned int pp1 = p + 1;
  double r, c, s, y;

  for(auto [j, a, b] = std::tuple{0u, L.begin(), X.begin()}; j < p; j++, a += pp1, ++b){
    
    for(auto [k, xk] = std::tuple{0u, b}; k < K; k++, xk += p){
      
      r = sqrt(*a * *a + *xk * *xk);
      c = *a / r;
      s = *xk / r;
      
      // Set first dimension
      *a = r;
      
      // Calculate for remaining dimensions
      for(auto [u, v] = std::tuple{xk + 1, a + 1}; u != X.end_col(k); ++u, ++v){
        
        y = *v;
        *v = c * y + s * *u;
        *u = s * y - c * *u;
        
      }
      
    }
    
  }
  
}

// [[Rcpp::export]]
void cholupU_Rcpp(arma::mat& U,
                  arma::vec& x){
  
  // It would be best to use this function with a wrapper, as it modifies BOTH
  // U and x in-place.
  
  unsigned int p = x.n_elem;
  unsigned int pp1 = p + 1;
  double r, c, s, y;

  for(auto [j, a, b] = std::tuple{0u, U.begin(), x.begin()}; j < p; j++, a += pp1, ++b){
    
    r = sqrt(*a * *a + *b * *b);
    c = *a / r;
    s = *b / r;
    
    // Set first dimension
    *a = r;
    
    // Calculate for second dimension
    for(auto [u, v] = std::tuple{b + 1, a + p}; u != x.end(); ++u, v += p){
      
      // The variables iterate over...
      // u: iterates over x
      // v: iterates over U's row
      y = *v;
      *v = c * y + s * *u;
      *u = s * y - c * *u;
      
    }
    
  }
  
}

// [[Rcpp::export]]
void cholupU2_Rcpp(arma::mat& U,
                   arma::vec& x){
  
  // It would be best to use this function with a wrapper, as it modifies BOTH
  // U and x in-place.
  
  auto u = U.begin();
  auto xptr = x.memptr();
  unsigned int p = x.n_elem;
  unsigned int pp1 = p + 1;
  double r, c, s, y, z;
  
  for(auto [j, a, b] = std::tuple{0u, U.begin(), x.begin()}; j < p; j++, a += pp1, ++b){
    
    r = sqrt(*a * *a + *b * *b);
    c = *a / r;
    s = *b / r;
    
    // Set first dimension
    *a = r;
    
    // Calculate for second dimension
    for(auto [u, v] = std::tuple{b + 1, a + p}; u != x.end(); ++u, v += p){
      
      // The variables iterate over...
      // u: iterates over x
      // v: iterates over U's row
      y = *v;
      *v = c * y + s * *u;
      *u = s * y - c * *u;
      
    }
    
  }
  
}

// [[Rcpp::export]]
void cholupKU_Rcpp(arma::mat& U,
                   arma::mat& X){
  
  // It would be best to use this function with a wrapper, as it modifies BOTH
  // U and x in-place.
  
  unsigned int p = X.n_rows;
  unsigned int K = X.n_cols;
  unsigned int pp1 = p + 1;
  double r, c, s, y;
  
  for(auto [j, a, b] = std::tuple{0u, U.begin(), X.begin()}; j < p; j++, a += pp1, ++b){
    
    for(auto [k, xk] = std::tuple{0u, b}; k < K; k++, xk += p){
      
      r = sqrt(*a * *a + *xk * *xk);
      c = *a / r;
      s = *xk / r;
      
      // Set first dimension
      *a = r;
      
      // Calculate for second dimension
      for(auto [u, v] = std::tuple{xk + 1, a + p}; u != X.end_col(k); ++u, v += p){
        
        // The variables iterate over...
        // u: iterates over X's column
        // v: iterates over U's row
        y = *v;
        *v = c * y + s * *u;
        *u = s * y - c * *u;
        
      }
      
    }
    
  }
  
}

// [[Rcpp::export]]
void choldownL_Rcpp(arma::mat& L,
                    arma::vec& x){
  
  // It would be best to use this function with a wrapper, as it modifies BOTH
  // L and x in-place.
  
  unsigned int p = x.n_elem;
  unsigned int pp1 = p + 1;
  double r, c, s, y;
  
  for(auto [j, a, b] = std::tuple{0u, L.begin(), x.begin()}; j < p; j++, a += pp1, ++b){
    
    r = sqrt(*a * *a - *b * *b);
    c = *a / r;
    s = *b / r;
    
    for(auto [u, v] = std::tuple{b, a}; u != x.end(); ++u, ++v){
      
      // The variables iterate over...
      // v: iterates over x
      // u: iterates over L's column
      y = *v;
      *v = c * y - s * *u;
      *u = c * *u - s * y;
      
    }
    
  }
  
}

// [[Rcpp::export]]
void choldownU_Rcpp(arma::mat& U,
                    arma::vec& x){
  
  // It would be best to use this function with a wrapper, as it modifies BOTH
  // U and x in-place.
  
  unsigned int p = x.n_elem;
  unsigned int pp1 = p + 1;
  double r, c, s, y;
  
  for(auto [j, a, b] = std::tuple{0u, U.begin(), x.begin()}; j < p; j++, a += pp1, ++b){
    
    r = sqrt(*a * *a - *b * *b);
    c = *a / r;
    s = *b / r;
    
    for(auto [u, v] = std::tuple{b, a}; u != x.end(); ++u, v += p){
      
      // The variables iterate over...
      // u: iterates over x
      // v: iterates over L's column
      y = *v;
      *v = c * y - s * *u;
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
  
  unsigned int km1 = k - 1;
  unsigned int p = L.n_cols;
  unsigned int pp1 = p + 1;
  unsigned int pm1 = p - 1;
  
  if(km1 == pm1){
    
    // If k is the last row/column, we copy up until that:
    for(unsigned int j = 0; j < pm1; j++){
      
      for(auto [u, v] = std::tuple{L.begin_col(j) + j, L0.begin_col(j) + j}; v != L0.end_col(j); ++u, ++v){
        
        *v = *u;
        
      }
      
    }
    
  } else if(km1 == 0){
    
    // Perform rank-1 update on lower-right corner
    unsigned int pk = p - k;
    arma::colvec x(pk);
    for(auto [i, j] = std::tuple{L.begin_col(km1) + k, x.begin()}; j != x.end(); ++i, ++j){
      
      *j = *i;
      
    }
    double r, c, s;
    unsigned int jm1;
    
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
    for(unsigned int j = 0; j < km1; j++){
      
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
    unsigned int pk = p - k;
    arma::colvec x(pk);
    for(auto [i, j] = std::tuple{L.begin_col(km1) + k, x.begin()}; j != x.end(); ++i, ++j){
      
      *j = *i;
      
    }
    double r, c, s;
    unsigned int jm1;
    
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
  
  unsigned int km1 = k - 1;
  unsigned int p = U.n_cols;
  unsigned int pp1 = p + 1;
  unsigned int pm1 = p - 1;
  unsigned int jm1;
  
  if(km1 == pm1){
    
    // If k is the last row/column, we copy up until that:
    for(unsigned int j = 0; j < pm1; j++){
      
      for(auto [u, v, w] = std::tuple{U.begin_col(j), U0.begin_col(j), 0u}; w <= j; ++u, ++v, w++){
        
        *v = *u;
        
      }
      
    }
    
  } else if(km1 == 0){
    
    // Perform rank-1 update on lower-right corner
    unsigned int pk = p - k;
    arma::colvec x(pk);
    for(auto [i, j] = std::tuple{U.begin_col(k) + km1, x.begin()}; j != x.end(); i += p, ++j){
      
      *j = *i;
      
    }
    double r, c, s;
    
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
    unsigned int pk = p - k;
    for(unsigned int j = 0; j < km1; j++){
      
      // Copy the top-left (above k) corner
      for(auto [u, v, w] = std::tuple{U.begin_col(j), U0.begin_col(j), 0u}; w <= j; ++u, ++v, w++){
        
        *v = *u;
        
      }
      
      // Copy the top-right columns after kth column
      for(auto [u, v, w] = std::tuple{U.begin_col(k) + j, U0.begin_col(km1) + j, 0u}; w < pk; u += p, v += pm1, w++){
        
        *v = *u;
        
      }
      
    }
    
    // Perform rank-1 update on lower-right corner
    arma::colvec x(pk);
    for(auto [i, j] = std::tuple{U.begin_col(k) + km1, x.begin()}; j != x.end(); i += p, ++j){
      
      *j = *i;
      
    }
    double r, c, s;
    
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
  
  unsigned int pin = Lin.n_cols;
  unsigned int pinp1 = pin + 1;
  unsigned int pout = z.n_elem;
  unsigned int poutp1 = pout + 1;
  unsigned int km1 = k - 1;
  
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
    double r, c, s;
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
    for(auto [j, a, k0, k1, kk] = std::tuple{0u, x.begin(), Lin.begin(), Lout.begin(), Lout.begin() + km1}; a != x.end(); j++, ++a, k0 += pinp1, k1 += poutp1, kk += pout){
      
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
    double r, c, s;
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
    for(auto [j, a, k0, k1, kk] = std::tuple{0u, x.begin(), Lin.begin(), Lout.begin(), Lout.begin() + km1}; a != x.end(); j++, ++a, k0 += pinp1, k1 += poutp1, kk += pout){
      
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

// [[Rcpp::export]]
void choladdU_Rcpp(arma::mat& Uout,
                   const arma::mat& Uin,
                   const arma::vec& z,
                   const int& k){
  
  unsigned int pin = Uin.n_cols;
  unsigned int pinp1 = pin + 1;
  unsigned int pout = z.n_elem;
  unsigned int poutp1 = pout + 1;
  unsigned int km1 = k - 1;
  
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
    Uout(km1, km1) = omega;
    
    // Copy U12, and calculate/fix y <-- x * U12
    for(auto [b, k1] = std::tuple{y.begin(), Uout.begin_col(k) + km1}; b != y.end(); ++b, k1 += pout){
      
      *b /= omega;
      *k1 = *b;
      
    }
    
    // choldown of lower-right corner
    double r, c, s;
    for(auto [b, k0, k1] = std::tuple{y.begin(), Uin.begin_col(km1) + km1, Uout.begin_col(k) + k}; b != y.end(); ++b, k0 += pinp1, k1 += poutp1){
      
      // j iterates over the number of entries; 
      r = sqrt(*k0 * *k0 - *b * *b);
      c = *k0 / r;
      s = *b / r;
      
      for(auto [u, v0, v1] = std::tuple{b, k0, k1}; u != y.end(); ++u, v0 += pin, v1 += pout){
        
        // The variables iterate over...
        // u: iterates over y
        // v0: iterates over Uin's row
        // v1: iterates over Uout's row
        *v1 = c * *v0 - s * *u;
        *u = c * *u - s * *v0;
        
      }
      
    }
    
  } else if(km1 < pin){
    
    arma::vec x = z.head(km1);
    arma::vec y = z.tail(pout - k);
    
    // Copy U11 and calculate x <-- U11^{-1}x
    for(auto [j, a, k0, k1, kk] = std::tuple{0u, x.begin(), Uin.begin(), Uout.begin(), Uout.begin_col(km1)}; a != x.end(); j++, ++a, k0 += pinp1, k1 += poutp1, ++kk){
      
      // Fix current entry of x:
      *a /= *k0;
      *k1 = *k0;
      *kk = *a;
      
      // Loop through remaining entries in the row
      if(j < km1){
        
        for(auto [u, v0, v1] = std::tuple{a + 1, k0 + pin, k1 + pout}; u != x.end(); ++u, v0 += pin, v1 += pout){
          
          // u iterates through x
          // v0 iterates through Uin's column
          // v1 iterates through Uout's column
          *u -= *v0 * *a;
          *v1 = *v0;
          
        }
        
      }
      
    }
    
    // Fix the kth diagonal entry
    double omega = sqrt(w - arma::dot(x, x));
    Uout(km1, km1) = omega;
    
    // Copy U12, and calculate/fix y <-- x * U12
    
    // Outer loop will be over y (columns of Uin, Uout)
    for(auto [b, j] = std::tuple{y.begin(), km1}; b != y.end(); ++b, j++){
      
      // Inner loop over x (rows of Uin, Uout)
      for(auto [a, k0, k1] = std::tuple{x.begin(), Uin.begin_col(j), Uout.begin_col(j + 1)}; a != x.end(); ++a, ++k0, ++k1){
        
        *b -= *k0 * *a;
        *k1 = *k0;
        
      }
      
    }
    for(auto [b, k1] = std::tuple{y.begin(), Uout.begin_col(k) + km1}; b != y.end(); ++b, k1 += pout){
      
      *b /= omega;
      *k1 = *b;
      
    }
    
    double r, c, s;
    for(auto [b, k0, k1] = std::tuple{y.begin(), Uin.begin_col(km1) + km1, Uout.begin_col(k) + k}; b != y.end(); ++b, k0 += pinp1, k1 += poutp1){
      
      // j iterates over the number of entries; 
      r = sqrt(*k0 * *k0 - *b * *b);
      c = *k0 / r;
      s = *b / r;
      
      for(auto [u, v0, v1] = std::tuple{b, k0, k1}; u != y.end(); ++u, v0 += pin, v1 += pout){
        
        // The variables iterate over...
        // u: iterates over y
        // v0: iterates over Uin's row
        // v1: iterates over Uout's row
        *v1 = c * *v0 - s * *u;
        *u = c * *u - s * *v0;
        
      }
      
    }
    
  } else {
    
    arma::vec x = z.head(km1);
    
    // Copy U11 and calculate x <-- U11^{-1}x
    for(auto [j, a, k0, k1, kk] = std::tuple{0u, x.begin(), Uin.begin(), Uout.begin(), Uout.begin_col(km1)}; a != x.end(); j++, ++a, k0 += pinp1, k1 += poutp1, ++kk){
      
      // Fix current entry of x:
      *a /= *k0;
      *k1 = *k0;
      *kk = *a;
      
      // Loop through remaining entries in the row
      if(j < km1){
        
        for(auto [u, v0, v1] = std::tuple{a + 1, k0 + pin, k1 + pout}; u != x.end(); ++u, v0 += pin, v1 += pout){
          
          // u iterates through x
          // v0 iterates through Uin's column
          // v1 iterates through Uout's column
          *u -= *v0 * *a;
          *v1 = *v0;
          
        }
        
      }
      
    }
    
    // Fix the kth diagonal entry
    double omega = sqrt(w - arma::dot(x, x));
    Uout(km1, km1) = omega;
    
  }
  
}

// [[Rcpp::export]]
void Lsolvex_Rcpp(const arma::mat& L,
                  arma::vec& x){
  
  // This function will modify x in-place, so it is best to use this with a
  // function wrapper.
  
  auto a = x.memptr();
  auto l = L.begin();
  unsigned int p = x.n_elem;
  unsigned int j, k;
  double b;
  
  for(j = 0; j < p; j++){
    
    l += j;
    // Fix current entry of x:
    a[j] /= *l;
    b = a[j];
    ++l;
    
    for(k = j + 1; k < p; k++, ++l){
      
      a[k] -= b * *l;
      
    }
    
  }
  
}

// [[Rcpp::export]]
void LsolveX_Rcpp(const arma::mat& L,
                  arma::mat& X){
  
  // This function will modify x in-place, so it is best to use this with a
  // function wrapper.
  auto x = X.memptr();
  auto l = L.begin();
  unsigned int n = X.n_rows;
  unsigned int p = X.n_cols;
  unsigned int kn = -n;
  unsigned int i, j, k;
  double b;
  
  for(k = 0; k < p; k++){
    
    l = L.begin();
    kn += n;
    
    for(i = 0; i < n; i++){
      
      l += i;
      // Fix current entry of x:
      x[kn + i] /= *l;
      b = x[kn + i];
      ++l;
      
      for(j = i + 1; j < n; j++, ++l){
        
        x[kn + j] -= b * *l;
        
      }
      
    }
    
  }
  
}

// [[Rcpp::export]]
void Usolvex_Rcpp(const arma::mat& U,
                  arma::vec& x){
  
  // This function will modify x in-place, so it is best to use this with a
  // function wrapper.
  
  auto xptr = x.memptr();
  auto u = U.end();
  --u;
  unsigned int n = x.n_elem;
  unsigned int nm1 = n - 1;
  unsigned int ni;
  unsigned int i;
  int j;
  double a;
  
  for(i = 0; i < n; i++){
    
    u -= i;
    ni = nm1 - i;
    
    // Fix current entry of x
    xptr[ni] /= *u;
    a = xptr[ni];
    --u;
    
    for(j = ni - 1; j >= 0; j--, --u){
      
      xptr[j] -= a * *u;
      
    }
    
  }
  
}

// [[Rcpp::export]]
void UsolveX_Rcpp(const arma::mat& U,
                  arma::mat& X){
  
  // This function will modify x in-place, so it is best to use this with a
  // function wrapper.
  
  auto xptr = X.memptr();
  auto u = U.end();
  unsigned int n = X.n_rows;
  unsigned int nm1 = n - 1;
  unsigned int p = X.n_cols;
  unsigned int ni;
  unsigned int i, k;
  unsigned int maxk = n * p;
  int j, minj;
  double a;
  
  for(k = nm1; k < maxk; k += n){
    
    u = U.end();
    --u;
    minj = k - nm1;
    
    for(i = 0; i < n; i++){
      
      u -= i;
      ni = k - i;
      
      // Fix current entry of x
      xptr[ni] /= *u;
      a = xptr[ni];
      --u;
      
      for(j = ni - 1; j >= minj; j--, --u){
        
        xptr[j] -= a * *u;
        
      }
      
    }
    
  }
  
}
