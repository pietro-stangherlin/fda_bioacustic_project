// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


// explicit kroenecker here
// [[Rcpp::export]]
arma::mat ComputeBvar(const arma::mat& A,
                      const arma::mat& B,
                      const arma::mat& C) {
  // Compute Kronecker product of B and C
  arma::mat K = arma::kron(B, C);
  
  // Matrix multiplication A %*% (B ⊗ C)
  arma::mat result = A * K * A.t();
  
  return result;
}
