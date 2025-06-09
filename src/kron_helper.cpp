// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include<Rcpp>
using namespace Rcpp;
using namespace RcppParallel;


// [[Rcpp::export]]
NumericMatrix KroneckerProdByBlocks(const NumericMatrix& A,
                                    const NumericMatrix& B,
                                    const NumericMatrix& C) {
  int nrowA = A.nrow();
  int ncolA = A.ncol();
  int ncolB = B.ncol();
  int nrowB = B.nrow();
  int ncolC = C.ncol();
  int nrowC = C.nrow();
  
  // Output matrix dimensions
  int nrowRes = nrowA;
  int ncolRes = ncolB * ncolC;
  NumericMatrix res_matr(nrowRes, ncolRes);
  
  int col_idx = 0;
  
  for (int jB = 0; jB < ncolB; ++jB) {
    
    for (int jC = 0; jC < ncolC; ++jC) {
      // Construct right vector of size nrowB * nrowC
      NumericVector right(nrowB * nrowC);
      for (int iB = 0; iB < nrowB; ++iB) {
        double b_val = B(iB, jB);
        for (int iC = 0; iC < nrowC; ++iC) {
          right[iB * nrowC + iC] = b_val * C(iC, jC);
        }
      }
      
      // Compute each row of result
      for (int i = 0; i < nrowA; ++i) {
        double sum = 0.0;
        for (int k = 0; k < ncolA; ++k) {
          sum += A(i, k) * right[k];
        }
        res_matr(i, col_idx) = sum;
      }
      
      col_idx++;
      
      if (col_idx % 100 == 0) {
        Rcpp::checkUserInterrupt();
        Rcpp::Rcout << "Processed column " << col_idx << std::endl;
      }
    }
  }
  
  return res_matr;
}

struct KroneckerWorker : public Worker {
  const RMatrix<double> A;
  const RMatrix<double> B;
  const RMatrix<double> C;
  
  RMatrix<double> res;
  
  KroneckerWorker(const NumericMatrix& A,
                  const NumericMatrix& B,
                  const NumericMatrix& C,
                  NumericMatrix& res)
    : A(A), B(B), C(C), res(res) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    int nrowA = A.nrow();
    int ncolA = A.ncol();
    int nrowB = B.nrow();
    int ncolC = C.ncol();
    int nrowC = C.nrow();
    
    for (std::size_t col_idx = begin; col_idx < end; ++col_idx) {
      int jB = col_idx / ncolC;
      int jC = col_idx % ncolC;
      
      std::vector<double> right(nrowB * nrowC);
      
      for (int iB = 0; iB < nrowB; ++iB) {
        double b_val = B(iB, jB);
        for (int iC = 0; iC < nrowC; ++iC) {
          right[iB * nrowC + iC] = b_val * C(iC, jC);
        }
      }
      
      for (int i = 0; i < nrowA; ++i) {
        double sum = 0.0;
        for (int k = 0; k < ncolA; ++k) {
          sum += A(i, k) * right[k];
        }
        res(i, col_idx) = sum;
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix KroneckerProdByBlocksParallel(const NumericMatrix& A,
                                            const NumericMatrix& B,
                                            const NumericMatrix& C) {
  int ncolB = B.ncol();
  int ncolC = C.ncol();
  int ncolRes = ncolB * ncolC;
  int nrowRes = A.nrow();
  
  NumericMatrix res(nrowRes, ncolRes);
  
  KroneckerWorker worker(A, B, C, res);
  parallelFor(0, ncolRes, worker);
  
  return res;
}





