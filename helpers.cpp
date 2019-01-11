// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::field<uvec> gen_k_ind(arma::vec x, int K) {
  field<uvec> ans(K);
  // The 0 case needs to be explicitly handled
  // because when we re-scale our covariates 
  // we will have a value which will be exactly 0. 
  ans(0) = find(K * x <= 1 && K * x >= 0);
  
  for(int i = 1; i < K; i = i + 1) {
    ans(i) = find(K * x <= (i+1) && K * x > i);
  }
  return ans;
}

// [[Rcpp::export]]
arma::vec gen_k_vec(arma::vec x, int K) {
  vec ans(x.n_elem);
  
  // The 0 case needs to be explicitly handled
  // because when we re-scale our covariates 
  // we will have a value which will be exactly 0. 
  uvec temp = find(K * x <= 1 && K * x >= 0);
  ans(temp) = zeros(temp.n_elem);
  
  for(int i = 1; i < K; i = i + 1) {
    temp = find(K * x <= (i+1) && K * x > i);
    ans(temp) = i * ones(temp.n_elem);
  }
  return ans;
}
// [[Rcpp::export]]
arma::mat gen_r(arma::vec k_vec, arma::vec x,
                int K) {
  // We create two vectors and stroe them in a matrix
  // This matrix will be n * K, where n = length(x).
  //
  //  Args:
  //    k_vec: A vector denoting membership in intervals.
  //    x: Data vector, assume to be sorted.
  //    K: Integer specifying the number of intervals.
  
  vec a = 1 - (K*x - k_vec);
  vec b = K*x - k_vec;
  
  uvec ind = find(k_vec == 0);
  a.elem(ind) = ones<vec>(ind.n_elem);
  b.elem(ind) = zeros<vec>(ind.n_elem);
  
  return join_rows(a, b);
}

// [[Rcpp::export]]
arma::vec rty(arma::mat sp_R, arma::vec y,
              arma::field<uvec> k_ind, int K) {
  // We create the vector given by t(R) * y
  // where R is stored in special sparse form , i.e. as two vectors.
  //
  //  Args:
  //    sp_R: n * 2 matrix used to represent sparse matrix R.
  //    y: Data vector, assume to be sorted according to x.
  //    k_vec: A vector denoting membership in intervals.
  //    K: Integer specifying the number of intervals.
  
  vec ans(K,fill::zeros);
  vec ay = sp_R.col(0) % y;
  vec by = sp_R.col(1) % y;
  ans(0) = accu(ay.elem(k_ind(0))) + accu(ay.elem(k_ind(1)));
  for( int i = 1; i < K-1; i = i + 1 ) {
    ans(i) = accu(by.elem(k_ind(i))) + 
      accu(ay.elem(k_ind(i+1)));
  }
  ans(K-1) = accu(by.elem(k_ind(K-1)));
  return ans;
}

// [[Rcpp::export]]
arma::mat rtr(arma::mat sp_R,
               arma::field<uvec> k_ind, int K) {
  // We create the vector given by t(R) * y
  // where R is stored in special sparse form , i.e. as two vectors.
  //
  //  Args:
  //    sp_R: n * 2 matrix used to represent sparse matrix R.
  //    k_ind: A field with indices denoting membership in intervals.
  //    K: Integer specifying the number of intervals.
  
  vec ans(K);
  vec ans2(K);
  
  vec a = sp_R.col(0);
  vec b = sp_R.col(1);
  ans(0) = accu(square(a.elem(k_ind(0)))) + 
    accu(square(a.elem(k_ind(1))));
  ans2(0) = 0;
  for( int i = 1; i < K - 1; i = i + 1 ) {
    ans(i) = accu(square(b.elem( k_ind(i) ))) + 
      accu(square(a.elem( k_ind(i + 1) )));
    
    ans2(i) = accu(a.elem(k_ind(i)) % b.elem(k_ind(i)));
  }
  ans(K - 1) = accu(square(b.elem( k_ind(K - 1) )));
  ans2(K - 1) = accu(a.elem(k_ind(K - 1)) % b.elem(k_ind(K - 1)));
  return join_rows(ans, ans2);
}

// [[Rcpp::export]]
arma::vec max_eigen(arma::mat rtr_mat, int K) {
  sp_mat ans(K,K);
  ans.diag() = rtr_mat.col(0);
  
  vec temp = rtr_mat.col(1);
  vec temp2 = temp.tail(K - 1);
  ans.diag(1) = temp2;
  ans.diag(-1) = temp2;
 
 return eigs_sym(ans, 1); 
}

// [[Rcpp::export]]
arma::vec tridiag_y(arma::mat tri_diag, arma::vec y, int K) {
  vec ans(K);
  ans(0) = tri_diag(0, 0) * y(0) + tri_diag(1, 1) * y(1);
  
  for(int i = 1; i < K - 1; i = i + 1) {
    ans(i) = tri_diag(i+1 - 1, 1) * y(i-1) + 
      tri_diag(i, 0) * y(i) + tri_diag(i+1, 1) * y(i+1);
  }
  
  ans(K - 1) = tri_diag(K - 1, 1) * y(K-2) + 
    tri_diag(K-1, 0) * y(K-1);
  return ans;
}

// [[Rcpp::export]]
arma::mat rtheta(arma::mat spR, arma::mat theta,
                 arma::field<uvec> k_ind, int K, int n) {
  mat ans(n, theta.n_cols);
  uvec ind0 = k_ind(0);
  ans.rows(ind0) = ones(ind0.n_elem) * theta.row(0);
  
  for(int i = 1; i < K; i = i + 1) {
    uvec indi = k_ind(i);
    mat temp = spR.rows(indi);

    ans.rows(indi) = temp * theta.rows(i - 1, i);
  }
  return ans;
}
  


