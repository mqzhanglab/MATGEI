// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List get_V3_1(arma::mat X, arma::mat G, arma::mat Y, arma::mat mu) {
  int J = mu.n_cols;
  int n = mu.n_rows;
  arma::mat Fu_matrix(n * J, n * J, arma::fill::zeros);
  
  for (int l = 0; l < J; ++l) {
    for (int t = 0; t <= l; ++t) {
      if (l == t) {
        Fu_matrix.submat(l * n, l * n, (l + 1) * n - 1, (l + 1) * n - 1).diag() = mu.col(l) % (1 - mu.col(l));
      } else {
        Fu_matrix.submat(l * n, t * n, (l + 1) * n - 1, (t + 1) * n - 1).diag() = -mu.col(l) % mu.col(t);
        Fu_matrix.submat(t * n, l * n, (t + 1) * n - 1, (l + 1) * n - 1).diag() = -mu.col(l) % mu.col(t);
      }
    }
  }
  
  
  Rcpp::List V_list(J);
  arma::mat X_bar = arma::kron(arma::eye(J - 1, J - 1), X);
  arma::mat G_bar = arma::kron(arma::eye(J - 1, J - 1), G);
  
  for (int i_j = 0; i_j < J; ++i_j) {
    arma::mat Fu;
    if (i_j == 0) {
      Fu = Fu_matrix(arma::span(n, J*n-1), arma::span(n, J*n-1));
    } else if (i_j == J - 1) {
      Fu = Fu_matrix(arma::span(0, ((J - 1) * n - 1)),arma::span(0, (J-1)*n-1));
    } else {
      arma::uvec ind1 = arma::linspace<arma::uvec>(0, (i_j*n-1),i_j*n);
      arma::uvec ind2 = arma::linspace<arma::uvec>((i_j+1)*n, (J*n-1),n*(J-i_j-1));
      arma::uvec indc = arma::join_cols(ind1, ind2);// Concatenate the indices
      // Rcpp::Rcout<<indc;
      Fu = Fu_matrix(indc,indc);
    }
    
    double epsilon = 1e-16;
    Fu.diag() += epsilon;
    
    //Rcpp::Rcout << "Is positive finite? "<< (X_bar.t() * Fu * X_bar).is_sympd()<<std::endl;
    //Rcpp::Rcout << "Is negative finite? "<<(-X_bar.t() * Fu * X_bar).is_sympd()<<std::endl;
    arma::mat L = arma::chol(Fu);
    arma::mat M = L * X_bar;
    arma::mat N = L * G_bar;
    
    arma::mat Hat = M * inv_sympd(M.t() * M) * M.t();
    Hat.diag()=Hat.diag()-1;
    arma::mat cur_V = - N.t() * Hat * N;
    
    V_list[i_j] = cur_V;
    
  }
  
  return V_list;
}