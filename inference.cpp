#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

arma::ivec colSums(const arma::imat & X){
  int nCols = X.n_cols;
  arma::ivec out(nCols);
  for(int i = 0; i < nCols; i++){
    out(i) = sum(X.col(i));
  }
  return(out);
}

// [[Rcpp::export]]
bool calculate_accept_prob(int g, int k, arma::mat& Z, arma::mat& W, arma::mat& Y, NumericVector gamma, NumericVector mu, NumericVector lambda, int G, int K, int N, int count, double alpha){
  arma::rowvec pred0(N), pred1(N);
  // leave out latent factor k from prediction Y_pred(g, _)
  Z(g, k) = 0;
  pred0 = Z.row(g) * W + mu[g];
  // compute \sum_{i=1}^N W_{ki}^2
  double sum_W2 = sum(square(W.row(k)));
  double denom = lambda[k] + gamma[g] * sum_W2;
  double log_det = log(lambda[k]) - log(denom);
  double mean = gamma[g] / denom * dot(Y.row(g) - pred0, W.row(k));
  double precision = denom;
  double log_lik = precision * mean * mean;
  double logz = log((count + alpha/K) / (G - count));
  double z = logz + 0.5*log_det + 0.5*log_lik;
  double accept_prob = 1.0 / (1.0 + exp(-z));
  //printf("accept prob: %1.3f, mean: %1.3f, precision: %1.3f\n", accept_prob, mean, precision);
  
  double u = R::runif(0, 1);
  bool accepted = (u < accept_prob);
  if(accepted){
    Z(g, k) = R::rnorm(mean, 1/sqrt(precision));
  }
  return accepted;
}

// [[Rcpp::export]]
arma::mat update_Z(arma::mat& Y, arma::imat Z, arma::mat ZZ, arma::mat & W, NumericVector gamma, NumericVector mu, NumericVector lambda, int G, int K, int N, double alpha) {
  //arma::imat Z(Z0.begin(), G, K, false);
  //arma::ivec genes_per_pathway = colSums(Z);
  
  arma::ivec genes_per_pathway = colSums(Z);
  for(int k=0; k<K; k++){
    for(int g=0; g<G; g++){
      // leave out gene g from the count
      int count = genes_per_pathway[k] - Z(g, k);
      bool accepted =  calculate_accept_prob(g, k, ZZ, W, Y, gamma, mu, lambda, G, K, N, count, alpha);
      genes_per_pathway[k] = count;
      // set Z[g, k] to 1 with probability accept_prob
      if(accepted){
        Z(g, k) = 1;
        genes_per_pathway[k] += 1;
      }
    }
  }
  return ZZ;
}

// [[Rcpp::export]]
NumericMatrix fast_outer_product_sum(NumericVector gamma, NumericMatrix A, NumericMatrix B, int nrows, int ncols, int length){
  NumericMatrix out(nrows, ncols);
  for(int j=0; j<ncols; j++){
    for(int i=0; i<nrows; i++){
      for(int g=0; g<length; g++){
        out(i, j) += gamma[g] * A(g, i) * B(g, j);
      }
    }
  }
  return out;
}
