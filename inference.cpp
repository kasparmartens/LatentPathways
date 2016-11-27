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
double calculate_accept_prob(int g, int k, arma::imat& Z, arma::mat& W, arma::mat& Y, NumericVector gamma, int G, int K, int N, int count, double alpha){
  arma::rowvec pred0(N), pred1(N);
  // leave out latent factor k from prediction Y_pred(g, _)
  if(Z(g, k) == 0){
    pred0 = Z.row(g) * W;
    pred1 = pred0 + W.row(k);
  } else{
    pred1 = Z.row(g) * W;
    pred0 = pred1 - W.row(k);
  }
  arma::rowvec y = Y.row(g);
  double loglik0 = sum(-square(y-pred0)/gamma[g]);
  double loglik1 = sum(-square(y-pred1)/gamma[g]);
  double logz = log((count + alpha/K) / (G - count));
  double z = logz + loglik1 - loglik0;
  double accept_prob = 1.0 / (1.0 + exp(-z));
  //printf("logz: %1.2f, loglik0: %1.2f, loglik1: %1.2f, z: %1.2f, prob: %1.3f\n", logz, loglik0, loglik1, z, accept_prob);
  return accept_prob;
}

// [[Rcpp::export]]
arma::imat update_Z(arma::mat& Y, arma::imat Z, arma::mat & W, NumericVector gamma, int G, int K, int N, double alpha) {
  //arma::imat Z(Z0.begin(), G, K, false);
  //arma::ivec genes_per_pathway = colSums(Z);
  
  arma::ivec genes_per_pathway = colSums(Z);
  for(int k=0; k<K; k++){
    for(int g=0; g<G; g++){
      // leave out gene g from the count
      int count = genes_per_pathway[k] - Z(g, k);
      double accept_prob =  calculate_accept_prob(g, k, Z, W, Y, gamma, G, K, N, count, alpha);
      genes_per_pathway[k] = count;
      // set Z[g, k] to 1 with probability accept_prob
      double u = R::runif(0, 1);
      if(u < accept_prob){
        Z(g, k) = 1;
        genes_per_pathway[k] += 1;
      } else{
        Z(g, k) = 0;
      }
    }
  }
  return Z;
}



