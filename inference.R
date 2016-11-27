library(mvtnorm)
library(reshape2)

source("helpers.R")
Rcpp::sourceCpp("inference.cpp")

infer_latent_factors = function(Y, X, k_latent, alpha, max_iter = 100, burnin = max_iter/2, keep_Y_trace = TRUE, which_genes = NULL, which_samples = NULL){
  K = k_latent
  G = nrow(Y)
  N = ncol(Y)
  Ypred_mean = matrix(0, nrow(Y), ncol(Y))
  Z_mean = matrix(0, G, K)
  W_mean = matrix(0, K, N)
  beta_mean = rep(0, K)
  loglik = rep(0, max_iter - burnin)
  Ypred_trace = NULL
  if(keep_Y_trace){
    if(is.null(which_genes)) which_genes = sample(1:G, 5)
    if(is.null(which_samples)) which_samples = sample(1:N, 5)
    Ypred_trace = matrix(0, length(which_genes)*length(which_samples), max_iter - burnin)
  }
  
  a0 = 0.01
  b0 = 0.01
  # initialise precision parameters
  rho = 0.01
  lambda = 1
  gamma = rep(1, G)
  # initialise Z
  Z = matrix(0, G, K)
  for(k in 1:K){
    Z[, k] = rbinom(G, 1, 0.05)
  }
  # initialise W
  W = matrix(rnorm(K*N), K, N)
  # initialise beta
  beta = rnorm(K, 0, 1)
  
  for(i in 1:max_iter){
    # update Z
    Z = update_Z(Y, Z, W, gamma, G, K, N, alpha)
    # a quick fix: numerical problems occur when Z has one empty latent factor
    if(any(colSums(Z) == 0)){
      for(kk in which(colSums(Z) == 0))
        Z[sample(1:G, 1), kk] = 1
    }
    
    # update W
    sigma_W_inv = lambda * diag(K) + matrix_list_sum(lapply(1:G, function(g){
      gamma[g] * outer(Z[g, ], Z[g, ])
    }))
    sigma_W = solve(sigma_W_inv)
    temp_W = lambda * outer(beta, colSums(X)) + matrix_list_sum(lapply(1:G, function(g){
      gamma[g] * outer(Z[g, ], Y[g, ])
    }))
    mu_W = sigma_W %*% temp_W
    noise = t(rmvnorm(N, rep(0, K), sigma_W))
    W = mu_W + noise
    
    # Ypred
    Ypred0 = Z %*% W
    
    # update beta
    x_colsums = colSums(X)
    sigma_beta = 1 / (rho + sum(lambda * x_colsums**2))
    mu_beta = sigma_beta * colSums(lambda * x_colsums * t(W))
    beta = mu_beta + rnorm(K, 0, sqrt(sigma_beta))
    
    # update precision parameters
    a_lambda = a0 + 0.5*K*N
    b_lambda = b0 + 0.5*sum((W - outer(beta, x_colsums))**2)
    lambda = rgamma(1, a_lambda, b_lambda)
    a_gamma = a0 + 0.5*N
    b_gamma = b0 + 0.5*rowSums((Y - Ypred0)**2)
    gamma = rgamma(G, a_gamma, b_gamma)
    
    # loglikelihood
    if(i > burnin){
      m_k = colSums(Z)
      loglik_Z = sum(log(alpha/K) + lgamma(m_k + alpha/K) + lgamma(G - m_k + 1) - lgamma(N + 1 + alpha/K))
      loglik_W = sum(dnorm(W, outer(beta, x_colsums), 1/sqrt(lambda), log=TRUE))
      loglik_Y = sum(dnorm(Y, Ypred0, 1/sqrt(gamma), log=TRUE))
      loglik[i-burnin] = loglik_Z + loglik_W + loglik_Y
    }
    
    # save traces
    if(i > burnin){
      Z_mean = Z_mean + 1/(max_iter-burnin) * Z
      W_mean = W_mean + 1/(max_iter-burnin) * W
      beta_mean = beta_mean + 1/(max_iter-burnin) * beta
      Ypred_mean = Ypred_mean + 1/(max_iter-burnin) * Ypred0
      if(keep_Y_trace) Ypred_trace[, i-burnin] = melt(Ypred0[which_genes, which_samples, drop=FALSE])$value
      # cat(sprintf("range gamma (%1.3f, %1.3f), mean(gamma): %1.3f, lambda: %1.3f\n", min(gamma), max(gamma), mean(gamma), lambda))
    }
  }
  if(keep_Y_trace) Ypred_trace = postprocess_Y_trace(Y, Ypred_trace, which_genes, which_samples)
  
  return(list(Yobs = Y, Z = Z_mean, W = W_mean, beta = beta, Ypred = Ypred_mean, Ypred_trace = Ypred_trace, loglik = loglik, iter = (burnin+1):max_iter))
}
