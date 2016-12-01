library(mvtnorm)
library(reshape2)

source("helpers.R")
Rcpp::sourceCpp("inference.cpp")

calculate_IBP_prob = function(alpha, m_k, K, G, log = FALSE){
  logprob = sum(log(alpha/K) + lgamma(m_k + alpha/K) + lgamma(G - m_k + 1) - lgamma(G + 1 + alpha/K))
  if(log) return(logprob) else return(exp(logprob))
}

infer_latent_factors = function(Y, X, Y_test, X_test, K, alpha = 5, max_iter = 100, burnin = max_iter/2, keep_Y_trace = TRUE, which_genes = NULL){
  if(is.data.frame(Y)) Y = as.matrix(Y)
  if(is.data.frame(X)) X = as.matrix(X)
  if(is.data.frame(Y_test)) Y_test = as.matrix(Y_test)
  if(is.data.frame(X_test)) X_test = as.matrix(X_test)
  G = nrow(Y)
  N = ncol(Y)
  Ypred_mean = matrix(0, G, N)
  Z_mean = matrix(0, G, K)
  ZZ_mean = matrix(0, G, K)
  W_mean = matrix(0, K, N)
  beta_mean = rep(0, K)
  beta0_mean = rep(0, K)
  intercept_mean = rep(0, G)
  loglik_trace = rep(0, max_iter - burnin)
  loglik_test_trace = rep(0, max_iter - burnin)
  alpha_trace = rep(0, max_iter - burnin)
  Ypred_trace = NULL
  # for the test set
  Ntest = ncol(Y_test)
  Ypred_test_mean = matrix(0, G, Ntest)
  Ypred_test_trace = NULL
  
  if(keep_Y_trace){
    if(is.null(which_genes)) which_genes = sample(1:G, 5)
    which_samples1 = sample(1:N, 5)
    which_samples2 = sample(1:Ntest, 5)
    Ypred_trace = matrix(0, length(which_genes)*length(which_samples1), max_iter - burnin)
    Ypred_test_trace = matrix(0, length(which_genes)*length(which_samples2), max_iter - burnin)
  }
  
  a0 = 0.1
  b0 = 0.1
  rho = 0.01
  # initialise precision parameters
  lambda_W = 1
  gamma = rep(1, G)
  lambda = rep(1, K)
  # initialise Z
  Z = matrix(0, G, K)
  for(k in 1:K){
    Z[, k] = rbinom(G, 1, 0.05)
  }
  ZZ = Z * rnorm(G*K)
  # initialise W
  W = matrix(rnorm(K*N), K, N)
  # initialise beta
  beta = rnorm(K, 0, 1)
  beta0 = rnorm(K, 0, 1)
  # initialise mu_Y
  intercept_g = rep(0, G)
  
  x_colsums = colSums(X)
  x_test_colsums = colSums(X_test)

  for(i in 1:max_iter){
    # update Z
    ZZ = update_Z(Y, Z, ZZ, W, gamma, intercept_g, lambda, G, K, N, alpha)
    Z = 1*(ZZ != 0)
    
    # a quick fix: numerical problems occur when Z has one empty latent factor
    if(any(colSums(Z) == 0)){
      for(kk in which(colSums(Z) == 0))
        Z[sample(1:G, 1), kk] = 1
    }
    m_k = colSums(Z)
    
    # update alpha via RW-MH
    alpha_proposed = exp(log(alpha) + rnorm(1, 0, 0.1))
    accept_num = log(alpha_proposed) + dgamma(alpha_proposed, a0, b0, log = TRUE) + calculate_IBP_prob(alpha_proposed, m_k, K, G, TRUE)
    accept_denom = log(alpha) + dgamma(alpha, a0, b0, log = TRUE) + calculate_IBP_prob(alpha, m_k, K, G, TRUE)
    acceptance_prob = exp(accept_num - accept_denom)
    if(runif(1) < acceptance_prob){
      alpha = alpha_proposed
    }
    
    # update lambda
    a_lambda = a0 + 0.5*N
    b_lambda = b0 + 0.5*colSums(ZZ**2)
    lambda = rgamma(K, a_lambda, b_lambda)
    
    # update W
    sigma_W_inv = lambda_W * diag(K) + fast_outer_product_sum(gamma, Z, Z, K, K, G)
    sigma_W = solve(sigma_W_inv)
    temp_W = lambda_W * outer(beta0, rep(1, N)) + lambda_W * outer(beta, colSums(X)) + fast_outer_product_sum(gamma, Z, Y, K, N, G)
    mu_W = sigma_W %*% temp_W
    noise = t(rmvnorm(N, rep(0, K), sigma_W))
    W = mu_W + noise
    
    # intercept_g
    resid = Y - ZZ %*% W # note without the intercept
    sigma_int = 1/(rho + N * gamma)
    mu_int = sigma_int * rowSums(gamma * resid)
    intercept_g = rnorm(G, mu_int, sqrt(sigma_int))

    # Ypred without noise
    Ypred0 = ZZ %*% W + intercept_g
    Ypred = Ypred0 + rnorm(G*N, 0, 1/sqrt(gamma))
    
    # update beta
    sigma_beta = 1 / (rho + sum(lambda_W * x_colsums**2))
    mu_beta = sigma_beta * colSums(lambda_W * x_colsums * t(W))
    beta = mu_beta + rnorm(K, 0, sqrt(sigma_beta))
    # beta0
    sigma_beta = 1 / (rho + lambda_W * N)
    mu_beta = sigma_beta * rowSums(lambda_W * (W - outer(beta, x_colsums)))
    beta0 = mu_beta + rnorm(K, 0, sqrt(sigma_beta))
    
    # update precision parameters
    a_lambda_W = a0 + 0.5*K*N
    b_lambda_W = b0 + 0.5*sum((W - beta0 - outer(beta, x_colsums))**2)
    lambda_W = rgamma(1, a_lambda_W, b_lambda_W)
    a_gamma = a0 + 0.5*N
    b_gamma = b0 + 0.5*rowSums((Y - Ypred0)**2)
    gamma = rgamma(G, a_gamma, b_gamma)
    
    # save traces
    if(i > burnin){
      alpha_trace[i-burnin] = alpha
      Z_mean = Z_mean + 1/(max_iter-burnin) * Z
      ZZ_mean = ZZ_mean + 1/(max_iter-burnin) * ZZ
      W_mean = W_mean + 1/(max_iter-burnin) * W
      intercept_mean = intercept_mean + 1/(max_iter-burnin) * intercept_g
      beta_mean = beta_mean + 1/(max_iter-burnin) * beta
      beta0_mean = beta0_mean + 1/(max_iter-burnin) * beta0
      Ypred_mean = Ypred_mean + 1/(max_iter-burnin) * Ypred
      if(keep_Y_trace) Ypred_trace[, i-burnin] = melt(Ypred[which_genes, which_samples1, drop=FALSE])$value
      # cat(sprintf("range gamma (%1.3f, %1.3f), mean(gamma): %1.3f, lambda_W: %1.3f\n", min(gamma), max(gamma), mean(gamma), lambda_W))
      
      # loglikelihood
      loglik_Z = calculate_IBP_prob(alpha, m_k, K, G, log = TRUE)
      loglik_W = sum(dnorm(W, beta0 + outer(beta, x_colsums), 1/sqrt(lambda_W), log=TRUE))
      loglik_Y = sum(dnorm(Y, Ypred0, 1/sqrt(gamma), log=TRUE))
      loglik_trace[i-burnin] = loglik_Z + loglik_W + loglik_Y
      
      
      # predict Y_test
      W_test = beta0 + outer(beta, x_test_colsums) + rnorm(K*Ntest, 0, 1/sqrt(lambda_W))
      Ypred0_test = intercept_g + ZZ %*% W_test
      Ypred_test = Ypred0_test + rnorm(G*Ntest, 0, 1/sqrt(gamma))
      Ypred_test_mean = Ypred_test_mean + 1/(max_iter-burnin) * Ypred_test
      if(keep_Y_trace) Ypred_test_trace[, i-burnin] = melt(Ypred_test[which_genes, which_samples2, drop=FALSE])$value
      
      # test loglikelihood
      loglik_W = sum(dnorm(W_test, beta0 + outer(beta, x_test_colsums), 1/sqrt(lambda_W), log=TRUE))
      loglik_Y = sum(dnorm(Y_test, Ypred0_test, 1/sqrt(gamma), log=TRUE))
      loglik_test_trace[i-burnin] = loglik_Z + loglik_W + loglik_Y
    }
    
    if(i %% 100 == 0){
      cat("Iter", i, "\n")
      if(i > burnin){
        train_acc = cor(as.numeric(Ypred), as.numeric(Y))
        test_acc = cor(as.numeric(Ypred_test), as.numeric(Y_test))
        cat("\tTrain acc:", train_acc, "Test acc", test_acc, "\n")
      }
      flush.console()
    }
  }
  if(keep_Y_trace){
    Ypred_trace = postprocess_Y_trace(Y, Ypred_trace, which_genes, which_samples1)
    Ypred_test_trace = postprocess_Y_trace(Y, Ypred_test_trace, which_genes, which_samples2)
  }
  
  acc_train = cor(as.numeric(Ypred), as.numeric(Y))
  acc_test = cor(as.numeric(Ypred_test), as.numeric(Y_test))
  
  return(list(Yobs = Y, Yobs_test = Y_test, 
              acc_train = acc_train, acc_test = acc_test, 
              Z = Z_mean, ZZ = ZZ_mean, 
              W = W_mean, 
              beta = beta_mean, beta0 = beta0_mean, 
              Ypred = Ypred_mean, Ypred_trace = Ypred_trace, 
              Ypred_test = Ypred_test_mean, Ypred_test_trace = Ypred_test_trace, 
              loglik = loglik_trace, loglik_test = loglik_test_trace, 
              alpha = alpha_trace, 
              intercept_g = intercept_g, 
              iter = (burnin+1):max_iter))
}
