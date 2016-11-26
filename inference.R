library(mvtnorm)

matrix_list_sum = function(lst) Reduce("+", lst)

infer_latent_factors = function(Y, X, k_latent, gamma, max_iter = 100, burnin = max_iter/2){
  K = k_latent
  G = nrow(Y)
  N = ncol(Y)
  Ypred = matrix(0, nrow(Y), ncol(Y))
  
  # gamma = rep(0.1, G)
  rho = 0.01
  lambda = 1
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
    Z = update_Z(Y, Z, W, gamma, G, K, N)
    
    # update W
    sigma_W_inv = lambda * diag(K) + matrix_list_sum(lapply(1:G, function(g){
      gamma[g] * outer(Z[g, ], Z[g, ])
    }))
    sigma_W = solve(sigma_W_inv)
    temp_W = lambda * beta * colSums(X) + matrix_list_sum(lapply(1:G, function(g){
      gamma[g] * outer(Z[g, ], Y[g, ])
    }))
    mu_W = sigma_W %*% temp_W
    noise = t(rmvnorm(N, rep(0, K), sigma_W))
    W = mu_W + noise
    
    # update beta
    x_colsums = colSums(X)
    sigma_beta = 1 / (rho + sum(lambda * x_colsums**2))
    mu_beta = sigma_beta * colSums(lambda * x_colsums * t(W))
    beta = mu_beta + rnorm(K, 0, sqrt(sigma_beta))
    
    # update precision parameters
    
    # prediction
    if(i > burnin) Ypred = Ypred + 1/(max_iter-burnin) * Z %*% W
  }
  return(list(Z = Z, W = W, beta = beta, Ypred = Ypred))
}
