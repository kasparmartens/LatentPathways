generate_data1 = function(N, G, K, alpha = 1){
  X = matrix(rbinom(G*N, 1, 0.1), G, N)
  pi = rbeta(K, alpha/K, 1)
  Z = matrix(0, G, K) 
  for(k in 1:K){
    Z[, k] = rbinom(G, 1, pi[k])
  }
  beta = rnorm(K, 0, 1)
  rho = 0.01
  lambda = 1
  gamma = rgamma(G, 1, 1)
  W = matrix(0, K, N)
  for(k in 1:K){
    W[k, ] = rnorm(N, beta[k] * colSums(X), sqrt(1/lambda))
  }
  Y = Z %*% W + rnorm(N*G, 0, sqrt(1/gamma))
  return(list(X = X, Y = Y, Z = Z, W = W, gamma = gamma))
}
