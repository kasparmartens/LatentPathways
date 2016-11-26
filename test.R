source("generate_data.R")
source("inference.R")

Rcpp::sourceCpp("inference.cpp")

data = generate_data1(N = 100, G = 500, K = 20, alpha = 1)
res = infer_latent_factors(data$Y, data$X, 20, data$gamma, max_iter = 300)



library(reshape2)
library(ggplot2)

y.m = melt(data$Y)
ypred.m = melt(res$Ypred)
df = data.frame(gene = y.m$Var1, y = y.m$value, ypred = ypred.m$value)

ggplot(subset(df, gene %in% sample(unique(gene), 100)), aes(y, ypred)) + 
  geom_point() + 
  geom_abline(intercept=0, slope=1, col="red", linetype="dashed") +
  theme_bw()

ggplot(subset(df, gene %in% sample(unique(gene), 100)), aes(y, ypred)) + 
  geom_point() + 
  geom_abline(intercept=0, slope=1, col="red", linetype="dashed") +
  facet_wrap(~ gene, scales="free") + 
  theme_bw()
