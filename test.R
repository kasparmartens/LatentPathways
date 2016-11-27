source("generate_data.R")
source("inference.R")

data = generate_data1(N = 100, G = 250, K = 10, alpha = 2)
res = infer_latent_factors(data$Y, data$X, data$K, data$alpha, max_iter = 5000, which_genes = 1:4, which_samples = 1:3)

library(reshape2)
library(dplyr)
library(ggplot2)

y.m = melt(data$Y)
ypred.m = melt(res$Ypred)
df = data.frame(gene = y.m$Var1, y = y.m$value, ypred = ypred.m$value)

ggplot(subset(df, gene %in% sample(unique(gene), 100)), aes(y, ypred)) + 
  geom_point() + 
  geom_abline(intercept=0, slope=1, col="red", linetype="dashed") +
  theme_bw()



ggplot(res$Ypred_trace) + 
  # geom_abline(intercept=0, slope=1, col="red", linetype="dashed") +
  facet_wrap(~ paste("gene", gene) + paste("ind", sample), scales="free") + 
  geom_line(aes(x = iter, y = ypred))+
  geom_hline(aes(yintercept = yobs), col="red", linetype="dashed") + 
  theme_bw()


# Ypred_trace = do.call("rbind", lapply(seq_along(res$Ypred_trace), function(i){
#   mat = res$Ypred_trace[[i]]
#   melt(mat) %>%
#     mutate(iter = i, y = y.m[["value"]]) %>%
#     rename(gene = Var1, ind = Var2, ypred = value)
# }))
