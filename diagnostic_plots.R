library(reshape2)
library(ggplot2)

plot_log_posterior = function(obj){
  df = data.frame(iter = obj$iter, loglik = obj$loglik)
  p = ggplot(df, aes(iter, loglik)) + 
    geom_line() + 
    theme_bw() + labs(x = "iteration", y = "log posterior")
  return(p)
}

scatterplot_observed_vs_predicted = function(obj, alpha = 1){
  y.m = melt(obj$Yobs)
  ypred.m = melt(obj$Ypred)
  df = data.frame(gene = y.m$Var1, y = y.m$value, ypred = ypred.m$value)
  
  p = ggplot(df, aes(y, ypred)) + 
    geom_point(alpha = alpha) + 
    geom_abline(intercept=0, slope=1, col="red", linetype="dashed") +
    theme_bw() + 
    labs(x = "Observed y", y = "Predicted y")
  return(p)
}

traceplot_observed_vs_predicted = function(obj){
  p = ggplot(obj$Ypred_trace) + 
    facet_grid(paste("gene", gene) ~ paste("sample", sample), scales="free") + 
    geom_line(aes(x = iter, y = ypred, col = factor(gene)))+
    geom_hline(aes(yintercept = yobs), col="red", linetype="dashed") + 
    theme_bw() + theme(legend.position = "none") +
    labs(x = "Iteration", y = "Predicted (and observed) y")
  return(p)
}


