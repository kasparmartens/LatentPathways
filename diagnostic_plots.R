library(reshape2)
library(dplyr)
library(ggplot2)

plot_log_posterior = function(obj){
  df = data.frame(iter = obj$iter, loglik = obj$loglik)
  p = ggplot(df, aes(iter, loglik)) + 
    geom_line() + 
    theme_bw() + labs(x = "iteration", y = "log posterior", title="log posterior for training data")
  return(p)
}

scatterplot_observed_vs_predicted = function(obj, alpha = 1, type = "train", faceting = FALSE, n_genes = 50){
  if(type == "train"){
    y.m = melt(obj$Yobs)
    ypred.m = melt(obj$Ypred)
  } else{
    y.m = melt(obj$Yobs_test)
    ypred.m = melt(obj$Ypred_test)
  }

  df = data.frame(gene = y.m$Var1, y = y.m$value, ypred = ypred.m$value)
  if(nrow(df) > 1e5) df = sample_n(df, 1e5)
  
  if(faceting){
    # show the whole (sub)dataset on the background
    if(nrow(df) > 10000){
      df_helper = df %>% sample_n(10000) %>% transform(gene=NULL)
    } else{
      df_helper = df %>% transform(gene=NULL)
    }
    # pick randomly 50 genes
    df = df %>% 
      filter(gene %in% sample(unique(gene), n_genes)) %>%
      group_by(gene) %>%
      mutate(corrcoef = cor(y, ypred), 
             label = sprintf("cor = %1.2f (gene %s)", corrcoef, gene)) %>%
      arrange(desc(corrcoef)) 
    
    p = ggplot(df, aes(y, ypred)) + 
      geom_point(alpha = 0.2, data = df_helper, col="grey80") + 
      geom_point(alpha = alpha) + 
      geom_abline(intercept=0, slope=1, col="red", linetype="dashed") +
      theme_bw() + facet_wrap(~ label) + 
      labs(x = "Observed y", y = "Predicted y")
  } else{
    p = ggplot(df, aes(y, ypred)) + 
      geom_point(alpha = alpha) + 
      geom_abline(intercept=0, slope=1, col="red", linetype="dashed") +
      theme_bw() + 
      labs(x = "Observed y", y = "Predicted y")
  }
  return(p)
}

traceplot_observed_vs_predicted = function(obj, type = "train"){
  if(type == "train"){
    Ypred_trace = obj$Ypred_trace
  } else{
    Ypred_trace = obj$Ypred_test_trace
  }
  
  p = ggplot(Ypred_trace) + 
    facet_grid(paste("gene", gene) ~ paste("sample", sample), scales="free") + 
    geom_line(aes(x = iter, y = ypred, col = factor(gene)))+
    geom_hline(aes(yintercept = yobs), col="red", linetype="dashed") + 
    theme_bw() + theme(legend.position = "none") +
    labs(x = "Iteration", y = "Predicted (and observed) y")
  return(p)
}

histogram_corr_coef = function(obj, type = "train"){
  if(type == "train"){
    Y = obj$Yobs
    Ypred = obj$Ypred
  } else{
    Y = obj$Yobs_test
    Ypred = obj$Ypred_test
  }
  
  n_genes = nrow(Y)
  coefs = rep(NA, n_genes)
  for(i in 1:n_genes){
    coefs[i] = cor(Y[i, ], Ypred[i, ])
  }
  p = ggplot(data.frame(x = coefs), aes(x)) + geom_histogram() + 
    theme_bw() + labs(x = "cor(y_observed, y_predicted) for every gene")
  return(p)
}
