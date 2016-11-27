library(dplyr)

matrix_list_sum = function(lst) Reduce("+", lst)

postprocess_Y_trace = function(Yobs, Ypred_trace, which_genes, which_samples){
  temp_trace = data.frame(gene = rep(which_genes, length(which_samples)), 
                          sample = rep(which_samples, each=length(which_genes)), 
                          yobs = melt(Yobs[which_genes, which_samples, drop=FALSE])$value, 
                          Ypred_trace)
  df = melt(temp_trace, 
            id.vars = c("gene", "sample", "yobs"), 
            value.name = "ypred", 
            variable.name = "iter") %>%
    mutate(iter = as.numeric(substr(iter, 2, 10)))
  return(df)
}
