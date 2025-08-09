
get_pvalue <- function(df, down_ratio="down_ratio", up_ratio="up_ratio", methods = "two.sided"){
  
  if (!(down_ratio %in% names(df)) || !(up_ratio %in% names(df))) {
    stop("Specified columns not found in the dataframe.")
  }
  
  p_calculate <- df %>% dplyr::select(all_of(c(down_ratio, up_ratio))) %>%
    apply(.,1,function(x) (x / sum(x)) * 100) %>% t() %>% as.data.frame()
  
  pvalue <- t.test(p_calculate[[down_ratio]], p_calculate[[up_ratio]], paired = TRUE, alternative = methods)
  
  # df$down_ratio <- p_calculate$down_ratio
  # df$up_ratio <- p_calculate$up_ratio
  
  return(pvalue$p.value)
}