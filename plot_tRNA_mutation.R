
plot_tRNA_mutation <- function(input_file="./7.plot/all_tRNA_mutation.xlsx", sheet=6, treat = "t6mo_mean_diff"){
  library(ggplot2)
  library(viridis)
  library(dplyr)
  library(readxl)
  library(stringr)
  
  # data clean
  df <- read_excel(input_file, sheet = sheet) %>%
    dplyr::select(c(amino_acid_codon, t6mo_mean_diff ,t61mo_mean_diff, t6t61mo_mean_diff)) %>% 
    na.omit()
  df$codon <- gsub("_","-",df$amino_acid_codon)
  
  df <- df %>% dplyr::select(c(codon, all_of(treat))) %>% 
    rename(mean_diff = all_of(treat)) %>% 
    mutate(mean_diff = mean_diff * 100) %>% 
    dplyr::filter(mean_diff > 0)
  
  df$mean_diff <- -df$mean_diff
  
  name <- stringr::str_split(treat, pattern = "_") %>% unlist() %>% .[1]
  
  # plot
  p1 <- ggplot(df, 
               aes(x = reorder(codon, 
                               mean_diff), 
                   y = mean_diff)) +
    geom_vline(aes(xintercept = as.numeric(factor(codon))),
               linetype = "dashed",
               color = "gray", 
               alpha = 0.6) +
    geom_point(aes(fill = mean_diff), 
               shape = 21, 
               alpha = 1,
               size = 2.8) +
    scale_fill_viridis(option = "D",
                       limits = c(min(df$mean_diff),
                                  max(df$mean_diff)),
                       name = "Diff",
                       direction = -1) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7, color = "black", face = "bold"),
      axis.text.y = element_text(size = 12, color = "black"),
      legend.title = element_text(size = 10)
    ) +
    labs(title = paste0("m1A change of ", name, " KO"),
         x = "tRNA",
         y = "tRNA m1A level Diff (KO-WT) %") +
    theme(panel.border = element_rect(color = "black", 
                                      fill = NA, 
                                      linewidth = 1))
  
  print(p1)
  ggsave(paste0("./7.plot/merge_isoform_mean_diff_", name, ".pdf"), p1, width = 10, height = 5)
}
