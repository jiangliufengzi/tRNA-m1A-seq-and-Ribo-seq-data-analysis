
# set colors
mycol <- setNames(c("#EB4232", "#2E9FDF", "#d8d8d8"), c("up", "down", "none"))

# set theme
mytheme <- theme_classic() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.margin = margin(15,5.5,5.5,5.5))

##################################################
# plot function
##################################################
plot_diff_tRNA_exp <- function(input_file="./7.plot/all_tRNA_mutation.xlsx", sheet=7, 
                               control1="ctr_1", control2="ctr_2", 
                               treat1="t6mo_1", treat2="t6mo_2",
                               cut_log2FC=1, cut_pvalue=0.05,
                               out_file="tRNA_mutation"){
  library(ggplot2)
  library(ggrepel)
  library(readxl)
  library(tidyverse)
  
  # 检查目录是否存在
  if (!dir.exists("./7.plot/")) {
    dir.create("./7.plot/", recursive = TRUE)
  }
  
  # 检查输入文件是否存在
  if (!file.exists(input_file)) {
    stop("Input file not found!")
  }
  
  # read file and filter low expression tRNA
  df <- read_excel(input_file, sheet = sheet) %>% mutate(sum=rowSums(.[-1])) %>%
    dplyr::filter(sum > 1000) %>% dplyr::select(-sum)
  
  df_sub <- df %>% dplyr::select(c("seq_name", all_of(control1), all_of(control2), all_of(treat1), all_of(treat2)))
  
  # caculate log2FoldChange and pvalue
  df_sub$logFC <- log2((df_sub[[control1]] + df_sub[[control2]]) / (df_sub[[treat1]] + df_sub[[treat2]]))
  
  df_sub$pvalue <- apply(df_sub %>% dplyr::select(all_of(c(control1, control2, treat1, treat2))), 1, function(x) {
    t.test(x[1:2], x[3:4])$p.value
  })
  
  # filter up and down regulated tRNA
  df_sub$group <- ifelse(df_sub$pvalue < cut_pvalue & df_sub$logFC > cut_log2FC, "up",
                         ifelse(df_sub$pvalue < cut_pvalue & df_sub$logFC < -cut_log2FC, "down", "none"))
  df_sub$group <- factor(df_sub$group, levels = c("up","down","none"))
  
  # num of up down gene
  up_num <- table(df_sub$group) %>% as.data.frame() %>% dplyr::filter(Var1 == "up") %>% dplyr::pull(Freq)
  down_num <- table(df_sub$group) %>% as.data.frame() %>% dplyr::filter(Var1 == "down") %>% dplyr::pull(Freq)
  
  # ggplot2
  # x_lim <- c(min(df_sub$logFC)*1.1, max(df_sub$logFC)*1.1)
  x_lim <- c(-1.5, 1.5)
  p <- ggplot(data = df_sub, aes(x = logFC, y = -log10(pvalue), color = group)) +
    geom_point(size = 1.2) +
    scale_colour_manual(name = "", values = alpha(mycol, 0.7)) +
    scale_x_continuous(limits = x_lim,
                       breaks = c(-1.5, seq(-1.5, 1.5, by = 0.5)),
                       labels = c("-1.5", as.character(seq(-1.5, 1.5, by = 0.5)))) +
    scale_y_continuous(expand = expansion(add = c(0.15, 0)),
                       limits = c(0, 5),
                       breaks = c(0, seq(0, 5, by = 1)),
                       labels = c("0", as.character(seq(0, 5, by = 1)))) +
    geom_hline(yintercept = -log10(cut_pvalue), linewidth = 0.5, color = "black", lty = "dashed") +
    geom_vline(xintercept = c(-cut_log2FC, cut_log2FC), size = 0.5, color = "black", lty = "dashed") +
    # annotate("text",x= -0.7,y= 4,label=down_num,size=6) +
    # annotate("text",x= 0.7,y= 4,label=up_num,size=6) +
    mytheme
  
  print(p)
  
  # save plot
  ggsave(paste0("./7.plot/", "volcano_",out_file, ".pdf"), plot = p, width = 3, height = 2)
}


##################################################
# plot function for alkb and no alkb
##################################################
plot_diff_tRNA_exp_multi <- function(input_file="./7.plot/all_tRNA_mutation.xlsx", sheet=7, 
                               control1="ctr_1", control2="ctr_2", control3="ctr_1a", control4="ctr_2a",
                               treat1="t6mo_1", treat2="t6mo_2", treat3="t6mo_1a", treat4="t6mo_2a",
                               cut_log2FC=1, cut_pvalue=0.05,
                               out_file="tRNA_mutation"){
  library(ggplot2)
  library(ggrepel)
  library(readxl)
  library(tidyverse)

  # 检查目录是否存在
  if (!dir.exists("./7.plot/")) {
    dir.create("./7.plot/", recursive = TRUE)
  }
  
  # 检查输入文件是否存在
  if (!file.exists(input_file)) {
    stop("Input file not found!")
  }
  
  # read file and filter low expression tRNA
  df <- read_excel(input_file, sheet = sheet) %>% mutate(sum=rowSums(.[-1])) %>% 
    dplyr::filter(sum > 1000) %>% dplyr::select(-sum)
  
  df_sub <- df %>% dplyr::select(c("seq_name", all_of(control1), all_of(control2), all_of(control3), all_of(control4),
                                   all_of(treat1), all_of(treat2), all_of(treat3), all_of(treat4)))
  
  # caculate log2FoldChange and pvalue
  df_sub$logFC <- log2((df_sub[[control1]] + df_sub[[control2]] + df_sub[[control3]] + df_sub[[control4]]) / 
                         (df_sub[[treat1]] + df_sub[[treat2]] + df_sub[[treat3]] + df_sub[[treat4]]))
  
  df_sub$pvalue <- apply(df_sub %>% dplyr::select(all_of(c(control1,control2,control3,control4,treat1,treat2,treat3,treat4))), 1, function(x) {
    t.test(x[1:4], x[5:8])$p.value
  })
  
  # filter up and down regulated tRNA
  df_sub$group <- ifelse(df_sub$pvalue < cut_pvalue & df_sub$logFC > cut_log2FC, "up",
                         ifelse(df_sub$pvalue < cut_pvalue & df_sub$logFC < -cut_log2FC, "down", "none"))
  df_sub$group <- factor(df_sub$group, levels = c("up","down","none"))
  
  # num of up down gene
  up_num <- table(df_sub$group) %>% as.data.frame() %>% dplyr::filter(Var1 == "up") %>% dplyr::pull(Freq)
  down_num <- table(df_sub$group) %>% as.data.frame() %>% dplyr::filter(Var1 == "down") %>% dplyr::pull(Freq)
  
  # ggplot2
  # x_lim <- c(min(df_sub$logFC)*1.1, max(df_sub$logFC)*1.1)
  x_lim <- c(-1.5, 1.5)
  p <- ggplot(data = df_sub, aes(x = logFC, y = -log10(pvalue), color = group)) +
    geom_point(size = 1.2) +
    scale_colour_manual(name = "", values = alpha(mycol, 0.7)) +
    scale_x_continuous(limits = x_lim,
                       breaks = c(-1.5, seq(-1.5, 1.5, by = 0.5)),
                       labels = c("-1.5", as.character(seq(-1.5, 1.5, by = 0.5)))) +
    scale_y_continuous(expand = expansion(add = c(0.15, 0)),
                       limits = c(0, 5),
                       breaks = c(0, seq(0, 5, by = 1)),
                       labels = c("0", as.character(seq(0, 5, by = 1)))) +
    geom_hline(yintercept = -log10(cut_pvalue), linewidth = 0.5, color = "black", lty = "dashed") +
    geom_vline(xintercept = c(-cut_log2FC, cut_log2FC), size = 0.5, color = "black", lty = "dashed") +
    # annotate("text",x= -0.7,y= 4,label=down_num,size=6) +
    # annotate("text",x= 0.7,y= 4,label=up_num,size=6) +
    mytheme
  
  print(p)
  
  # save plot
  ggsave(paste0("./7.plot/", "volcano_",out_file, ".pdf"), plot = p, width = 3, height = 2)
}



