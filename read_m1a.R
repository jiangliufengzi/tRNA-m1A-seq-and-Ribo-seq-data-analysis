library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(viridis)

read_m1a <- function(file="./0_tRNA/all_tRNA_mutation.xlsx"){
  m1a <- readxl::read_excel(file, sheet = 6)
  m1a$seqname <- m1a$codon
  m1a <- m1a %>% tidyr::separate(codon, into = c("animo", "codon"), sep = "_")
  m1a$codon <- sapply(m1a$codon, reverse_complement)
  m1a$seqname <- paste(m1a$animo, m1a$codon, sep = "_")
  m1a <- m1a %>% dplyr::filter(animo != "Stop")
}