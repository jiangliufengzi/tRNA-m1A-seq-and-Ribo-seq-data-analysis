reverse_complement <- function(dna_sequence) {
  # 创建一个命名向量来映射碱基到其互补碱基
  complement_map <- c(A = "T", T = "A", C = "G", G = "C")
  
  # 将输入序列转换为字符向量
  dna_chars <- unlist(strsplit(dna_sequence, split = ""))
  
  # 获取互补碱基
  complement_chars <- complement_map[dna_chars]
  
  # 反转序列
  reverse_complement_chars <- rev(complement_chars)
  
  # 将字符向量合并回一个字符串
  reverse_complement_sequence <- paste(reverse_complement_chars, collapse = "")
  
  return(reverse_complement_sequence)
}
