
import pandas as pd
from Bio.Data import CodonTable
from Bio import SeqIO
from collections import defaultdict

print("Usage: from gene_codon_usage import get_codon_use\nget_codon_use(gene_list=None, input_fasta='GRCz11_longest_cds_transcripts.fasta')")

# 获取标准遗传密码表
codon_table = CodonTable.unambiguous_dna_by_id[1].forward_table
codon_table.update({"TAA": "Stop", "TAG": "Stop", "TGA": "Stop"})

# 定义输出函数
def get_codon_use(gene_list=None, input_fasta="GRCz11_longest_cds_transcripts.fasta"):
    """
    计算指定基因列表或整个fasta文件的密码子使用频率
    
    参数:
    - gene_list: 基因ID列表，如果为None则分析整个fasta文件
    - input_fasta: fasta文件路径
    
    返回:
    - DataFrame: 包含密码子、氨基酸和计数的数据框
    """
    # 初始化密码子计数器
    codon_count = defaultdict(int)
    
    # 计算fasta的密码子使用
    for record in SeqIO.parse(input_fasta, "fasta"):
        sequence = str(record.seq).upper()
        id = str(record.id)
        
        # 如果gene_list为None或者当前基因在gene_list中
        if gene_list is None or id in gene_list:
            for i in range(0, len(sequence) - 2, 3):
                print(i)
                codon = sequence[i:i+3]
                if len(codon) == 3:
                    codon_count[codon] += 1
  
    # 打印出密码子使用情况         
    for codon, count in codon_count.items():
        amino_acid = codon_table.get(codon, "Unknown")
        print(f"{codon} ({amino_acid}): {count}")
    
    # 转换为数据框输出  
    data = {
        "codon": list(codon_count.keys()),
        "Amino Acid": [codon_table.get(codon, "Unknown") for codon in codon_count.keys()],
        "count": list(codon_count.values())
    }
    
    df = pd.DataFrame(data)
    return df


