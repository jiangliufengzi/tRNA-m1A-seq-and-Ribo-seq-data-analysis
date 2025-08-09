from Bio.Data import CodonTable
from collections import defaultdict
from pyfaidx import Fasta
import glob
import os
import pandas as pd

############################################
# 获取标准遗传密码表
############################################
def get_codon_dict():
  codon_table = CodonTable.unambiguous_dna_by_id[1].forward_table
  codon_table.update({"TAA": "Stop", "TAG": "Stop", "TGA": "Stop"})
  
  def init_codon_table():
      return {key: 0 for key in codon_table.keys()}
  
  # 初始化 codon_init 为 defaultdict
  codon_init = defaultdict(float, init_codon_table())
  
  return codon_init

############################################
# 读取覆盖度文件，生成{"chrom|pos": coverage}字典
############################################
def get_coverage(bam_cov_file):
    bam_cov = defaultdict(float)
    
    try:
        with open(bam_cov_file, 'r') as cov_file:
            for line_number, line in enumerate(cov_file, start=1):  # 添加行号追踪
                try:
                    fields = line.strip().split('\t')
                    if len(fields) < 16:  # 检查列数是否足够
                        raise ValueError(f"Line {line_number} has insufficient columns: {line.strip()}")

                    key = f'{fields[0]}|{fields[1]}'
                    value = float(fields[15])  # 尝试转换为浮点数
                    bam_cov[key] = value
                except ValueError as e:
                    print(f"Error parsing line {line_number}: {e}")  # 打印错误信息
                except IndexError:
                    print(f"Error: Missing fields on line {line_number}: {line.strip()}")
    except FileNotFoundError:
        print(f"Error: File '{bam_cov_file}' not found.")
    except IOError as e:
        print(f"Error reading file '{bam_cov_file}': {e}")
    
    return bam_cov

############################################
# 定义反向互补函数
############################################
def reverse_complement(seq):
  rev = seq[::-1]
  complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
  seq = ''.join(complement.get(base, base) for base in rev)
  return seq

############################################
# 读取CDS信息文件-计算密码子覆盖度
############################################
def calculate_codon_coverage(bam_cov_file,
                             gene_list,
                             fasta_file='GRCz11.fasta',
                             cds_info_file='GRCz11_longest_cds_transcripts.csv',
                             peak_high=True
                            ):
  fasta = Fasta(fasta_file)
  codon_init = get_codon_dict()
  bam_cov = get_coverage(bam_cov_file)
  prefix = os.path.basename(bam_cov_file).split(".")[0]
  
  with open(cds_info_file, 'r') as cds_info:
    next(cds_info)
    
    for line in cds_info:
      # 初始化字典
      cds_base = {}
      
      gene_name,cds_count,_,cds_starts,cds_ends,_,strand,chrom,sequence = line.strip().split(',')

      if gene_list == None or (gene_name in gene_list):
      
        if int(cds_count) == 0:
          print(f"The transcript is not codon gene")
          continue
      
        if int(cds_count) > 0:
          cds_starts = sorted(map(int, [s.strip('"\'') for s in cds_starts.split('_')]))
          cds_ends = sorted(map(int, [s.strip('"\'') for s in cds_ends.split('_')]))
          
          # 生成cds的{"chrom|pos": base}字典
          for start, end in zip(cds_starts, cds_ends):
            for pos in range(start, end + 1):
              key = f'{chrom}|{pos}'
              base = str(fasta[chrom][pos - 1]).upper()
              cds_base[key] = base
          
          # 对 cds_base 的键进行排序
          sorted_cds_base = sorted(cds_base.items(), key=lambda x: int(x[0].split('|')[1]))
          
          # 每三个碱基为一个密码子
          for i in range(0, len(sorted_cds_base), 3):
            if len(sorted_cds_base[i:i+3]) < 3:
              continue
              
            codon = ''.join([x[1] for x in sorted_cds_base[i:i+3]])
            if strand == '-':
              codon = str(reverse_complement(codon))

            # 是否加权计算,不加劝的情况下，有覆盖度就加1
            if peak_high:
              coverage = sum([bam_cov[x[0]] for x in sorted_cds_base[i:i+3]]) / 3
              codon_init[codon] += coverage
            else:
              coverage = sum([bam_cov[x[0]] for x in sorted_cds_base[i:i+3]]) / 3
              if coverage > 0:
                 codon_init[codon] += 1
          
  return codon_init, prefix

# 筛选标准密码子
def get_codon_table():
    codon_table = CodonTable.unambiguous_dna_by_id[1].forward_table
    codon_table.update({"TAA": "Stop", "TAG": "Stop", "TGA": "Stop"})
    return codon_table

def output_codon(input_dir,
                 gene_list,
                 fasta_file='GRCz11.fasta',
                 cds_info_file='GRCz11_longest_cds_transcripts.csv',
                 out_codon_usage_file='codon_usage.xlsx',
                 peak_high=True):
    # 定义包含所有文件密码子使用的字典
    all_dict = {}

    # 获取标准密码子表
    codon_table = get_codon_table()

    # 批量处理
    file_list = glob.glob(os.path.join(input_dir, '*.cov.txt'))
    if len(file_list) > 0:
      print(f"正在处理{len(file_list)}个文件")

    for file in file_list:
        codon_init, prefix = calculate_codon_coverage(file, 
                                                      gene_list,
                                                      fasta_file=fasta_file,
                                                      cds_info_file=cds_info_file,
                                                      peak_high=peak_high)
        all_dict[prefix] = codon_init

    for key in all_dict.keys():
        print(key)
        # 使用list()创建副本来避免在迭代时修改字典
        for codon in list(all_dict[key].keys()):
            print(codon, all_dict[key][codon])
            if codon not in codon_table.keys():
                all_dict[key].pop(codon)

    df = pd.DataFrame(all_dict)

    # 输出为excel
    df.to_excel(out_codon_usage_file)

    return df
