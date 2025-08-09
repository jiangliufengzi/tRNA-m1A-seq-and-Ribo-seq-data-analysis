from Bio.Data import CodonTable
from collections import defaultdict
import os
from pyfaidx import Fasta
import pandas as pd
import glob

########################################################
# 读取覆盖度文件，生成{"chrom|pos": coverage}字典
########################################################
def get_coverage(bam_cov_file, min_count):
    bam_cov = defaultdict(float)
  
    with open(bam_cov_file, 'r') as cov_file:
        for line in cov_file:
            fields = line.strip().split('\t')
            if len(fields) > 15:
                key = f'{fields[0]}|{fields[1]}'
                value = float(fields[15])
                sum = int(fields[13])
                if sum > min_count:
                    bam_cov[key] = value

    return bam_cov

########################################################
# 读取最长转录本文件
########################################################
def calculate_codon_coverage(bam_cov_file,
                             gene_list=None,
                             fasta_file='GRCz11.fasta',
                             cds_info_file='GRCz11_longest_cds_transcripts.csv',
                             min_count=0,
                             window=100
                            ):
    fasta = Fasta(fasta_file)
    bam_cov = get_coverage(bam_cov_file, min_count)
    prefix = os.path.basename(bam_cov_file).split(".")[0]
    pos_dict_new = {}
  
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
                    sorted_cds_base = dict(sorted(cds_base.items(), key=lambda x: int(x[0].split('|')[1])))
                    sorted_cds_keys = list(sorted_cds_base.keys())
                    max_idx = len(sorted_cds_keys) - 1
            
                    # 计算转录本上每个位置的分数
                    for i in range(0, len(sorted_cds_keys)):

                        # 计算窗口边界，确保不超过列表长度
                        pos_n = min(i + int(window), max_idx)
                        pos_n_mid = min(int(i + int(window)/2), max_idx)
                        pos_n_extend = min(int(i + 1.5*int(window)), max_idx)

                        # 安全的切片操作
                        window1 = sorted_cds_keys[i:pos_n + 1]
                        window2 = sorted_cds_keys[pos_n_mid:pos_n_extend + 1]

                        # 确保窗口非空
                        if not window1 or not window2:
                            continue

                        sum_1 = sum([bam_cov[x] for x in window1])
                        sum_2 = sum([bam_cov[x] for x in window2])

                        if sum_1*sum_2 == 0:
                            si = 0
                        else:
                            si = (int(window)/2)*bam_cov[window1[-1]]*((sum_1 + sum_2)/(sum_1*sum_2))

                        # save in dict
                        pos_dict_new[sorted_cds_keys[i]] = si

    return pos_dict_new

########################################################
# 输出为数据框
########################################################
def get_score(bam_cov_file,
                gene_list=None,
                fasta_file='GRCz11.fasta',
                cds_info_file='GRCz11_longest_cds_transcripts.csv',
                min_count=0,
                window=100):
    
    # file_list = glob.glob(os.path.join(input_dir, '*.cov.txt'))

    score_dict = calculate_codon_coverage(bam_cov_file=bam_cov_file,
                                        gene_list=gene_list,
                                        fasta_file=fasta_file,
                                        cds_info_file=cds_info_file,
                                        min_count=min_count,
                                        window=window
                                        )

    # 创建数据框
    score_keys = list(score_dict.keys())
    data = {
        "chrom": [p.split('|')[0] for p in score_keys],
        "pos": [int(p.split('|')[1]) for p in score_keys],
        "value": [float(score_dict[p]) for p in score_keys]
    }

    df = pd.DataFrame(data)

    return df

