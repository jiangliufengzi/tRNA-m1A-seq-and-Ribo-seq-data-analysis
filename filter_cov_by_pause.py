##############################################################
# 加载包
##############################################################
import pandas as pd
from collections import defaultdict
import os

##############################################################
# 解析pause位点文件，生成字典{"chrom|pos": value}
##############################################################
def parse_pause_site(pause_site_file):
    pause_site_dict = {}
    with open(pause_site_file, "r") as infile:
        next(infile)
        for line in infile:
            chrom, pos, value = line.strip().split("\t")
            pause_key = f"{chrom}|{pos}"
            pause_value = float(value)

            pause_site_dict[pause_key] = pause_value
    return pause_site_dict

##############################################################
# 解析cov文件，生成字典{"chrom|pos": norm}
##############################################################
def get_coverage(bam_cov_file):
  bam_cov = defaultdict(float)
  
  with open(bam_cov_file, 'r') as cov_file:
    for line in cov_file:
      fields = line.strip().split('\t')
      key = f'{fields[0]}|{fields[1]}'
      value = float(fields[15])
      bam_cov[key] = value
  
  return bam_cov

##############################################################
# 筛选位于pause位点的cov信息
##############################################################
def get_dict_key_intersection(cov_dict, pause_dict):
    # 将两个字典的键转换为集合并计算交集
    common_keys = set(cov_dict.keys()) & set(pause_dict.keys())
    
    cov_new = {k: cov_dict[k] for k in common_keys}

    cov_keys = cov_new.keys()
    data = {
       "chrom": [str(x.split("|")[0]) for x in cov_keys],
       "pos": [int(x.split("|")[1]) for x in cov_keys],
       "value": [float(cov_new[x]) for x in cov_keys]
    }

    df=pd.DataFrame(data)

    # 填充为16列
    df_filled = df.copy()  # 复制数据框
    df_filled.insert(2, '_1', '_')  # 插入第3列
    df_filled.insert(3, '_2', '_')  # 插入第4列
    df_filled.insert(4, '_3', '_')  # 插入第5列
    df_filled.insert(5, '_4', '_')  # 插入第6列
    df_filled.insert(6, '_5', '_')  # 插入第7列
    df_filled.insert(7, '_6', '_')  # 插入第8列
    df_filled.insert(8, '_7', '_')  # 插入第9列
    df_filled.insert(9, '_8', '_')  # 插入第10列
    df_filled.insert(10, '_9', '_')  # 插入第11列
    df_filled.insert(11, '_10', '_')  # 插入第12列
    df_filled.insert(12, '_11', '_')  # 插入第13列
    df_filled.insert(13, '_12', '_')  # 插入第14列
    df_filled.insert(14, '_13', '_')  # 插入第15列

    # 将原来的 value 列（第3列）移动到第16列
    df_filled['value_16'] = df_filled.pop('value')

    return df_filled

##############################################################
# 输出筛选后的cov文件为txt格式
##############################################################
def get_pause_site_cov(bam_cov_file,
                       pause_site_file,
                       output_file=None):
    cov_dict = get_coverage(bam_cov_file)
    pause_dict = parse_pause_site(pause_site_file)

    df = get_dict_key_intersection(cov_dict, pause_dict)

    # 如果没有指定输出文件，则自动生成
    if output_file is None:
        base_name = os.path.basename(bam_cov_file).split('.')[0]
        output_file = f"{base_name}_pause_sites.txt"
    
    # 将结果输出到文件
    df.to_csv(output_file, sep='\t', index=False, header=False)
    
    print(f"Found {len(df)} positions at pause sites")
    print(f"Results saved to {output_file}")
    
    return df, output_file
   
   

