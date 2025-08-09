
"""
转录本分析工具 (transcript_analyzer.py)

此脚本提供从基因组数据中分析基因转录本、提取CDS序列和识别最长转录本的功能。
既可以作为脚本直接运行，也可以作为模块导入到其他Python程序中。

作为脚本使用:
    python transcript_analyzer.py --db GRCz11_gtf.db --fasta GRCz11.fasta

作为模块导入:
    from transcript_analyzer import enhanced_transcript_analysis
    results = enhanced_transcript_analysis(by="cds", fasta_file="GRCz11.fasta", db_file="GRCz11_gtf.db")
"""

import gffutils
import pandas as pd
from pyfaidx import Fasta
import argparse
import sys
import os

# 定义反向互补函数
def reverse_complement(sequence):
    """
    获取DNA序列的反向互补序列
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                 'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    return ''.join(complement.get(base, base) for base in reversed(sequence))

# 提取CDS序列
def get_cds_sequence(fasta, cds_list, strand):
    """
    获取CDS序列并根据链的方向进行处理
    
    参数：
    - fasta: Fasta对象
    - cds_list: CDS特征列表
    - strand: 链的方向
    
    返回：
    - sequence: 拼接后的CDS序列
    """
    sequences = []
    for cds in cds_list:
        # 获取染色体名称
        chrom = cds.seqid
        # 提取序列
        seq = str(fasta[chrom][cds.start-1:cds.end])
        sequences.append(seq)
    
    # 拼接所有CDS序列
    combined_seq = ''.join(sequences)
    
    # 如果是负链，返回反向互补序列
    if strand == '-':
        return reverse_complement(combined_seq)
    return combined_seq

# 获取最长转录本信息
def get_longest_transcript(db, gene_id, by="cds"):
    """
    获取最长转录本
    
    参数：
    - db: gffutils数据库对象
    - gene_id: 基因ID
    - by: 'cds'表示按CDS长度筛选，'exon'表示按exon长度筛选
    
    返回：
    - transcript_id: 最长转录本ID
    - length: 对应长度
    - gene_name: 基因名称
    """
    
    # 检查参数合法性
    if by.lower() not in ['cds', 'exon']:
        raise ValueError("Parameter 'by' must be either 'cds' or 'exon'")
    
    # 获取基因的所有转录本
    gene = db[gene_id]
    
    # 获取基因名称
    gene_name = None
    if 'gene_name' in gene.attributes:
        gene_name = gene.attributes['gene_name'][0]
    
    transcript_lengths = {}
    
    # 遍历每个转录本
    for transcript in db.children(gene, featuretype='transcript'):
        feature_type = 'CDS' if by.lower() == 'cds' else 'exon'
        
        # 计算长度
        length = sum(
            feature.end - feature.start + 1
            for feature in db.children(transcript, featuretype=feature_type)
        )
        
        transcript_lengths[transcript.id] = length
    
    # 如果没有找到任何符合条件的转录本
    if not transcript_lengths:
        return None, 0, gene_name
    
    # 获取最长的转录本
    longest_transcript_id = max(transcript_lengths.items(), key=lambda x: x[1])[0]
    return longest_transcript_id, transcript_lengths[longest_transcript_id], gene_name

# 提取最长CDS信息及序列
def enhanced_transcript_analysis(by="cds", 
                                 fasta_file="GRCz11.fasta", 
                                 db_file="GRCz11_gtf.db", 
                                 out_fasta="GRCz11_longest_cds_transcripts.fasta",
                                 out_info="GRCz11_longest_cds_transcripts.csv"):
    """
    增强版转录本分析
    
    参数：
    - gtf_file: GTF文件路径
    - fasta_file: 基因组fasta文件路径
    - by: 'cds'或'exon'
    """
    # db = gffutils.create_db(gtf_file, ":memory:", merge_strategy='merge')
    db = gffutils.FeatureDB(db_file, keep_order=True)
    
    # 读取fasta文件
    fasta = Fasta(fasta_file)
    
    results = []
    fasta_info = {}
    
    for gene in db.features_of_type('gene'):
        # gene = "ENSDARG00000009657"
        # transcript_id = longest_transcript_id
        # length = transcript_lengths[longest_transcript_id]
        transcript_id, _, gene_name = get_longest_transcript(db, gene.id, by)
        
        if transcript_id:
            transcript = db[transcript_id]
            cds_list = list(db.children(transcript, featuretype='CDS'))
            
            # 跳过没有CDS的转录本
            if len(cds_list) == 0:
                continue
               
            # 按位置排序CDS
            cds_list.sort(key=lambda x: x.start)
               
            cds_info = {
                'transcript_id': gene_name,
                'cds_count': len(cds_list),
                'total_cds_length': sum(cds.end - cds.start + 1 for cds in cds_list),
                'cds_starts': '_'.join(str(cds.start) for cds in cds_list),
                'cds_ends': '_'.join(str(cds.end) for cds in cds_list),
                'cds_lengths': '_'.join(str(cds.end - cds.start + 1) for cds in cds_list),
                'strand': transcript.strand,
                'chromosome': transcript.seqid
            }
            
            # 如果提供了fasta文件，获取CDS序列
            if fasta is not None:
                cds_info['cds_sequence'] = get_cds_sequence(fasta, cds_list, transcript.strand)
            
            results.append(cds_info)

            # 输出fasta信息
            fasta_key = f">{gene_name}"
            value = cds_info['cds_sequence']
            fasta_info[fasta_key] = value

    # 保存fasta文件
    with open(out_fasta, 'w') as out_fasta_file:
        for key, value in fasta_info.items():
            out_fasta_file.write(f"{key}\n{value}\n")
    print(f"CDS序列已保存至 {out_fasta}")

    # 转换为DataFrame
    df = pd.DataFrame(results)
    df.to_csv(out_info, index=False)
    print(f"CDS信息已保存至 {out_info}")
    
    # 关闭FASTA文件
    fasta.close()
    
    return df

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='分析基因转录本并提取CDS序列')
    
    parser.add_argument('--db', default="GRCz11_gtf.db",
                        help='gffutils数据库文件路径 (默认: GRCz11_gtf.db)')
    parser.add_argument('--fasta', default="GRCz11.fasta",
                        help='基因组fasta文件路径 (默认: GRCz11.fasta)')
    parser.add_argument('--by', choices=['cds', 'exon'], default='cds',
                        help='确定最长转录本的方法 (默认: cds)')
    parser.add_argument('--output', default='transcript_analysis_results.csv',
                        help='输出文件路径 (默认: transcript_analysis_results.csv)')
    parser.add_argument('--out_fasta', default='GRCz11_longest_cds_transcripts.fasta',
                        help='输出fasta文件路径 (默认: GRCz11_longest_cds_transcripts.fasta)')
    
    return parser.parse_args()

def main():
    """主函数"""
    args = parse_arguments()
    
    # 检查文件是否存在
    if not os.path.exists(args.db):
        sys.exit(f"错误: 数据库文件 '{args.db}' 未找到")
    
    if not os.path.exists(args.fasta):
        sys.exit(f"错误: Fasta文件 '{args.fasta}' 未找到")
    
    try:
        print(f"正在使用 {args.by.upper()} 方法分析转录本...")
        results_df = enhanced_transcript_analysis(
            by=args.by,
            fasta_file=args.fasta,
            db_file=args.db,
            out_fasta=args.out_fasta
        )
        
        # 保存结果
        results_df.to_csv(args.output, index=False)
        print(f"分析完成。结果已保存至 {args.output}")
        print(f"共处理了 {len(results_df)} 个转录本")
        
    except Exception as e:
        sys.exit(f"分析过程中出错: {e}")

# 当作为脚本运行时执行main函数
if __name__ == "__main__":
    main()
