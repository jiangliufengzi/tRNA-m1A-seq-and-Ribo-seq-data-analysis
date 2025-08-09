
"""
RNA覆盖度计算工具
此模块提供计算BAM/SAM文件中RNA覆盖度或计数的功能
既可作为脚本直接运行，也可作为模块导入使用
"""

print("python rna_coverage.py -i input.bam -o output.txt -t coverage")
print("import rna_coverage\nrna_coverage.calculate_rna_coverage("input.bam", "output.txt", type="counts")")

import pysam
from collections import defaultdict
import os
import logging
import argparse
import sys

# 设置日志
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def is_simple_cigar(read):
    """
    检查读段是否具有简单的CIGAR字符串（仅包含M/=/X/N操作）
    
    参数:
        read: pysam的比对片段对象
        
    返回:
        bool: 如果CIGAR仅包含M(0)/=(7)/X(8)/N(3)则为True，否则为False
    """
    # 只保留包含 M(0)/=(7)/X(8)/N(3) 的CIGAR
    if not read.cigartuples:  # 检查是否有CIGAR
        return False
    
    allowed_ops = {0, 3, 7, 8}  # M, N, =, X
    for op, length in read.cigartuples:
        if op not in allowed_ops:
            return False
    return True


def calculate_rna_coverage(input_file, output_file, type="counts"):
    """
    计算BAM/SAM文件中的RNA覆盖度或计数
    
    参数:
        input_file (str): 输入BAM/SAM文件路径
        output_file (str): 输出文件路径
        type (str): 'counts'或'coverage'，默认为'counts'
    
    返回:
        None: 结果写入输出文件
    
    异常:
        ValueError: 如果type不是'counts'或'coverage'
    """
    logger.info(f"处理文件: {input_file}")
    logger.info(f"分析类型: {type}")
    
    # 定义一个函数，返回初始化的字典
    def init_inner_dict():
        return {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0}

    # 创建嵌套字典
    position_dict = defaultdict(init_inner_dict)
    
    # 打开BAM/SAM文件
    with pysam.AlignmentFile(input_file, "rb") as samfile:
        # 遍历每条比对记录
        read_count = 0
        processed_count = 0
        
        for read in samfile:
            read_count += 1
            # 检查是否成功比对
            if not read.is_unmapped and is_simple_cigar(read):
                processed_count += 1
                # 获取参考序列名和比对位置
                refname = samfile.get_reference_name(read.reference_id)
                align_pos = read.reference_start
                query_seq = read.query_sequence
                
                if type == "coverage":
                    blocks = read.get_blocks()
                    query_pos = 0 # 记录query序列的位置

                    for block_start, block_end in blocks:
                        for pos in range(block_start, block_end):
                            key = f"{refname}|{pos + 1}"
                            if query_pos < len(query_seq):
                                base = query_seq[query_pos]
                                position_dict[key][base] += 1
                            query_pos += 1
                        
                elif type == "counts":
                    # 只记录起始位置的计数
                    key = f"{refname}|{align_pos + 1}"
                    base = query_seq[0]
                    position_dict[key][base] += 1
                    
                else:
                    raise ValueError("type参数必须是'counts'或'coverage'")
                
            if read_count % 1000000 == 0:
                logger.info(f"已处理 {read_count} 条读段...")
    
    logger.info(f"总读段数: {read_count}, 处理的读段数: {processed_count}")
    
    # 计算总覆盖度用于RPM标准化
    total_density = sum(sum(x.values()) for x in position_dict.values())
    logger.info(f"总密度: {total_density}")
    
    # 将结果写入输出文件
    with open(output_file, 'w') as f:
        # 写入表头
        header = "chromosome\tposition\tA\tA_count\tT\tT_count\tC\tC_count\tG\tG_count\tN\tN_count\tsum\ttotal\tnorm\tnormalized_total\n"
        f.write(header)
        
        # 对位置进行排序
        sorted_positions = sorted(position_dict.keys(), key=lambda x: (
            x.split('|')[0],  # 首先按染色体名排序
            int(x.split('|')[1])  # 然后按位置数字排序
        ))
        
        # 遍历排序后的位置
        for pos_key in sorted_positions:
            # 分解位置键
            chrom, pos = pos_key.split('|')
            counts = position_dict[pos_key]
            
            # 计算总和
            total = sum(counts.values())
            normalized_total = (total / total_density) * 1000000 if total_density > 0 else 0
            
            # 格式化输出行
            output_line = (f"{chrom}\t{pos}\t"
                         f"A\t{counts['A']}\t"
                         f"T\t{counts['T']}\t"
                         f"C\t{counts['C']}\t"
                         f"G\t{counts['G']}\t"
                         f"N\t{counts['N']}\t"
                         f"sum\t{total}\t"
                         f"norm\t{normalized_total:.6f}\n")
            
            f.write(output_line)
    
    logger.info(f"结果已写入: {output_file}")


def main():
    """
    主函数，处理命令行参数并执行相应操作
    """
    # 设置参数解析器
    parser = argparse.ArgumentParser(description='计算BAM/SAM文件中的RNA覆盖度或计数')
    
    # 输入选项
    parser.add_argument('-i', '--input', required=True, help='输入BAM/SAM文件')
    
    # 输出选项
    parser.add_argument('-o', '--output', help='输出文件')
    
    # 分析类型
    parser.add_argument('-t', '--type', choices=['counts', 'coverage'], default='counts',
                        help='分析类型: counts(默认)或coverage')
    
    # 详细输出
    parser.add_argument('-v', '--verbose', action='store_true', help='启用详细输出')
    
    # 解析参数
    args = parser.parse_args()
    
    # 设置日志级别
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, 
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # 处理输入
    try:
        # 如果未指定输出文件，根据输入文件名生成
        if not args.output:
            base_name = os.path.basename(args.input)
            args.output = f"{os.path.splitext(base_name)[0]}_{args.type}.txt"
            logger.info(f"未指定输出文件，使用: {args.output}")
        
        calculate_rna_coverage(args.input, args.output, args.type)
        logger.info(f"处理完成。结果已写入: {args.output}")
    
    except Exception as e:
        logger.error(f"错误: {str(e)}")
        return 1
    
    return 0


# 当作为脚本运行时执行main函数
if __name__ == "__main__":
    sys.exit(main())
