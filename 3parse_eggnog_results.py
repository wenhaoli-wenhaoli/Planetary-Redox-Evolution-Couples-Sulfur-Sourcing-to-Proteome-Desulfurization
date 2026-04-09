#!/usr/bin/env python3
"""解析 eggNOG 输出并构建 COG 存在/缺失矩阵"""

import pandas as pd
import re
from collections import defaultdict
from pathlib import Path

def parse_eggnog_to_cog_matrix(emapper_file, mapping_file, output_file):
    """将 eggNOG 注释转换为 COG 二元矩阵"""
    
    # 1. 构建基因到基因组的映射
    gene_to_genome = {}
    with open(mapping_file) as f:
        for line in f:
            if line.strip():
                parts = line.strip().split()
                if len(parts) >= 2:
                    genome, gene = parts[0], parts[1]
                    gene_to_genome[gene] = genome
    
    # 2. 解析注释文件
    counts = defaultdict(lambda: defaultdict(int))
    all_cogs = set()
    
    with open(emapper_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            query_id = parts[0]
            egg_ogs = parts[4]  # eggNOG_OGs 列
            
            # 提取 COG 编号
            cog_match = re.search(r'COG\d+', egg_ogs)
            
            if cog_match and query_id in gene_to_genome:
                cog_id = cog_match.group()
                genome_id = gene_to_genome[query_id]
                counts[genome_id][cog_id] += 1
                all_cogs.add(cog_id)
    
    # 3. 转换为 DataFrame
    df = pd.DataFrame.from_dict(counts, orient='index')
    df = df.reindex(columns=sorted(all_cogs)).fillna(0).astype(int)
    
    # 4. 转换为二元矩阵（存在/缺失）
    binary_df = (df > 0).astype(int)
    
    # 5. 保存结果
    binary_df.to_csv(output_file)
    print(f"输出文件: {output_file}")
    print(f"基因组数: {len(binary_df)}, COG 家族数: {len(all_cogs)}")
    
    return binary_df

if __name__ == "__main__":
    # 配置路径
    EGGNOG_OUT = "data/processed/eggnog_output/all_annotations.tsv"
    MAPPING_FILE = "data/processed/order_ids.tsv"  # 基因组-基因ID映射
    OUTPUT_MATRIX = "data/processed/genome_cog_matrix.csv"
    
    # 执行解析
    cog_matrix = parse_eggnog_to_cog_matrix(EGGNOG_OUT, MAPPING_FILE, OUTPUT_MATRIX)