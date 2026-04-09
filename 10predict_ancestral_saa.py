#!/usr/bin/env python3
"""
predict_ancestral_saa.py

Predict ancestral SAA frequencies using trained XGBoost model
Automatically runs both separate and merged predictions

Usage:
    python predict_ancestral_saa.py

Input files:
    cog_cir.tsv      - COG presence probabilities under CIR molecular clock model
    cog_ln.tsv       - COG presence probabilities under LN molecular clock model
    cog_ugam.tsv     - COG presence probabilities under UGAM molecular clock model

Input file format:
    - Each row represents a node in the phylogenetic tree
    - First column: node ID (digits = internal nodes, letters = tips/outgroup)
    - Other columns: presence_CXXXXX = probability (0-1) of COG XXXX being present
    - Columns starting with 'presence_C' or 'C' are COG features

Output files:
    - Separate files for each model (for supplementary materials)
    - Merged file (for main text figure)

输入文件:
    cog_cir.tsv      - CIR分子钟模型下，各节点的COG存在概率
    cog_ln.tsv       - LN分子钟模型下，各节点的COG存在概率
    cog_ugam.tsv     - UGAM分子钟模型下，各节点的COG存在概率

输入文件格式:
    - 每行代表系统发育树中的一个节点
    - 第一列: 节点ID (数字=内部节点, 字母=末端节点/外群)
    - 其他列: presence_CXXXXX = 该COG存在的概率 (0-1)
    - 以 'presence_C' 或 'C' 开头的列为COG特征

输出文件:
    - 每个模型分开的文件（用于附件）
    - 合并的文件（用于正文主图）
"""

import pandas as pd
import numpy as np
import joblib
import os
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

print("=" * 70)
print("Ancestral SAA Frequency Prediction")
print("=" * 70)

# Load model
print("\nLoading model...")
selected_features = joblib.load('results/selected_features.pkl')
scaler = joblib.load('results/scaler.pkl')
best_model = joblib.load('results/best_model.pkl')
print("Model loaded successfully")

# Input files
input_files = ['cog_cir.tsv', 'cog_ln.tsv', 'cog_ugam.tsv']
n_samples_per_row = 100
output_dir = 'results'
os.makedirs(output_dir, exist_ok=True)

# ==========================================
# PART 1: Separate processing (for supplementary materials)
# ==========================================
print("\n" + "#" * 70)
print("# PART 1: Separate processing (for supplementary materials)")
print("#" * 70)

for file in input_files:
    if not os.path.exists(file):
        print(f"\nSkipping {file} (file not found)")
        continue
    
    print(f"\n{'=' * 60}")
    print(f"Processing {file}")
    print(f"{'=' * 60}")
    
    # Read file
    df = pd.read_csv(file, sep='\t', low_memory=False)
    print(f"Raw rows: {len(df)}")
    
    node_col = 'node' if 'node' in df.columns else df.columns[0]
    
    # Filter internal nodes (starting with digits)
    df_filtered = df[df[node_col].astype(str).str[0].str.isdigit()].copy()
    print(f"Internal nodes: {len(df_filtered)}")
    
    # Data augmentation
    print(f"Generating {n_samples_per_row} augmented samples per row...")
    cog_columns = [col for col in df_filtered.columns 
                   if col.startswith('presence_C') or col.startswith('C')]
    
    augmented_rows = []
    for idx, row in tqdm(df_filtered.iterrows(), total=len(df_filtered), desc="Augmenting"):
        augmented_rows.append(row.to_dict())
        
        original_probs = row[cog_columns].values.astype(float)
        original_probs = np.nan_to_num(original_probs, nan=0.0)
        
        random_binary = (np.random.rand(n_samples_per_row, len(cog_columns)) < original_probs)
        random_df = pd.DataFrame(random_binary, columns=cog_columns)
        
        for col in df_filtered.columns:
            if col not in cog_columns:
                random_df[col] = row[col]
        
        augmented_rows.extend(random_df.to_dict('records'))
    
    df_aug = pd.DataFrame(augmented_rows)
    print(f"Total rows after augmentation: {len(df_aug)}")
    
    # Rename columns
    rename_dict = {col: col.replace('presence_', '') 
                   for col in df_aug.columns if col.startswith('presence_')}
    df_aug = df_aug.rename(columns=rename_dict)
    
    # Predict
    X_predict = df_aug.reindex(columns=selected_features, fill_value=0.0)
    X_predict_scaled = scaler.transform(X_predict)
    df_aug['predicted_s_aa_freq'] = best_model.predict(X_predict_scaled)
    
    # Save
    prefix = os.path.basename(file).replace('.tsv', '')
    df_aug.to_csv(f"{output_dir}/predicted_{prefix}_all.tsv", sep='\t', index=False)
    
    original_indices = range(len(df_filtered))
    predicted_original = df_aug.iloc[original_indices][[node_col, 'predicted_s_aa_freq']].copy()
    predicted_original.to_csv(f"{output_dir}/predicted_{prefix}_original_only.tsv", sep='\t', index=False)
    
    print(f"Saved: predicted_{prefix}_original_only.tsv ({len(predicted_original)} nodes)")

# ==========================================
# PART 2: Merged processing (for main text figure)
# ==========================================
print("\n" + "#" * 70)
print("# PART 2: Merged processing (for main text figure)")
print("#" * 70)

existing_files = [f for f in input_files if os.path.exists(f)]

if existing_files:
    print(f"\nMerging {len(existing_files)} files...")
    
    # Load and merge all files
    dfs = []
    for file in existing_files:
        print(f"  Reading {file}")
        dfs.append(pd.read_csv(file, sep='\t', low_memory=False))
    
    x_merged = pd.concat(dfs, ignore_index=True)
    print(f"Total merged rows: {len(x_merged)}")
    
    # Filter internal nodes
    node_col = 'node' if 'node' in x_merged.columns else x_merged.columns[0]
    x_filtered = x_merged[x_merged[node_col].astype(str).str[0].str.isdigit()].copy()
    print(f"Internal nodes: {len(x_filtered)}")
    
    # Data augmentation
    print(f"Generating {n_samples_per_row} augmented samples per row...")
    cog_columns = [col for col in x_filtered.columns 
                   if col.startswith('presence_C') or col.startswith('C')]
    
    augmented_rows = []
    for idx, row in tqdm(x_filtered.iterrows(), total=len(x_filtered), desc="Augmenting"):
        augmented_rows.append(row.to_dict())
        
        original_probs = row[cog_columns].values.astype(float)
        original_probs = np.nan_to_num(original_probs, nan=0.0)
        
        random_binary = (np.random.rand(n_samples_per_row, len(cog_columns)) < original_probs)
        random_df = pd.DataFrame(random_binary, columns=cog_columns)
        
        for col in x_filtered.columns:
            if col not in cog_columns:
                random_df[col] = row[col]
        
        augmented_rows.extend(random_df.to_dict('records'))
    
    x_augmented = pd.DataFrame(augmented_rows)
    print(f"Total rows after augmentation: {len(x_augmented)}")
    
    # Rename columns
    rename_dict = {col: col.replace('presence_', '') 
                   for col in x_augmented.columns if col.startswith('presence_')}
    x_augmented = x_augmented.rename(columns=rename_dict)
    
    # Predict
    X_predict = x_augmented.reindex(columns=selected_features, fill_value=0.0)
    X_predict_scaled = scaler.transform(X_predict)
    x_augmented['predicted_s_aa_freq'] = best_model.predict(X_predict_scaled)
    
    # Save all results
    x_augmented.to_csv(f"{output_dir}/predicted_merged_all_rows.tsv", sep='\t', index=False)
    
    # Save original nodes only (for main text figure)
    original_indices = range(len(x_filtered))
    predicted_original = x_augmented.iloc[original_indices][[node_col, 'predicted_s_aa_freq']].copy()
    predicted_original.to_csv(f"{output_dir}/predicted_merged_original_only.tsv", sep='\t', index=False)
    
    print(f"\nSaved: predicted_merged_original_only.tsv ({len(predicted_original)} nodes)")
    print("\nPreview (first 10 rows):")
    print(predicted_original.head(10))
else:
    print("Error: No input files found")

# Final summary
print("\n" + "=" * 70)
print("PREDICTION COMPLETED")
print("=" * 70)
print("\nOutput files:")
print("  [Supplementary Materials]")
for file in input_files:
    prefix = os.path.basename(file).replace('.tsv', '') if os.path.exists(file) else file.replace('.tsv', '')
    if os.path.exists(file):
        print(f"    - results/predicted_{prefix}_original_only.tsv")
print("  [Main Text Figure]")
print("    - results/predicted_merged_original_only.tsv")
print("=" * 70)