# collect_all_presence.py （通用版：不限制文件夹名以 K 开头）
import os
import sys
import pandas as pd
from collections import defaultdict
from pathlib import Path

def main(input_dir, output_file):
    base_dir = Path(input_dir).resolve()
    if not base_dir.is_dir():
        print(f"错误：输入路径不是目录: {base_dir}")
        sys.exit(1)

    all_nodes = set()
    data = defaultdict(dict)          # node -> {family: presence}

    # 遍历所有子文件夹（不限制名字）
    for item in base_dir.iterdir():
        if not item.is_dir():
            continue
        
        family_name = item.name.strip()  # 直接用文件夹名作为标识（如 K00012、sat 等）
        
        recon_dir = item / "reconciliations"
        if not recon_dir.is_dir():
            print(f"跳过 {family_name}：未找到 reconciliations 目录")
            continue

        file_path = recon_dir / "perspecies_eventcount.txt"
        if not file_path.is_file():
            print(f"跳过 {family_name}：未找到 {file_path}")
            continue

        print(f"正在读取 {family_name} ...")

        with open(file_path, encoding='utf-8') as f:
            header = f.readline().strip()  # 跳过表头
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = [p.strip() for p in line.split(',')]
                if len(parts) < 6:
                    continue

                node = parts[0]
                try:
                    presence = float(parts[5])  # presence 在第6列（索引5）
                except (ValueError, IndexError):
                    presence = 0.0

                data[node][family_name] = presence
                all_nodes.add(node)

    if not all_nodes:
        print("未找到任何有效的节点数据")
        return

    print(f"\n共收集到 {len(all_nodes)} 个节点，{len(set(fam for d in data.values() for fam in d))} 个 family/基因")

    sorted_nodes = sorted(all_nodes)  # 可按需自定义排序
    family_list = sorted({fam for node_dict in data.values() for fam in node_dict})

    rows = []
    for node in sorted_nodes:
        row = {'node': node}
        for fam in family_list:
            row[f'presence_{fam}'] = data[node].get(fam, 0.0)
        rows.append(row)

    df = pd.DataFrame(rows)
    columns = ['node'] + [f'presence_{fam}' for fam in family_list]
    df = df[columns]

    df.to_csv(output_file, sep='\t', index=False, float_format='%.6f')
    print(f"\n已生成文件：{output_file}")
    print(f"总行数（节点数）：{len(df)}")
    print(f"总列数（1 + family数）：{len(df.columns)}")
    print("\n前 5 行预览：")
    print(df.head().to_string(index=False))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法:")
        print("  python collect_all_presence.py  <输入文件夹>  <输出tsv文件名>")
        print("示例:")
        print("  python collect_all_presence.py cir/ all_presence.tsv")
        print("  python collect_all_presence.py output_reconciliation/cir/ presence_matrix.tsv")
        sys.exit(1)

    input_folder = sys.argv[1]
    output_tsv = sys.argv[2]
    main(input_folder, output_tsv)