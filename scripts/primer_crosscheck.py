#!/usr/bin/env python3
"""
GCAプライマー クロスチェックスクリプト
設計したプライマーの特異性を近縁種・動植物配列と照合して検証

機能:
- プライマーと非標的配列のアライメント
- ミスマッチ数の計算（全体および3'末端）
- 交差反応リスクの評価
- 詳細レポート生成
"""

import os
import sys
import re
import csv
import argparse
from datetime import datetime
from pathlib import Path
from Bio import SeqIO

# 定数（このスクリプトの親の親ディレクトリを基準）
SCRIPT_DIR = Path(__file__).parent.absolute()
PROJECT_ROOT = SCRIPT_DIR.parent

# デフォルトのディレクトリ
DEFAULT_SEQUENCES_DIR = PROJECT_ROOT / "sequences"
DEFAULT_PRIMERS_DIR = PROJECT_ROOT / "primer_candidates"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "crosscheck_results"

# グローバル変数（mainで上書き可能）
BASE_DIR = PROJECT_ROOT
SEQUENCES_DIR = DEFAULT_SEQUENCES_DIR
PRIMERS_DIR = DEFAULT_PRIMERS_DIR
OUTPUT_DIR = DEFAULT_OUTPUT_DIR

# ========================================
# 動的データ読み込み関数
# ========================================

def load_reference_sequences():
    """FASTAファイルから参照配列を読み込む"""
    references = {}
    
    # カテゴリとファイル名のマッピング
    files_to_load = [
        ("target", os.path.join(SEQUENCES_DIR, "target_sequence.fasta")),
        ("related", os.path.join(SEQUENCES_DIR, "related_sequences.fasta")),
        ("nontarget", os.path.join(SEQUENCES_DIR, "nontarget_sequences.fasta"))
    ]
    
    for category, filepath in files_to_load:
        if os.path.exists(filepath):
            try:
                for record in SeqIO.parse(filepath, "fasta"):
                    # IDから種名らしき部分を抽出（簡易的な処理。必要に応じて拡張）
                    # "NC_001620.1 Artemia franciscana mitochondrion..." などの形式を想定
                    parts = record.description.split(" ", 2)
                    if len(parts) >= 3:
                        org_name = f"{parts[1]}_{parts[2].split()[0]}"
                    else:
                        org_name = record.id.replace(".", "_")
                        
                    references[org_name] = {
                        "accession": record.id,
                        "category": category,
                        "sequence": str(record.seq)
                    }
                print(f"Loaded {category} sequences from {filepath}")
            except Exception as e:
                print(f"Error loading {filepath}: {e}")
        else:
            print(f"Warning: Sequence file not found: {filepath}")
            
    return references

def load_primer_candidates():
    """CSVファイルからプライマー候補を読み込む"""
    primers = {}
    primer_csv = os.path.join(PRIMERS_DIR, "primer_candidates.csv")
    
    if os.path.exists(primer_csv):
        try:
            with open(primer_csv, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    primer_id = row.get("ID", row.get("id", "Unknown"))
                    if primer_id == "Unknown":
                        continue
                    
                    primers[primer_id] = {
                        "forward": row.get("Forward", row.get("Forward_Sequence", "")),
                        "reverse": row.get("Reverse", row.get("Reverse_Sequence", "")),
                        "product_size": int(row.get("Product_Size", row.get("product_size", "0")))
                    }
            print(f"Loaded {len(primers)} primer candidates from {primer_csv}")
        except Exception as e:
            print(f"Error loading {primer_csv}: {e}")
    else:
        print(f"Warning: Primer candidates file not found: {primer_csv}")
        
    return primers

# ========================================
# 解析ロジック
# ========================================

def reverse_complement(seq):
    """逆相補鎖を返す"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return ''.join([complement.get(base, 'N') for base in seq[::-1]])

def find_primer_binding_sites(primer_seq, target_seq, max_mismatches=5):
    """
    プライマー結合部位を検索
    ミスマッチ許容でアライメント
    """
    primer_seq = primer_seq.upper()
    target_seq = target_seq.upper()
    primer_len = len(primer_seq)

    binding_sites = []

    # Forward方向
    for i in range(len(target_seq) - primer_len + 1):
        segment = target_seq[i:i + primer_len]
        mismatches = sum(1 for a, b in zip(primer_seq, segment) if a != b)

        if mismatches <= max_mismatches:
            # 3'末端のミスマッチをカウント
            three_prime_mismatches = sum(
                1 for a, b in zip(primer_seq[-5:], segment[-5:]) if a != b
            )
            binding_sites.append({
                'position': i,
                'direction': 'forward',
                'mismatches': mismatches,
                'three_prime_mismatches': three_prime_mismatches,
                'aligned_seq': segment
            })

    # Reverse complement方向
    primer_rc = reverse_complement(primer_seq)
    for i in range(len(target_seq) - primer_len + 1):
        segment = target_seq[i:i + primer_len]
        mismatches = sum(1 for a, b in zip(primer_rc, segment) if a != b)

        if mismatches <= max_mismatches:
            three_prime_mismatches = sum(
                1 for a, b in zip(primer_rc[-5:], segment[-5:]) if a != b
            )
            binding_sites.append({
                'position': i,
                'direction': 'reverse_complement',
                'mismatches': mismatches,
                'three_prime_mismatches': three_prime_mismatches,
                'aligned_seq': segment
            })

    return binding_sites


def calculate_crossreaction_risk(forward_sites, reverse_sites):
    """
    交差反応リスクを計算
    Forward と Reverse の両方が結合する場合のリスクを評価
    """
    if not forward_sites or not reverse_sites:
        return 0.0, "極めて低い", []

    # 最良の結合サイト（最小ミスマッチ）
    best_forward = min(forward_sites, key=lambda x: x['mismatches'])
    best_reverse = min(reverse_sites, key=lambda x: x['mismatches'])

    total_mismatches = best_forward['mismatches'] + best_reverse['mismatches']
    three_prime_mismatches = (best_forward['three_prime_mismatches'] +
                              best_reverse['three_prime_mismatches'])

    # リスクスコア計算
    # 3'末端ミスマッチは重要度が高い（係数2.0）
    risk_score = total_mismatches + (three_prime_mismatches * 1.5)

    # リスクレベル判定
    if risk_score >= 8:
        risk_level = "極めて低い (<1%)"
        risk_percent = 0.5
    elif risk_score >= 6:
        risk_level = "低い (1-5%)"
        risk_percent = 3.0
    elif risk_score >= 4:
        risk_level = "中程度 (5-20%)"
        risk_percent = 12.0
    elif risk_score >= 2:
        risk_level = "高い (20-50%)"
        risk_percent = 35.0
    else:
        risk_level = "極めて高い (>50%)"
        risk_percent = 70.0

    details = {
        'forward_mismatches': best_forward['mismatches'],
        'reverse_mismatches': best_reverse['mismatches'],
        'forward_3prime': best_forward['three_prime_mismatches'],
        'reverse_3prime': best_reverse['three_prime_mismatches'],
        'total_mismatches': total_mismatches,
        'risk_score': risk_score
    }

    return risk_percent, risk_level, details


def crosscheck_primer_set(primer_name, primer_data, reference_sequences):
    """
    1つのプライマーセットを全参照配列とクロスチェック
    """
    results = {
        'primer_name': primer_name,
        'forward': primer_data['forward'],
        'reverse': primer_data['reverse'],
        'product_size': primer_data['product_size'],
        'organisms': {}
    }

    for org_name, org_data in reference_sequences.items():
        # Forward プライマーの結合サイト検索
        forward_sites = find_primer_binding_sites(
            primer_data['forward'],
            org_data['sequence']
        )

        # Reverse プライマーの結合サイト検索（逆相補で検索）
        reverse_sites = find_primer_binding_sites(
            primer_data['reverse'],
            org_data['sequence']
        )

        # リスク計算
        risk_percent, risk_level, details = calculate_crossreaction_risk(
            forward_sites, reverse_sites
        )

        results['organisms'][org_name] = {
            'category': org_data['category'],
            'accession': org_data['accession'],
            'forward_sites': len(forward_sites),
            'reverse_sites': len(reverse_sites),
            'best_forward_mismatch': min([s['mismatches'] for s in forward_sites]) if forward_sites else 'N/A',
            'best_reverse_mismatch': min([s['mismatches'] for s in reverse_sites]) if reverse_sites else 'N/A',
            'risk_percent': risk_percent,
            'risk_level': risk_level,
            'details': details
        }

    return results


def generate_crosscheck_report(all_results):
    """
    クロスチェック結果のレポートを生成
    """
    report = []
    report.append("=" * 80)
    report.append("GCAプライマー クロスチェックレポート")
    report.append(f"生成日時: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("=" * 80)

    # カテゴリ別の日本語名
    category_names = {
        'target': '主標的',
        'related': '近縁種',
        'nontarget': '非標的生物',
        'artemia_other': 'Artemia他種',
        'crustacea': '甲殻類',
        'human': 'ヒト',
        'insect': '昆虫類',
        'plant': '植物',
        'nematode': '線虫'
    }

    for result in all_results:
        report.append("")
        report.append("=" * 80)
        report.append(f"プライマーセット: {result['primer_name']}")
        report.append("=" * 80)
        report.append(f"Forward: 5'-{result['forward']}-3' ({len(result['forward'])} bp)")
        report.append(f"Reverse: 5'-{result['reverse']}-3' ({len(result['reverse'])} bp)")
        report.append(f"産物サイズ: {result['product_size']} bp")
        report.append("")

        # カテゴリ別にまとめる
        categories = {}
        for org_name, org_result in result['organisms'].items():
            cat = org_result['category']
            if cat not in categories:
                categories[cat] = []
            categories[cat].append((org_name, org_result))

        report.append("-" * 80)
        report.append("交差反応性評価サマリー")
        report.append("-" * 80)
        report.append("")
        report.append(f"{'生物群':<20} {'代表種':<30} {'ミスマッチ(F/R)':<15} {'リスク':<20}")
        report.append("-" * 80)

        for cat in ['target', 'related', 'nontarget', 'artemia_other', 'crustacea', 'insect', 'plant', 'human', 'nematode']:
            if cat in categories:
                for org_name, org_result in categories[cat]:
                    f_mm = org_result['best_forward_mismatch']
                    r_mm = org_result['best_reverse_mismatch']
                    mm_str = f"{f_mm}/{r_mm}"

                    # リスクレベルに応じたマーク
                    risk = org_result['risk_percent']
                    if cat == 'target':
                        mark = "✓ 標的"
                    elif risk < 1:
                        mark = "✓✓ 極めて低い"
                    elif risk < 5:
                        mark = "✓ 低い"
                    elif risk < 20:
                        mark = "⚠ 中程度"
                    else:
                        mark = "✗ 高い"

                    report.append(f"{category_names.get(cat, cat):<20} {org_name:<30} {mm_str:<15} {mark:<20}")

        report.append("")
        report.append("-" * 80)
        report.append("詳細解析")
        report.append("-" * 80)

        for cat in ['related', 'nontarget', 'artemia_other', 'crustacea', 'insect', 'plant', 'human', 'nematode']:
            if cat in categories:
                report.append(f"\n【{category_names.get(cat, cat)}】")
                for org_name, org_result in categories[cat]:
                    details = org_result.get('details', {})
                    report.append(f"  {org_name} ({org_result['accession']})")
                    report.append(f"    Forward ミスマッチ: {org_result['best_forward_mismatch']}塩基")
                    report.append(f"    Reverse ミスマッチ: {org_result['best_reverse_mismatch']}塩基")
                    if details:
                        report.append(f"    3'末端ミスマッチ(F): {details.get('forward_3prime', 'N/A')}塩基")
                        report.append(f"    3'末端ミスマッチ(R): {details.get('reverse_3prime', 'N/A')}塩基")
                    report.append(f"    交差反応リスク: {org_result['risk_level']}")

    # 総合評価
    report.append("")
    report.append("=" * 80)
    report.append("総合評価・推奨事項")
    report.append("=" * 80)
    report.append("")

    # 各プライマーセットの総合スコア計算
    primer_scores = []
    for result in all_results:
        total_risk = 0
        high_risk_count = 0
        for org_name, org_result in result['organisms'].items():
            if org_result['category'] != 'target':
                total_risk += org_result['risk_percent']
                if org_result['risk_percent'] > 5:
                    high_risk_count += 1

        avg_risk = total_risk / (len(result['organisms']) - 1)  # 標的以外
        primer_scores.append({
            'name': result['primer_name'],
            'avg_risk': avg_risk,
            'high_risk_count': high_risk_count
        })

    # ソート（リスクが低い順）
    primer_scores.sort(key=lambda x: (x['high_risk_count'], x['avg_risk']))

    report.append("プライマー推奨順位:")
    for i, ps in enumerate(primer_scores, 1):
        if ps['high_risk_count'] == 0:
            rec = "強く推奨"
        elif ps['high_risk_count'] <= 1:
            rec = "推奨"
        else:
            rec = "要検討"
        report.append(f"  {i}. {ps['name']}: 平均リスク {ps['avg_risk']:.1f}%, "
                     f"高リスク生物 {ps['high_risk_count']}種 → {rec}")

    report.append("")
    report.append("-" * 80)
    report.append("注意事項:")
    report.append("- 本解析はIn-silico予測であり、実験的検証が必要です")
    report.append("- 3'末端のミスマッチは特に重要（PCR阻害効果大）")
    report.append("- 甲殻類（Daphnia, Penaeus）との交差反応に特に注意")
    report.append("- NCBI Primer-BLASTでの追加検証を推奨")
    report.append("=" * 80)

    return "\n".join(report)


def generate_csv_report(all_results, output_file):
    """CSV形式のレポートを生成"""
    with open(output_file, 'w', encoding='utf-8') as f:
        # ヘッダー
        f.write("Primer_Set,Organism,Category,Accession,Forward_Mismatch,Reverse_Mismatch,"
               "3prime_F_Mismatch,3prime_R_Mismatch,Risk_Percent,Risk_Level\n")

        for result in all_results:
            for org_name, org_result in result['organisms'].items():
                details = org_result.get('details', {})
                # detailsがリストの場合は空の辞書として扱う
                if isinstance(details, list):
                    details = {}
                f.write(f"{result['primer_name']},{org_name},{org_result['category']},"
                       f"{org_result['accession']},{org_result['best_forward_mismatch']},"
                       f"{org_result['best_reverse_mismatch']},"
                       f"{details.get('forward_3prime', 'N/A') if isinstance(details, dict) else 'N/A'},"
                       f"{details.get('reverse_3prime', 'N/A') if isinstance(details, dict) else 'N/A'},"
                       f"{org_result['risk_percent']},{org_result['risk_level']}\n")


def main():
    global SEQUENCES_DIR, PRIMERS_DIR, OUTPUT_DIR

    # コマンドライン引数パーサー
    parser = argparse.ArgumentParser(description="GCAプライマー クロスチェック解析")
    parser.add_argument("--sequences-dir", type=str, default=None,
                       help="配列ファイルのディレクトリ")
    parser.add_argument("--primers-dir", type=str, default=None,
                       help="プライマー候補のディレクトリ")
    parser.add_argument("--output-dir", type=str, default=None,
                       help="出力ディレクトリ")
    args = parser.parse_args()

    # ディレクトリ設定
    if args.sequences_dir:
        SEQUENCES_DIR = Path(args.sequences_dir)
    if args.primers_dir:
        PRIMERS_DIR = Path(args.primers_dir)
    if args.output_dir:
        OUTPUT_DIR = Path(args.output_dir)

    print("=" * 60)
    print("GCAプライマー クロスチェック解析")
    print("=" * 60)

    # 出力ディレクトリ作成
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # データ読み込み
    reference_sequences = load_reference_sequences()
    if not reference_sequences:
        print("エラー: 参照配列が見つかりませんでした。先に fetch_sequences.py を実行してください。")
        sys.exit(1)
        
    primers_to_check = load_primer_candidates()
    if not primers_to_check:
        print("エラー: プライマー候補が見つかりませんでした。先に design_primers_primer3.py を実行してください。")
        sys.exit(1)

    # 全プライマーセットをチェック
    all_results = []
    for primer_name, primer_data in primers_to_check.items():
        print(f"\n解析中: {primer_name}...")
        result = crosscheck_primer_set(primer_name, primer_data, reference_sequences)
        all_results.append(result)

        # ターゲット結果を抽出して表示（ターゲットが複数ある場合は最初のものを表示）
        target_results = [res for name, res in result['organisms'].items() if res['category'] == 'target']
        if target_results:
            target_result = target_results[0]
            print(f"  標的マッチ: F={target_result.get('best_forward_mismatch', 'N/A')}mm, "
                  f"R={target_result.get('best_reverse_mismatch', 'N/A')}mm")
        else:
            print("  警告: 標的配列の解析結果がありません。")

    # レポート生成
    print("\nレポート生成中...")

    # テキストレポート
    report_text = generate_crosscheck_report(all_results)
    report_file = os.path.join(OUTPUT_DIR, "crosscheck_report.txt")
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report_text)
    print(f"テキストレポート: {report_file}")

    # CSVレポート
    csv_file = os.path.join(OUTPUT_DIR, "crosscheck_results.csv")
    generate_csv_report(all_results, csv_file)
    print(f"CSVレポート: {csv_file}")

    # コンソールにも表示
    print("\n" + report_text)


if __name__ == "__main__":
    main()
