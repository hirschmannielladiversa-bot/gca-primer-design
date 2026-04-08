#!/usr/bin/env python3
"""
GCAプライマー最適化 包括的レポート生成スクリプト
- プライマー特性分析
- クロスチェック結果
- 推奨PCR条件
- 総合評価

改善点:
- ハードコードされたArtemiaデータを削除
- 実際に取得した配列・設計したプライマーのデータを使用
- セクション2（プライマー特性分析）を追加
"""

import os
import sys
import math
import argparse
import csv
from datetime import datetime
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:
    print("BioPythonがインストールされていません。")
    print("インストール: pip install biopython")
    sys.exit(1)

# パス設定（このスクリプトの親の親ディレクトリを基準）
SCRIPT_DIR = Path(__file__).parent.absolute()
PROJECT_ROOT = SCRIPT_DIR.parent

# デフォルトのディレクトリ
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "reports"
DEFAULT_PRIMER_DIR = PROJECT_ROOT / "primer_candidates"
DEFAULT_CROSSCHECK_DIR = PROJECT_ROOT / "crosscheck_results"
DEFAULT_SEQUENCES_DIR = PROJECT_ROOT / "sequences"

# グローバル変数（mainで上書き可能）
BASE_DIR = PROJECT_ROOT
OUTPUT_DIR = DEFAULT_OUTPUT_DIR
PRIMER_DIR = DEFAULT_PRIMER_DIR
CROSSCHECK_DIR = DEFAULT_CROSSCHECK_DIR
SEQUENCES_DIR = DEFAULT_SEQUENCES_DIR


def calculate_gc(seq):
    """GC含量を計算"""
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return (gc / len(seq)) * 100


def calculate_tm_nn(seq, na_conc=50, primer_conc=250):
    """Nearest-neighbor法によるTm計算"""
    nn_params = {
        'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2),
        'AT': (-7.2, -20.4), 'TA': (-7.2, -21.3),
        'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7),
        'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
        'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0),
        'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
        'CG': (-10.6, -27.2), 'GC': (-9.8, -24.4),
        'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9)
    }

    seq = seq.upper()
    dH = 0
    dS = 0

    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        if dinuc in nn_params:
            dH += nn_params[dinuc][0]
            dS += nn_params[dinuc][1]

    dS_corrected = dS + 0.368 * (len(seq) - 1) * math.log(na_conc / 1000)
    tm = (dH * 1000) / (dS_corrected + 1.987 * math.log(primer_conc / 4e9)) - 273.15
    return tm


def get_target_sequence_info():
    """標的配列の情報を取得"""
    target_fasta = SEQUENCES_DIR / "target_sequence.fasta"
    if not target_fasta.exists():
        return None

    try:
        records = list(SeqIO.parse(target_fasta, "fasta"))
        if records:
            first_record = records[0]
            return {
                'accession': first_record.id,
                'description': first_record.description,
                'length': len(first_record.seq),
                'sequence': str(first_record.seq)[:100] + "..." if len(first_record.seq) > 100 else str(first_record.seq),
                'num_sequences': len(records)
            }
    except Exception as e:
        print(f"Warning: Could not read target_sequence.fasta: {e}")
    return None


def load_primer_candidates():
    """プライマー候補を読み込む"""
    primer_csv = PRIMER_DIR / "primer_candidates.csv"
    primers = []

    if primer_csv.exists():
        try:
            with open(primer_csv, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    primers.append(row)
        except Exception as e:
            print(f"Warning: Could not read primer_candidates.csv: {e}")

    return primers


def load_crosscheck_results():
    """クロスチェック結果を読み込む"""
    crosscheck_csv = CROSSCHECK_DIR / "crosscheck_results.csv"
    results = []

    if crosscheck_csv.exists():
        try:
            with open(crosscheck_csv, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    results.append(row)
        except Exception as e:
            print(f"Warning: Could not read crosscheck_results.csv: {e}")

    return results


def generate_report(target_species, target_gene):
    """Markdownレポート生成"""
    report = []

    # target_speciesが指定されていない場合のデフォルト
    if not target_species:
        target_species = "Target Species"
    if not target_gene:
        target_gene = "Target Gene"

    # データ読み込み
    primer_candidates = load_primer_candidates()
    crosscheck_results = load_crosscheck_results()
    target_seq_info = get_target_sequence_info()

    # ヘッダー
    report.append(f"# {target_species} GCAプライマー({target_gene})最適化レポート")
    report.append("")
    report.append(f"**レポート生成日**: {datetime.now().strftime('%Y年%m月%d日 %H:%M')}")
    report.append("")
    report.append(f"**プロジェクト**: Gut Content Analysis (GCA) - {target_species} 検出用プライマー最適化")
    report.append("")
    report.append("---")
    report.append("")

    # 目次
    report.append("## 目次")
    report.append("")
    report.append("1. [エグゼクティブサマリー](#1-エグゼクティブサマリー)")
    report.append("2. [プライマー特性分析](#2-プライマー特性分析)")
    report.append("3. [クロスチェック結果](#3-クロスチェック結果)")
    report.append("4. [推奨PCR条件](#4-推奨pcr条件)")
    report.append("5. [総合評価と推奨事項](#5-総合評価と推奨事項)")
    report.append("6. [次のステップ](#6-次のステップ)")
    report.append("")
    report.append("---")
    report.append("")

    # 1. エグゼクティブサマリー
    report.append("## 1. エグゼクティブサマリー")
    report.append("")

    # 使用した配列情報
    report.append("### 使用した標的配列")
    report.append("")
    if target_seq_info:
        report.append(f"- **アクセッション番号**: {target_seq_info['accession']}")
        report.append(f"- **説明**: {target_seq_info['description']}")
        report.append(f"- **配列長**: {target_seq_info['length']} bp")
        report.append(f"- **取得配列数**: {target_seq_info['num_sequences']} 件")
    else:
        report.append("⚠ 標的配列ファイル(target_sequence.fasta)が見つかりませんでした。")
    report.append("")

    report.append("### プライマー設計結果サマリー")
    report.append("")

    if primer_candidates:
        report.append(f"{len(primer_candidates)}個のプライマー候補が設計されました。")
        report.append("")
        report.append("| 推奨順位 | ID | 産物サイズ | Tm (F/R) | 品質スコア |")
        report.append("|----------|----|-----------|----------|------------|")
        for i, p in enumerate(primer_candidates[:5]):  # 上位5つを表示
            pid = p.get('ID', p.get('id', 'N/A'))
            psize = p.get('Product_Size', p.get('product_size', 'N/A'))
            ftm = p.get('F_Tm', p.get('left_tm', 'N/A'))
            rtm = p.get('R_Tm', p.get('right_tm', 'N/A'))
            score = p.get('Quality_Score', p.get('quality_score', 'N/A'))
            report.append(f"| {i+1} | {pid} | {psize} bp | {ftm}/{rtm}°C | {score} |")
    else:
        report.append("プライマー候補の設計結果ファイル(primer_candidates.csv)が見つかりませんでした。")
        report.append("Primer3による設計が失敗したか、スキップされた可能性があります。")

    report.append("")
    report.append("### 推奨アクション")
    report.append("")
    report.append("```")
    report.append("上位候補プライマーの合成検討")
    report.append("  → アニーリング温度はTm値に基づき 55-60°C を基準に試験")
    report.append("  → In vitroでの特異性・感度検証を実施")
    report.append("```")
    report.append("")
    report.append("---")
    report.append("")

    # 2. プライマー特性分析（新規追加）
    report.append("## 2. プライマー特性分析")
    report.append("")

    if primer_candidates:
        report.append("### 2.1 設計プライマー一覧")
        report.append("")
        report.append("| ID | Forward (5'→3') | Reverse (5'→3') | 産物サイズ |")
        report.append("|-----|-----------------|-----------------|-----------|")

        for p in primer_candidates[:10]:  # 上位10個
            pid = p.get('ID', p.get('id', 'N/A'))
            fwd = p.get('Forward', p.get('forward', 'N/A'))
            rev = p.get('Reverse', p.get('reverse', 'N/A'))
            psize = p.get('Product_Size', p.get('product_size', 'N/A'))
            report.append(f"| {pid} | {fwd} | {rev} | {psize} bp |")

        report.append("")
        report.append("### 2.2 プライマー熱力学特性")
        report.append("")
        report.append("| ID | F長 | R長 | F Tm | R Tm | Tm差 | F GC% | R GC% | F 3'末端 | R 3'末端 |")
        report.append("|-----|-----|-----|------|------|------|-------|-------|----------|----------|")

        for p in primer_candidates[:10]:
            pid = p.get('ID', p.get('id', 'N/A'))
            fwd = p.get('Forward', p.get('forward', ''))
            rev = p.get('Reverse', p.get('reverse', ''))
            ftm = p.get('F_Tm', p.get('left_tm', 'N/A'))
            rtm = p.get('R_Tm', p.get('right_tm', 'N/A'))
            fgc = p.get('F_GC', p.get('left_gc', 'N/A'))
            rgc = p.get('R_GC', p.get('right_gc', 'N/A'))

            f_len = len(fwd) if fwd != 'N/A' else 'N/A'
            r_len = len(rev) if rev != 'N/A' else 'N/A'
            f_3prime = fwd[-1] if fwd and fwd != 'N/A' else 'N/A'
            r_3prime = rev[-1] if rev and rev != 'N/A' else 'N/A'

            try:
                tm_diff = abs(float(ftm) - float(rtm))
                tm_diff_str = f"{tm_diff:.1f}"
            except:
                tm_diff_str = "N/A"

            report.append(f"| {pid} | {f_len} | {r_len} | {ftm}°C | {rtm}°C | {tm_diff_str}°C | {fgc}% | {rgc}% | {f_3prime} | {r_3prime} |")

        report.append("")
        report.append("### 2.3 プライマー評価基準")
        report.append("")
        report.append("| 項目 | 理想値 | 許容範囲 |")
        report.append("|------|--------|----------|")
        report.append("| プライマー長 | 20 bp | 18-22 bp |")
        report.append("| Tm | 60°C | 58-63°C |")
        report.append("| Tm差（F/R） | 0°C | <2°C |")
        report.append("| GC含量 | 50% | 45-55% |")
        report.append("| 3'末端 | G/C | G/C推奨 |")
        report.append("| 産物サイズ（GCA用） | 100-120 bp | 100-150 bp |")
        report.append("")
    else:
        report.append("プライマー候補が見つかりませんでした。")
        report.append("")

    report.append("---")
    report.append("")

    # 3. クロスチェック結果
    report.append("## 3. クロスチェック結果")
    report.append("")

    if crosscheck_results:
        report.append("### 交差反応性評価サマリー")
        report.append("以下は標的以外の生物とのアライメント（ミスマッチ数）に基づく予測リスクです。")
        report.append("")
        report.append("| プライマー | 生物 | カテゴリ | リスク等級 |")
        report.append("|-----------|------|---------|----------|")

        # 標的以外の結果を表示
        nontargets = [r for r in crosscheck_results if r.get('Category', '') != 'target']
        for r in nontargets[:10]:  # 最大10件
            report.append(f"| {r.get('Primer_Set')} | {r.get('Organism')} | {r.get('Category')} | {r.get('Risk_Level')} |")

        report.append("")
        report.append("詳細なクロスチェック結果は `crosscheck_results/crosscheck_report.txt` または `crosscheck_results.csv` を確認してください。")
    else:
        report.append("クロスチェック結果のファイル(crosscheck_results.csv)が見つかりませんでした。")

    report.append("")
    report.append("---")
    report.append("")

    # 4. 推奨PCR条件
    report.append("## 4. 推奨PCR条件")
    report.append("")

    # Tm値に基づいてアニーリング温度を推奨
    recommended_annealing = "55-58"
    if primer_candidates:
        try:
            ftm = float(primer_candidates[0].get('F_Tm', primer_candidates[0].get('left_tm', 60)))
            rtm = float(primer_candidates[0].get('R_Tm', primer_candidates[0].get('right_tm', 60)))
            avg_tm = (ftm + rtm) / 2
            recommended_annealing = f"{avg_tm - 5:.0f}-{avg_tm - 2:.0f}"
        except:
            pass

    report.append(f"### 4.1 推奨アニーリング温度: {recommended_annealing}°C")
    report.append("")
    report.append("設計したプライマーのTm値に基づく推奨値です。")
    report.append("")
    report.append("#### 標準PCRプロトコル")
    report.append("")
    report.append("```")
    report.append("PCRサイクル条件（推奨）:")
    report.append("┌─────────────────────────────────────────┐")
    report.append("│ 1. 初期変性    95°C   3分              │")
    report.append("│ 2. 変性        95°C   30秒             │")
    report.append(f"│ 3. アニーリング {recommended_annealing}°C  30秒       │")
    report.append("│ 4. 伸長        72°C   30秒             │")
    report.append("│ 5. ステップ2-4を 35-40サイクル         │")
    report.append("│ 6. 最終伸長    72°C   5分              │")
    report.append("│ 7. 保持        4°C                     │")
    report.append("└─────────────────────────────────────────┘")
    report.append("```")
    report.append("")

    report.append("#### タッチダウンPCR（推奨）")
    report.append("")
    report.append("```")
    report.append("タッチダウンPCRサイクル条件:")
    report.append("┌─────────────────────────────────────────┐")
    report.append("│ Phase 1: 初期変性                       │")
    report.append("│   95°C  3分                             │")
    report.append("│                                         │")
    report.append("│ Phase 2: タッチダウン（10サイクル）      │")
    report.append("│   95°C  30秒                            │")
    report.append("│   65°C → 55°C  30秒（-1°C/サイクル）    │")
    report.append("│   72°C  30秒                            │")
    report.append("│                                         │")
    report.append("│ Phase 3: アニーリング固定（30サイクル）  │")
    report.append("│   95°C  30秒                            │")
    report.append("│   55°C  30秒                            │")
    report.append("│   72°C  30秒                            │")
    report.append("│                                         │")
    report.append("│ Phase 4: 最終伸長                       │")
    report.append("│   72°C  5分                             │")
    report.append("└─────────────────────────────────────────┘")
    report.append("```")
    report.append("")

    report.append("### 4.2 反応組成（20 μL系）")
    report.append("")
    report.append("| 試薬 | 濃度 | 体積 |")
    report.append("|------|------|------|")
    report.append("| 2× PCR Master Mix (Hot-start) | 1× | 10 μL |")
    report.append("| Forward Primer (10 μM) | 0.5 μM | 1 μL |")
    report.append("| Reverse Primer (10 μM) | 0.5 μM | 1 μL |")
    report.append("| Template DNA | 1-10 ng | 2 μL |")
    report.append("| **DMSO (オプション)** | **3%** | **0.6 μL** |")
    report.append("| H2O | - | 5.4 μL |")
    report.append("| **Total** | | **20 μL** |")
    report.append("")

    report.append("### 4.3 添加剤による改善")
    report.append("")
    report.append("| 添加剤 | 濃度 | 効果 | 推奨度 |")
    report.append("|--------|------|------|--------|")
    report.append("| **DMSO** | 3-5% | 二次構造解消、特異性向上 | ★★★ 強く推奨 |")
    report.append("| Betaine | 1M | AT-rich領域の増幅改善 | ★★ 推奨 |")
    report.append("| BSA | 0.1-0.5 mg/mL | 阻害物質対策（GCAサンプル用） | ★★ 推奨 |")
    report.append("| MgCl2追加 | +0.5-1.0 mM | プライマー結合安定化 | ★ 状況に応じて |")
    report.append("")
    report.append("---")
    report.append("")

    # 5. 総合評価と推奨事項
    report.append("## 5. 総合評価と推奨事項")
    report.append("")

    report.append("### 5.1 プライマー推奨順位")
    report.append("")
    if primer_candidates:
        report.append("| 順位 | プライマー | 品質スコア | 推奨理由 |")
        report.append("|------|-----------|----------|----------|")
        for i, p in enumerate(primer_candidates[:3]):  # 上位3つ
            pid = p.get('ID', p.get('id', 'N/A'))
            score = p.get('Quality_Score', p.get('quality_score', 'N/A'))

            # 推奨理由を生成
            reasons = []
            try:
                ftm = float(p.get('F_Tm', p.get('left_tm', 0)))
                rtm = float(p.get('R_Tm', p.get('right_tm', 0)))
                if 58 <= ftm <= 63 and 58 <= rtm <= 63:
                    reasons.append("適正Tm")
                if abs(ftm - rtm) < 2:
                    reasons.append("Tm差小")
            except:
                pass

            fwd = p.get('Forward', p.get('forward', ''))
            rev = p.get('Reverse', p.get('reverse', ''))
            if fwd and fwd[-1] in 'GC' and rev and rev[-1] in 'GC':
                reasons.append("3'末端GC")

            reason_str = ", ".join(reasons) if reasons else "高スコア候補"
            report.append(f"| {i+1} | {pid} | {score} | {reason_str} |")
    else:
        report.append("具体的な推奨順位はありません。新規設計の実施が必要です。")

    report.append("")

    report.append("### 5.2 推奨アクション")
    report.append("")
    report.append("#### 即時")
    report.append("")
    report.append("1. 上位プライマー候補の合成を発注")
    report.append("2. 温度勾配実験で最適アニーリング温度を決定")
    report.append("3. DMSO添加条件を併せて検証")
    report.append("")

    report.append("#### 短期")
    report.append("")
    report.append("1. 感度試験（DNA希釈系列）を実施")
    report.append("2. GCAサンプルでの実証実験")
    report.append("3. 近縁種との交差反応を実験的に確認")
    report.append("")

    report.append("### 5.3 注意事項")
    report.append("")
    report.append(f"- **近縁種との交差反応**: {target_species}の近縁種や共存生物との交差反応リスクを `crosscheck_results.csv` や `crosscheck_report.txt` で確認してください。")
    report.append("- **GCAサンプルのDNA品質**: 分解が進んでいる場合、短い産物（<150bp）が有利です。")
    report.append("- 本レポートはPrimer3によるIn-silico設計結果に基づいています。最終的には実験的検証（PCR）が必要です。")
    report.append("")
    report.append("---")
    report.append("")

    # 6. 次のステップ
    report.append("## 6. 次のステップ")
    report.append("")
    report.append("### 実験計画")
    report.append("")
    report.append("```")
    report.append("Phase 1: プライマー検証")
    report.append("  ├── 温度勾配実験")
    report.append("  ├── 感度試験")
    report.append("  └── 添加剤効果検証")
    report.append("")
    report.append("Phase 2: 特異性確認")
    report.append("  ├── 近縁種DNAでのPCR")
    report.append("  ├── ヒトDNAでのPCR")
    report.append("  └── 混合サンプルでのPCR")
    report.append("")
    report.append("Phase 3: 実サンプル検証")
    report.append("  ├── GCAサンプルでのPCR")
    report.append("  └── 定量性評価")
    report.append("```")
    report.append("")

    report.append("### 必要機材・試薬")
    report.append("")
    report.append("| 品目 | 用途 | 備考 |")
    report.append("|------|------|------|")
    report.append("| Hot-start Taq | PCR | 非特異的増幅抑制 |")
    report.append("| DMSO | 添加剤 | 3-5%使用 |")
    report.append("| 2%アガロースゲル | 電気泳動 | 100-150bp産物確認用 |")
    report.append("| 100bp DNAラダー | マーカー | - |")
    report.append("")
    report.append("---")
    report.append("")

    # 付録
    report.append("## 付録")
    report.append("")
    report.append("### プライマー配列一覧（発注用）")
    report.append("")
    report.append("```")

    if primer_candidates:
        for i, p in enumerate(primer_candidates[:5]):  # 上位5個
            pid = p.get('ID', p.get('id', i+1))
            fwd = p.get('Forward', p.get('forward', 'N/A'))
            rev = p.get('Reverse', p.get('reverse', 'N/A'))
            psize = p.get('Product_Size', p.get('product_size', 'N/A'))

            fwd_name = f"{target_species.replace(' ', '-')}-{target_gene}-{pid}F"
            rev_name = f"{target_species.replace(' ', '-')}-{target_gene}-{pid}R"

            report.append(f"Primer Set {pid}:")
            report.append(f"  {fwd_name}: 5'-{fwd}-3' ({len(fwd)} bp)")
            report.append(f"  {rev_name}: 5'-{rev}-3' ({len(rev)} bp)")
            report.append(f"  予測産物サイズ: {psize} bp")
            report.append("")
    else:
        report.append("プライマー候補がありません。")

    report.append("```")
    report.append("")

    report.append("### 使用した参照配列")
    report.append("")
    if target_seq_info:
        report.append(f"- **標的配列**: [{target_seq_info['accession']}](https://www.ncbi.nlm.nih.gov/nuccore/{target_seq_info['accession']})")
        report.append(f"  - {target_seq_info['description']}")
    else:
        report.append("- 標的配列情報なし")
    report.append("")

    report.append("### 参考リンク")
    report.append("")
    report.append("- NCBI Primer-BLAST: https://www.ncbi.nlm.nih.gov/tools/primer-blast/")
    report.append("- IDT OligoAnalyzer: https://www.idtdna.com/pages/tools/oligoanalyzer")
    report.append("- Primer3: https://primer3.ut.ee/")
    report.append("")
    report.append("---")
    report.append("")
    report.append(f"*レポート生成: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*")
    report.append("")
    report.append(f"*プロジェクトディレクトリ: `{PROJECT_ROOT}`*")

    return "\n".join(report)


def main():
    global OUTPUT_DIR, PRIMER_DIR, CROSSCHECK_DIR, SEQUENCES_DIR

    parser = argparse.ArgumentParser(description="GCAプライマー最適化 包括的レポート生成スクリプト")
    parser.add_argument("--target", default="Target Species", help="標的種")
    parser.add_argument("--gene", default="Target Gene", help="対象遺伝子 (例: COI, matK)")
    parser.add_argument("--project-dir", type=str, default=None,
                       help="プロジェクトディレクトリ（全データを含むフォルダ）")
    args = parser.parse_args()

    # プロジェクトディレクトリが指定されている場合、各ディレクトリを更新
    if args.project_dir:
        project_dir = Path(args.project_dir)
        SEQUENCES_DIR = project_dir / "sequences"
        PRIMER_DIR = project_dir / "primer_candidates"
        CROSSCHECK_DIR = project_dir / "crosscheck_results"
        OUTPUT_DIR = project_dir  # レポートはプロジェクトディレクトリ直下に出力

    print("=" * 60)
    print(f"GCAプライマー最適化 包括的レポート生成 ({args.target} - {args.gene})")
    print("=" * 60)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\nレポート生成中...")
    report = generate_report(args.target, args.gene)

    # Markdown保存（プロジェクトディレクトリ指定時は report.md として保存）
    if args.project_dir:
        output_file = OUTPUT_DIR / "report.md"
    else:
        filename_target = args.target.replace(" ", "_")
        output_file = OUTPUT_DIR / f"{filename_target}_Primer_Optimization_Report.md"

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(report)

    print(f"\n✓ レポート保存完了: {output_file}")

    # 確認用に最初の50行を表示
    print("\n" + "=" * 60)
    print("レポートプレビュー（最初の50行）:")
    print("=" * 60)
    for line in report.split('\n')[:50]:
        print(line)
    print("...")
    print(f"\n全文は {output_file} を参照してください。")


if __name__ == "__main__":
    main()
