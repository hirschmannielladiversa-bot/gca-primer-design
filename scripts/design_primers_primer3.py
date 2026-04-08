#!/usr/bin/env python3
"""
Primer3を使用したGCA用プライマー設計スクリプト
標的種のCOI/matK等の遺伝子をターゲット
"""

import os
import sys
import argparse
from pathlib import Path

try:
    import primer3
except ImportError:
    print("primer3-pyがインストールされていません。")
    print("インストール: pip install primer3-py")
    sys.exit(1)

try:
    from Bio import SeqIO, AlignIO
    from Bio.Seq import Seq
except ImportError:
    print("BioPythonがインストールされていません。")
    print("インストール: pip install biopython")
    sys.exit(1)

# パス設定（このスクリプトの親の親ディレクトリを基準）
SCRIPT_DIR = Path(__file__).parent.absolute()
PROJECT_ROOT = SCRIPT_DIR.parent

# デフォルトのディレクトリ
DEFAULT_SEQUENCES_DIR = PROJECT_ROOT / "sequences"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "primer_candidates"

# グローバル変数（mainで上書き可能）
SEQUENCES_DIR = DEFAULT_SEQUENCES_DIR
OUTPUT_DIR = DEFAULT_OUTPUT_DIR


def calculate_tm(seq):
    """Nearest-neighbor法によるTm計算（簡易版）"""
    import math

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

    # 塩濃度補正（50 mM Na+）
    dS_corrected = dS + 0.368 * (len(seq) - 1) * math.log(0.05)
    tm = (dH * 1000) / (dS_corrected + 1.987 * math.log(250e-9 / 4)) - 273.15

    return tm


def calculate_gc(seq):
    """GC含量を計算"""
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return (gc / len(seq)) * 100


def design_primers_with_primer3(sequence, seq_id="target", at_rich=False):
    """
    Primer3を使用してプライマーを設計
    GCA用に最適化されたパラメータを使用

    Parameters:
        sequence: テンプレート配列
        seq_id: 配列ID
        at_rich: AT-rich配列用パラメータを使用（Lepidoptera等のミトコンドリアDNA向け）
    """
    if at_rich:
        # AT-rich配列用パラメータ（Lepidoptera, Diptera等のミトコンドリアDNA）
        # 参考: LepF1/LepR1プライマー条件、PMC7379581
        primer3_params = {
            'SEQUENCE_ID': seq_id,
            'SEQUENCE_TEMPLATE': sequence,

            # プライマー設計パラメータ
            'PRIMER_TASK': 'generic',
            'PRIMER_PICK_LEFT_PRIMER': 1,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_PICK_INTERNAL_OLIGO': 0,

            # 産物サイズ（GCA用に短く）
            'PRIMER_PRODUCT_SIZE_RANGE': [[100, 140], [140, 180], [180, 220]],

            # プライマー長（長めでTmを補償）
            'PRIMER_MIN_SIZE': 20,
            'PRIMER_OPT_SIZE': 22,
            'PRIMER_MAX_SIZE': 25,

            # Tm（低めに設定 - AT-rich配列対応）
            'PRIMER_MIN_TM': 50.0,
            'PRIMER_OPT_TM': 55.0,
            'PRIMER_MAX_TM': 60.0,
            'PRIMER_TM_FORMULA': 1,  # SantaLucia法

            # GC含量（緩和 - AT-rich配列では30%程度が現実的）
            'PRIMER_MIN_GC': 25.0,
            'PRIMER_OPT_GC_PERCENT': 40.0,
            'PRIMER_MAX_GC': 55.0,

            # 3'末端のGC（AT-richでは困難なため緩和）
            'PRIMER_GC_CLAMP': 0,

            # 二次構造（緩和）
            'PRIMER_MAX_SELF_ANY': 8,
            'PRIMER_MAX_SELF_END': 4,
            'PRIMER_MAX_POLY_X': 5,  # AT-richではT/Aの連続が多い

            # ペア条件（緩和）
            'PRIMER_PAIR_MAX_DIFF_TM': 3.0,

            # 候補数
            'PRIMER_NUM_RETURN': 10,

            # 塩濃度
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 250.0,
        }
    else:
        # 標準パラメータ（GCA最適化）
        primer3_params = {
            'SEQUENCE_ID': seq_id,
            'SEQUENCE_TEMPLATE': sequence,

            # プライマー設計パラメータ
            'PRIMER_TASK': 'generic',
            'PRIMER_PICK_LEFT_PRIMER': 1,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_PICK_INTERNAL_OLIGO': 0,

            # 産物サイズ（GCA用に短く）
            'PRIMER_PRODUCT_SIZE_RANGE': [[100, 120], [120, 140], [140, 160]],

            # プライマー長
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MAX_SIZE': 22,

            # Tm（より高めに設定）
            'PRIMER_MIN_TM': 58.0,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_TM_FORMULA': 1,  # SantaLucia法

            # GC含量
            'PRIMER_MIN_GC': 45.0,
            'PRIMER_OPT_GC_PERCENT': 50.0,
            'PRIMER_MAX_GC': 55.0,

            # 3'末端のGC
            'PRIMER_GC_CLAMP': 1,  # 3'末端にG/Cを1つ以上

            # 二次構造
            'PRIMER_MAX_SELF_ANY': 6,
            'PRIMER_MAX_SELF_END': 3,
            'PRIMER_MAX_POLY_X': 3,

            # ペア条件
            'PRIMER_PAIR_MAX_DIFF_TM': 2.0,

            # 候補数
            'PRIMER_NUM_RETURN': 10,

            # 塩濃度
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 250.0,
        }

    # Primer3実行
    results = primer3.bindings.design_primers(
        seq_args={
            'SEQUENCE_ID': seq_id,
            'SEQUENCE_TEMPLATE': sequence,
        },
        global_args=primer3_params
    )

    return results


def calculate_pcr_conditions(primer, at_rich=False):
    """
    プライマー情報からPCR条件を計算

    ベストプラクティス参照:
    - Thermo Fisher: PCR Cycling Parameters
    - NEB: Guidelines for PCR Optimization
    - IDT: Annealing Temperature Calculation

    Args:
        primer: プライマー情報辞書
        at_rich: AT-rich配列モード（低Tm対応）

    Returns:
        PCR条件辞書
    """
    left_tm = primer.get('left_tm', 55)
    right_tm = primer.get('right_tm', 55)
    product_size = primer.get('product_size', 150)

    # アニーリング温度: min(Tm) - 5°C（ベストプラクティス）
    # AT-richモードでは -3°C に緩和
    min_tm = min(left_tm, right_tm)
    if at_rich:
        annealing_temp = round(min_tm - 3, 1)
        # AT-rich配列では45-55°Cが適切
        annealing_temp = max(45, min(55, annealing_temp))
    else:
        annealing_temp = round(min_tm - 5, 1)
        # 標準配列では50-65°Cが適切
        annealing_temp = max(50, min(65, annealing_temp))

    # 伸長時間: 1kb あたり 60秒（Taq）、短いアンプリコンは最低15秒
    # GCA用（100-200bp）は15-30秒で十分
    if product_size <= 200:
        extension_time = 15  # 秒
    elif product_size <= 500:
        extension_time = 30
    elif product_size <= 1000:
        extension_time = 60
    else:
        extension_time = max(60, int(product_size / 1000 * 60))

    # サイクル数: GCA用は30-35サイクル推奨
    cycles = 35

    # 初期変性: 95°C, 2分
    # 変性: 95°C, 30秒
    # アニーリング: 計算温度, 30秒
    # 伸長: 72°C, 計算時間
    # 最終伸長: 72°C, 5分

    pcr_conditions = {
        'initial_denaturation': {'temp': 95, 'time': 120},  # 秒
        'denaturation': {'temp': 95, 'time': 30},
        'annealing': {'temp': annealing_temp, 'time': 30},
        'extension': {'temp': 72, 'time': extension_time},
        'final_extension': {'temp': 72, 'time': 300},
        'cycles': cycles,
        'hold': {'temp': 4},
        # 追加情報
        'notes': [],
    }

    # 条件に応じた注意事項
    if at_rich:
        pcr_conditions['notes'].append("AT-rich配列: 低アニーリング温度を使用")
    if min_tm < 50:
        pcr_conditions['notes'].append("低Tm: タッチダウンPCRを検討")
    if abs(left_tm - right_tm) > 3:
        pcr_conditions['notes'].append(f"Tm差が大きい ({abs(left_tm - right_tm):.1f}°C): 勾配PCRを検討")

    return pcr_conditions


def parse_primer3_results(results, at_rich=False):
    """Primer3の結果をパース"""
    primers = []

    num_returned = results.get('PRIMER_PAIR_NUM_RETURNED', 0)

    for i in range(num_returned):
        primer = {
            'id': i + 1,
            'left_seq': results.get(f'PRIMER_LEFT_{i}_SEQUENCE', ''),
            'right_seq': results.get(f'PRIMER_RIGHT_{i}_SEQUENCE', ''),
            'left_tm': results.get(f'PRIMER_LEFT_{i}_TM', 0),
            'right_tm': results.get(f'PRIMER_RIGHT_{i}_TM', 0),
            'left_gc': results.get(f'PRIMER_LEFT_{i}_GC_PERCENT', 0),
            'right_gc': results.get(f'PRIMER_RIGHT_{i}_GC_PERCENT', 0),
            'product_size': results.get(f'PRIMER_PAIR_{i}_PRODUCT_SIZE', 0),
            'left_pos': results.get(f'PRIMER_LEFT_{i}', [0, 0]),
            'right_pos': results.get(f'PRIMER_RIGHT_{i}', [0, 0]),
            'penalty': results.get(f'PRIMER_PAIR_{i}_PENALTY', 0),
        }
        # PCR条件を計算して追加
        primer['pcr_conditions'] = calculate_pcr_conditions(primer, at_rich)
        primers.append(primer)

    return primers


def evaluate_primer(primer):
    """プライマーの品質評価"""
    score = 100

    # Tmチェック（58-62が理想）
    for tm in [primer['left_tm'], primer['right_tm']]:
        if tm < 58 or tm > 63:
            score -= 10
        if tm < 55 or tm > 65:
            score -= 20

    # Tm差（<2が理想）
    tm_diff = abs(primer['left_tm'] - primer['right_tm'])
    if tm_diff > 2:
        score -= 10 * (tm_diff - 2)

    # GC含量（45-55%が理想）
    for gc in [primer['left_gc'], primer['right_gc']]:
        if gc < 45 or gc > 55:
            score -= 5

    # 産物サイズ（100-120が最適）
    if primer['product_size'] > 140:
        score -= 5
    if primer['product_size'] > 160:
        score -= 10

    # 3'末端チェック
    for seq in [primer['left_seq'], primer['right_seq']]:
        if seq and seq[-1] not in 'GC':
            score -= 15

    primer['quality_score'] = max(0, score)
    return primer


def format_output(primers, output_file=None):
    """結果を整形して出力"""
    output = []
    output.append("=" * 70)
    output.append("Primer3 設計結果 - 標的種GCA用プライマー")
    output.append("=" * 70)
    output.append("")

    for p in primers:
        output.append(f"--- 候補 {p['id']} (品質スコア: {p.get('quality_score', 'N/A')}) ---")
        output.append(f"Forward: 5'-{p['left_seq']}-3' ({len(p['left_seq'])} bp)")
        output.append(f"Reverse: 5'-{p['right_seq']}-3' ({len(p['right_seq'])} bp)")
        output.append(f"Forward Tm: {p['left_tm']:.1f}°C, GC: {p['left_gc']:.1f}%")
        output.append(f"Reverse Tm: {p['right_tm']:.1f}°C, GC: {p['right_gc']:.1f}%")
        output.append(f"Tm差: {abs(p['left_tm'] - p['right_tm']):.1f}°C")
        output.append(f"産物サイズ: {p['product_size']} bp")
        output.append(f"3'末端: F={p['left_seq'][-1] if p['left_seq'] else 'N/A'}, "
                     f"R={p['right_seq'][-1] if p['right_seq'] else 'N/A'}")

        # PCR条件を追加
        if 'pcr_conditions' in p:
            pcr = p['pcr_conditions']
            output.append("")
            output.append("推奨PCR条件:")
            output.append(f"  初期変性: {pcr['initial_denaturation']['temp']}°C, {pcr['initial_denaturation']['time']//60}分")
            output.append(f"  [{pcr['cycles']}サイクル]")
            output.append(f"    変性: {pcr['denaturation']['temp']}°C, {pcr['denaturation']['time']}秒")
            output.append(f"    アニーリング: {pcr['annealing']['temp']}°C, {pcr['annealing']['time']}秒")
            output.append(f"    伸長: {pcr['extension']['temp']}°C, {pcr['extension']['time']}秒")
            output.append(f"  最終伸長: {pcr['final_extension']['temp']}°C, {pcr['final_extension']['time']//60}分")
            if pcr.get('notes'):
                output.append(f"  注意: {'; '.join(pcr['notes'])}")

        output.append("")

    result_text = "\n".join(output)

    if output_file:
        with open(output_file, 'w') as f:
            f.write(result_text)
        print(f"結果保存: {output_file}")

    return result_text


def export_json(primers, target_info, output_file):
    """
    プライマー情報をJSON形式でエクスポート
    HTML側で読み込み可能な形式

    Args:
        primers: プライマーリスト
        target_info: ターゲット情報辞書
        output_file: 出力ファイルパス
    """
    import json
    from datetime import datetime

    # エクスポート用データ構造
    export_data = {
        'metadata': {
            'generated_at': datetime.now().isoformat(),
            'generator': 'GCA Primer Auto Design Package',
            'version': '1.3',
        },
        'target': target_info,
        'primers': []
    }

    for p in primers:
        primer_data = {
            'id': p['id'],
            'name': f"{target_info.get('prefix', 'Primer')}-{p['id']}",
            'forward': {
                'sequence': p['left_seq'],
                'tm': round(p['left_tm'], 1),
                'gc': round(p['left_gc'], 1),
                'length': len(p['left_seq']),
            },
            'reverse': {
                'sequence': p['right_seq'],
                'tm': round(p['right_tm'], 1),
                'gc': round(p['right_gc'], 1),
                'length': len(p['right_seq']),
            },
            'product_size': p['product_size'],
            'quality_score': p.get('quality_score', 0),
            'pcr_conditions': p.get('pcr_conditions', {}),
        }
        export_data['primers'].append(primer_data)

    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(export_data, f, ensure_ascii=False, indent=2)

    print(f"JSON保存: {output_file}")
    return export_data


def main():
    global SEQUENCES_DIR, OUTPUT_DIR

    # コマンドライン引数パーサー
    parser = argparse.ArgumentParser(description="Primer3 GCA用プライマー設計")
    parser.add_argument("--input-dir", type=str, default=None,
                       help="配列ファイルのディレクトリ")
    parser.add_argument("--output-dir", type=str, default=None,
                       help="出力ディレクトリ")
    parser.add_argument("--at-rich", action="store_true",
                       help="AT-rich配列用パラメータを使用（Lepidoptera等のミトコンドリアDNA向け）")
    args = parser.parse_args()

    # ディレクトリ設定
    if args.input_dir:
        SEQUENCES_DIR = Path(args.input_dir)
    if args.output_dir:
        OUTPUT_DIR = Path(args.output_dir)

    print("=" * 60)
    print("Primer3 GCA用プライマー設計")
    if args.at_rich:
        print("【AT-rich モード】Lepidoptera等のミトコンドリアDNA用パラメータを使用")
    print("=" * 60)

    # 出力ディレクトリ作成
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 配列ファイルを読み込む
    af_fasta = SEQUENCES_DIR / "target_sequence.fasta"

    if af_fasta.exists():
        print(f"配列ファイル読み込み: {af_fasta}")
        records = list(SeqIO.parse(af_fasta, "fasta"))
        if records:
            sequence = str(records[0].seq)
            seq_id = records[0].id
            print(f"配列ID: {seq_id}")
        else:
            print("エラー: 配列ファイルが空です。")
            print("先にfetch_sequences.pyを実行して配列を取得してください。")
            sys.exit(1)
    else:
        print(f"エラー: 配列ファイルが見つかりません: {af_fasta}")
        print("先にfetch_sequences.pyを実行して配列を取得してください。")
        sys.exit(1)

    print(f"配列長: {len(sequence)} bp")

    # Primer3で設計
    print("\nPrimer3でプライマー設計中...")
    results = design_primers_with_primer3(sequence, seq_id=seq_id, at_rich=args.at_rich)

    # 結果パース（at_richフラグを渡してPCR条件計算に反映）
    primers = parse_primer3_results(results, at_rich=args.at_rich)

    if not primers:
        print("プライマー候補が見つかりませんでした。")
        if not args.at_rich:
            print("\n【ヒント】AT-rich配列（Lepidoptera等のミトコンドリアDNA）の場合:")
            print("  --at-rich オプションを使用して再試行してください。")
            print("  例: python design_primers_primer3.py --at-rich --input-dir <dir>")
        else:
            print("AT-richモードでも候補が見つかりませんでした。")
            print("配列を確認してください。")

        # エラー情報
        if 'PRIMER_ERROR' in results:
            print(f"\nエラー: {results['PRIMER_ERROR']}")
        if 'PRIMER_WARNING' in results:
            print(f"警告: {results['PRIMER_WARNING']}")
        sys.exit(1)

    print(f"\n{len(primers)}個のプライマー候補が見つかりました。")

    # 品質評価
    primers = [evaluate_primer(p) for p in primers]

    # スコアでソート
    primers.sort(key=lambda x: x.get('quality_score', 0), reverse=True)

    # 結果出力
    output_file = OUTPUT_DIR / "primer_candidates.txt"
    result_text = format_output(primers, str(output_file))
    print(result_text)

    # CSV形式でも保存
    csv_file = OUTPUT_DIR / "primer_candidates.csv"
    with open(csv_file, 'w') as f:
        f.write("ID,Forward,Reverse,F_Tm,R_Tm,F_GC,R_GC,Product_Size,Quality_Score,Annealing_Temp\n")
        for p in primers:
            annealing = p.get('pcr_conditions', {}).get('annealing', {}).get('temp', '')
            f.write(f"{p['id']},{p['left_seq']},{p['right_seq']},"
                   f"{p['left_tm']:.1f},{p['right_tm']:.1f},"
                   f"{p['left_gc']:.1f},{p['right_gc']:.1f},"
                   f"{p['product_size']},{p.get('quality_score', 0)},{annealing}\n")
    print(f"CSV保存: {csv_file}")

    # JSON形式でエクスポート（HTML連携用）
    # プレフィックスをディレクトリ名から生成
    input_dir_name = SEQUENCES_DIR.parent.name if SEQUENCES_DIR.parent else "unknown"
    parts = input_dir_name.replace('_', ' ').split()
    if len(parts) >= 2:
        prefix = parts[0][0].upper() + parts[1][0].lower() if len(parts[1]) > 0 else parts[0][:2]
        gene = parts[-1] if len(parts) > 2 else "COI"
        prefix = f"{prefix}-{gene}"
    else:
        prefix = "Primer"

    target_info = {
        'species': input_dir_name.replace('_', ' ').rsplit(' ', 1)[0] if '_' in input_dir_name else input_dir_name,
        'gene': parts[-1] if len(parts) > 1 else "unknown",
        'sequence_id': seq_id,
        'sequence_length': len(sequence),
        'prefix': prefix,
        'at_rich_mode': args.at_rich,
    }

    json_file = OUTPUT_DIR / "primer_candidates.json"
    export_json(primers, target_info, str(json_file))


if __name__ == "__main__":
    main()
