#!/usr/bin/env python3
"""
NCBI Entrez APIを使用して配列を取得するスクリプト
標的種および比較対象生物の遺伝子配列を取得

改善点:
- COI遺伝子の複数の名称バリエーションに対応
- mitochondrion検索によるフォールバック
- スペルミス候補の提案
- 複数のアクセッション番号を確実に取得
- 異常に異なる配列を自動除外（品質フィルタリング）
"""

import os
import sys
import time
import argparse
from pathlib import Path
from collections import defaultdict

try:
    from Bio import Entrez, SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("BioPythonがインストールされていません。")
    print("インストール: pip install biopython")
    sys.exit(1)

# NCBI Entrez設定
# NCBI 利用規約上、有効なメールアドレスを設定することが推奨されています。
# 環境変数 GCA_ENTREZ_EMAIL で上書きできます。デフォルトは example.com の汎用値。
import os as _os
Entrez.email = _os.environ.get("GCA_ENTREZ_EMAIL", "gca_primer_design@example.com")

# デフォルト出力ディレクトリ（このスクリプトの親の親ディレクトリ/sequences）
SCRIPT_DIR = Path(__file__).parent.absolute()
PROJECT_ROOT = SCRIPT_DIR.parent
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "sequences"

# グローバル変数（main()で上書き可能）
OUTPUT_DIR = DEFAULT_OUTPUT_DIR

# 初期化時に削除するファイル（古いデータを防ぐ）
FILES_TO_CLEAN = [
    "target_sequence.fasta",
    "related_sequences.fasta",
    "nontarget_sequences.fasta",
]

# 遺伝子名のバリエーション（NCBIでの登録名が様々なため）
GENE_ALIASES = {
    "COI": [
        "COI",
        "CO1",
        "COXI",
        "COX1",
        "COX",
        "cytochrome oxidase subunit I",
        "cytochrome oxidase subunit 1",
        "cytochrome c oxidase subunit I",
        "cytochrome c oxidase subunit 1",
        "cytochrome oxidase",
    ],
    "COII": [
        "COII",
        "CO2",
        "COXII",
        "COX2",
        "cytochrome oxidase subunit II",
        "cytochrome c oxidase subunit II",
    ],
    "CYTB": [
        "cytb",
        "cyt b",
        "cytochrome b",
    ],
    "16S": [
        "16S",
        "16S rRNA",
        "16S ribosomal RNA",
    ],
    "18S": [
        "18S",
        "18S rRNA",
        "18S ribosomal RNA",
    ],
    "ITS": [
        "ITS",
        "ITS1",
        "ITS2",
        "internal transcribed spacer",
    ],
    # rRNA入力時は16S, 18S, ITSすべてで検索
    "RRNA": [
        "16S",
        "16S rRNA",
        "18S",
        "18S rRNA",
        "ITS",
        "ITS1",
        "ITS2",
        "ribosomal RNA",
        "rRNA",
    ],
    "MATK": [
        "matK",
        "maturase K",
    ],
    "RBCL": [
        "rbcL",
        "ribulose-1,5-bisphosphate carboxylase",
    ],
}

# COI検索時にmitochondrionも検索するフラグ
MITOCHONDRIAL_GENES = ["COI", "CO1", "COXI", "COX1", "COX", "COII", "CYTB", "ND1", "ND2", "ATP6"]

# 配列フィルタリングの閾値
MIN_SEQUENCE_LENGTH = 200  # 最小配列長
MAX_SEQUENCE_LENGTH = 20000  # 最大配列長
MIN_SIMILARITY_THRESHOLD = 0.15  # 絶対最小類似度（これ以下は明らかな外れ値）
OUTLIER_THRESHOLD = 0.20  # 平均からこれ以上離れていたら外れ値
MIN_SEQUENCES_TO_KEEP = 3  # 最低限保持する配列数


def clean_old_files():
    """古い配列ファイルを削除"""
    for filename in FILES_TO_CLEAN:
        filepath = OUTPUT_DIR / filename
        if filepath.exists():
            filepath.unlink()
            print(f"  削除: {filename}")


def fetch_sequence(accession, max_retries=3):
    """NCBIから配列を取得（リトライ機能付き）"""
    for attempt in range(max_retries):
        try:
            print(f"    取得中: {accession}...")
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            time.sleep(0.4)  # NCBI APIレート制限対策
            return record
        except Exception as e:
            error_str = str(e)
            # サーバーエラーの場合はリトライ
            if "Backend failed" in error_str or "500" in error_str or "503" in error_str:
                if attempt < max_retries - 1:
                    wait_time = (attempt + 1) * 2
                    print(f"    サーバーエラー、{wait_time}秒後にリトライ...")
                    time.sleep(wait_time)
                    continue
            print(f"    エラー: {accession} - {e}")
            return None
    return None


def get_gene_aliases(gene):
    """遺伝子名のエイリアスリストを取得"""
    gene_upper = gene.upper()

    # 既知のエイリアスがあればそれを使用
    for key, aliases in GENE_ALIASES.items():
        if gene_upper == key or gene_upper in [a.upper() for a in aliases]:
            return aliases

    # なければ元の名前のみ
    return [gene]


def search_sequences(species, gene_query, max_results=10, max_retries=3):
    """単一のクエリで検索（リトライ機能付き）"""
    for attempt in range(max_retries):
        try:
            search_handle = Entrez.esearch(
                db="nucleotide",
                term=gene_query,
                retmax=max_results,
                idtype="acc"
            )
            search_record = Entrez.read(search_handle)
            search_handle.close()
            time.sleep(0.4)  # レート制限対策
            return search_record.get("IdList", [])
        except Exception as e:
            error_str = str(e)
            # サーバーエラーの場合はリトライ
            if "Backend failed" in error_str or "500" in error_str or "503" in error_str:
                if attempt < max_retries - 1:
                    wait_time = (attempt + 1) * 2  # 2秒、4秒、6秒と増加
                    print(f"    サーバーエラー、{wait_time}秒後にリトライ...")
                    time.sleep(wait_time)
                    continue
            print(f"    検索エラー: {e}")
            return []
    return []


def calculate_sequence_similarity(seq1, seq2):
    """
    2つの配列の類似度を計算（シンプルなk-mer類似度）
    """
    k = 10  # k-merサイズ

    # 短い方の配列に合わせる
    min_len = min(len(seq1), len(seq2))
    if min_len < k:
        return 0.0

    # k-merセットを作成
    kmers1 = set(seq1[i:i+k] for i in range(len(seq1) - k + 1))
    kmers2 = set(seq2[i:i+k] for i in range(len(seq2) - k + 1))

    # Jaccard類似度
    intersection = len(kmers1 & kmers2)
    union = len(kmers1 | kmers2)

    if union == 0:
        return 0.0

    return intersection / union


def extract_species_from_description(description):
    """配列の説明から種名を抽出"""
    # 一般的なパターン: "ACCESSION.1 Species name gene description"
    parts = description.split()
    if len(parts) >= 3:
        # アクセッション番号の後の2-3単語が種名
        # 例: "OR773366.1 Frankliniella intonsa voucher..."
        species_parts = []
        for i, part in enumerate(parts[1:4]):  # 最大3単語
            # voucher, isolate, strain などのキーワードで終了
            if part.lower() in ['voucher', 'isolate', 'strain', 'clone', 'gene', 'cytochrome', 'mitochondrion', 'complete']:
                break
            species_parts.append(part)
        if species_parts:
            return ' '.join(species_parts)
    return "Unknown species"


def filter_outlier_sequences(records, min_similarity=MIN_SIMILARITY_THRESHOLD, species_name=None):
    """
    異常に異なる配列を除外する

    アルゴリズム:
    1. 全配列間のペアワイズ類似度を計算
    2. 各配列の平均類似度を計算
    3. 全体の平均類似度から大きく外れる配列を除外
    4. 最低限の配列数は保持する
    """
    if len(records) <= 2:
        return records, []

    print(f"\n  --- 配列品質フィルタリング ---")
    print(f"  対象種: {species_name or '不明'}")
    print(f"  入力配列数: {len(records)}")

    # 配列を文字列に変換
    sequences = {r.id: str(r.seq).upper() for r in records}

    # 長さフィルタリング
    valid_records = []
    removed_by_length = []
    for r in records:
        seq_len = len(r.seq)
        if seq_len < MIN_SEQUENCE_LENGTH:
            removed_by_length.append((r, f"配列長が短すぎる ({seq_len} bp < {MIN_SEQUENCE_LENGTH} bp)"))
        elif seq_len > MAX_SEQUENCE_LENGTH:
            removed_by_length.append((r, f"配列長が長すぎる ({seq_len} bp > {MAX_SEQUENCE_LENGTH} bp)"))
        else:
            valid_records.append(r)

    for r, reason in removed_by_length:
        species = extract_species_from_description(r.description)
        print(f"    ✗ 除外: {r.id} ({species}) - {reason}")

    if len(valid_records) <= 2:
        return valid_records, [r for r, _ in removed_by_length]

    # ペアワイズ類似度を計算
    similarities = defaultdict(list)
    ids = [r.id for r in valid_records]

    print(f"  類似度計算中 ({len(ids)} 配列)...")

    for i, id1 in enumerate(ids):
        for j, id2 in enumerate(ids):
            if i < j:
                sim = calculate_sequence_similarity(sequences[id1], sequences[id2])
                similarities[id1].append(sim)
                similarities[id2].append(sim)

    # 各配列の平均類似度
    avg_similarities = {id: sum(sims) / len(sims) if sims else 0 for id, sims in similarities.items()}

    # 全体の平均類似度
    all_sims = [s for sims in similarities.values() for s in sims]
    if not all_sims:
        return valid_records, [r for r, _ in removed_by_length]

    global_avg = sum(all_sims) / len(all_sims)
    print(f"  全体平均類似度: {global_avg:.2%}")

    # 配列を類似度でソート（高い順）
    sorted_records = sorted(valid_records, key=lambda r: avg_similarities.get(r.id, 0), reverse=True)

    # 外れ値を検出（相対的なフィルタリング）
    filtered_records = []
    potential_outliers = []

    for r in sorted_records:
        avg_sim = avg_similarities.get(r.id, 0)
        species = extract_species_from_description(r.description)

        # 絶対的に類似度が低い（明らかに異なる配列）
        if avg_sim < min_similarity:
            potential_outliers.append((r, f"平均類似度が極端に低い ({avg_sim:.2%})"))
        # 全体平均から大きく外れている
        elif global_avg - avg_sim > OUTLIER_THRESHOLD:
            potential_outliers.append((r, f"外れ値 (平均類似度 {avg_sim:.2%}, 全体平均 {global_avg:.2%})"))
        else:
            filtered_records.append(r)
            print(f"    ✓ 保持: {r.id} ({species}) - 類似度: {avg_sim:.2%}")

    # 最低限の配列数を確保（外れ値でも類似度の高い順に追加）
    if len(filtered_records) < MIN_SEQUENCES_TO_KEEP and potential_outliers:
        # 外れ値を類似度の高い順にソート
        potential_outliers.sort(key=lambda x: avg_similarities.get(x[0].id, 0), reverse=True)

        for r, reason in potential_outliers:
            if len(filtered_records) >= MIN_SEQUENCES_TO_KEEP:
                break
            filtered_records.append(r)
            avg_sim = avg_similarities.get(r.id, 0)
            species = extract_species_from_description(r.description)
            print(f"    ⚠ 復活（最低数確保）: {r.id} ({species}) - 類似度: {avg_sim:.2%}")
            potential_outliers = [(rr, rr_reason) for rr, rr_reason in potential_outliers if rr.id != r.id]

    # 実際に除外される配列
    removed_by_similarity = potential_outliers
    for r, reason in removed_by_similarity:
        species = extract_species_from_description(r.description)
        print(f"    ✗ 除外: {r.id} ({species}) - {reason}")

    all_removed = [r for r, _ in removed_by_length] + [r for r, _ in removed_by_similarity]

    print(f"  フィルタリング後: {len(filtered_records)} 配列 ({len(all_removed)} 配列除外)")

    return filtered_records, all_removed


def search_and_fetch_sequences(species, gene, max_results=10, min_sequences=3):
    """
    複数のクエリパターンで検索して配列を取得

    改善点:
    - 複数のクエリパターンで検索を継続（最初にヒットしても終了しない）
    - 重複アクセッション番号を排除
    - 異常に異なる配列を除外
    """

    # 遺伝子名のエイリアスを取得
    gene_aliases = get_gene_aliases(gene)
    gene_upper = gene.upper()

    # 検索クエリのパターン
    query_patterns = []

    # パターン1: [Gene]フィールドで検索（全エイリアス）
    for alias in gene_aliases:
        query_patterns.append(f"{species}[Organism] AND {alias}[Gene]")

    # パターン2: [Title]フィールドで検索（全エイリアス）
    for alias in gene_aliases:
        query_patterns.append(f"{species}[Organism] AND {alias}[Title]")

    # パターン3: mitochondrion/mitochondrial検索（ミトコンドリア遺伝子の場合）
    if gene_upper in MITOCHONDRIAL_GENES or "COI" in gene_upper or "COX" in gene_upper:
        query_patterns.append(f"{species}[Organism] AND mitochondrion[Title] AND complete[Title]")
        query_patterns.append(f"{species}[Organism] AND mitochondrion[Title]")
        query_patterns.append(f"{species}[Organism] AND mitochondrial[Title]")
        query_patterns.append(f"{species}[Organism] AND mitochondria[Title]")
        # COI/COX関連の追加検索
        query_patterns.append(f"{species}[Organism] AND cytochrome[Title] AND oxidase[Title]")

    # パターン4: rRNA関連の場合はribosomal検索も追加
    if gene_upper in ["RRNA", "16S", "18S", "ITS"]:
        query_patterns.append(f"{species}[Organism] AND ribosomal[Title]")
        query_patterns.append(f"{species}[Organism] AND rRNA[Title]")

    # 各クエリパターンを試して、重複を排除しながら収集
    all_ids = []
    query_hits = {}  # どのクエリでヒットしたかを記録

    for query in query_patterns:
        if len(all_ids) >= max_results * 2:  # 十分な候補を収集
            break

        print(f"  検索: {query[:70]}...")
        ids = search_sequences(species, query, max_results=max_results)

        if ids:
            new_ids = [id for id in ids if id not in all_ids]
            if new_ids:
                print(f"    → {len(new_ids)}件の新規ID（既存: {len(ids) - len(new_ids)}件は重複）")
                for id in new_ids:
                    all_ids.append(id)
                    query_hits[id] = query
            else:
                print(f"    → 全て重複（{len(ids)}件）")

        # 最小数に達していなければ検索を継続
        if len(all_ids) < min_sequences:
            continue

    if not all_ids:
        print(f"  ✗ {species} の {gene} 配列が見つかりませんでした。")
        print(f"    試したクエリ数: {len(query_patterns)}")
        return []

    print(f"\n  合計 {len(all_ids)} 件のユニークなアクセッション番号を取得")

    # 配列を取得
    records = []
    for acc in all_ids[:max_results]:
        record = fetch_sequence(acc)
        if record:
            records.append(record)

    if not records:
        print(f"  ✗ 配列の取得に失敗しました")
        return []

    print(f"  ✓ {len(records)}件の配列を取得")

    # 複数配列がある場合、品質フィルタリングを実行
    if len(records) > 2:
        filtered_records, removed_records = filter_outlier_sequences(records, species_name=species)

        if removed_records:
            # 除外された配列を別ファイルに保存（参考用）- 種名を含める
            # セキュリティ: 許可文字 [A-Za-z0-9._-] のみのサニタイズ
            import re
            species_safe = re.sub(r"[^A-Za-z0-9._-]", "_", species.strip().replace(" ", "_"))
            species_safe = species_safe.lstrip(".").replace("..", "_") or "unknown"
            gene_safe = re.sub(r"[^A-Za-z0-9._-]", "_", str(gene).strip()) or "gene"
            excluded_file = OUTPUT_DIR / f"excluded_{species_safe}_{gene_safe}.fasta"
            SeqIO.write(removed_records, excluded_file, "fasta")
            print(f"  除外配列を保存: {excluded_file}")

        return filtered_records
    else:
        return records


def suggest_species_name(species):
    """よくあるスペルミスを検出して提案"""
    import re

    # 完全な単語の置換のみ行う（部分一致で誤置換を防ぐ）
    # キー: 正規表現パターン（単語境界を使用）
    corrections = [
        (r"\bfranklinniella\b", "Frankliniella"),
        (r"\bfrankliniela\b", "Frankliniella"),
        (r"\bfranklinella\b", "Frankliniella"),
        (r"\bthrip\b", "Thrips"),  # 単数形のみ置換（thripsは置換しない）
        (r"\baphis\b", "Aphis"),
        (r"\bbemesia\b", "Bemisia"),
        (r"\btetranychus\b", "Tetranychus"),
    ]

    corrected = species
    was_corrected = False

    for pattern, correct in corrections:
        if re.search(pattern, corrected, re.IGNORECASE):
            corrected = re.sub(pattern, correct, corrected, flags=re.IGNORECASE)
            was_corrected = True

    # 修正があった場合のみ返す
    if was_corrected:
        # 属名は大文字開始、種小名は小文字
        parts = corrected.split()
        if parts:
            parts[0] = parts[0].capitalize()
            if len(parts) > 1:
                parts[1:] = [p.lower() for p in parts[1:]]
        return ' '.join(parts)

    return None


def verify_species_exists(species):
    """種名がNCBIに存在するか確認"""
    try:
        query = f"{species}[Organism]"
        search_handle = Entrez.esearch(db="taxonomy", term=query, retmax=1)
        search_record = Entrez.read(search_handle)
        search_handle.close()
        time.sleep(0.3)

        if search_record.get("IdList"):
            return True
        return False
    except:
        return None  # 確認できない場合


def save_sequences(sequences, filename):
    """配列をFASTA形式で保存"""
    filepath = OUTPUT_DIR / filename
    SeqIO.write(sequences, filepath, "fasta")
    print(f"保存完了: {filepath}")


def print_sequence_summary(sequences):
    """取得した配列のサマリーを表示"""
    print("\n  取得配列一覧:")
    print("  " + "-" * 70)
    print(f"  {'アクセッション':<20} {'長さ':>8} {'説明':<40}")
    print("  " + "-" * 70)
    for seq in sequences:
        desc = seq.description[:40] + "..." if len(seq.description) > 40 else seq.description
        print(f"  {seq.id:<20} {len(seq.seq):>8} bp {desc}")
    print("  " + "-" * 70)


def main():
    global OUTPUT_DIR

    parser = argparse.ArgumentParser(description="NCBI Entrez APIを使用して配列を取得するスクリプト")
    parser.add_argument("--target", required=True, help="標的種")
    parser.add_argument("--gene", required=True, help="対象遺伝子 (例: COI, matK)")
    parser.add_argument("--related", default="", help="近縁種 (カンマ区切り)")
    parser.add_argument("--nontarget", default="", help="非標的生物 (カンマ区切り)")
    parser.add_argument("--max-seqs", type=int, default=10, help="最大取得配列数 (デフォルト: 10)")
    parser.add_argument("--min-seqs", type=int, default=3, help="最小取得配列数 (デフォルト: 3)")
    parser.add_argument("--output-dir", type=str, default=None, help="出力ディレクトリ")
    args = parser.parse_args()

    # 出力ディレクトリを設定
    if args.output_dir:
        OUTPUT_DIR = Path(args.output_dir)
    else:
        OUTPUT_DIR = DEFAULT_OUTPUT_DIR

    print("=" * 60)
    print("NCBI 配列取得スクリプト（改良版 v2）")
    print("=" * 60)
    print("新機能:")
    print("  - 複数アクセッション番号の確実な取得")
    print("  - 異常配列の自動除外（品質フィルタリング）")
    print(f"出力先: {OUTPUT_DIR}")
    print("=" * 60)

    # 出力ディレクトリ確認
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 古いファイルを削除（重要：古いデータの混入を防ぐ）
    print("\n--- 古いファイルをクリーンアップ ---")
    clean_old_files()

    all_sequences = {}
    target_species = args.target.strip()
    target_gene = args.gene.strip()

    # 種名のスペルチェック
    suggestion = suggest_species_name(target_species)
    if suggestion and suggestion.lower() != target_species.lower():
        print(f"\n⚠ スペル確認: '{target_species}' → '{suggestion}' ではありませんか？")

    # 種名がNCBIに存在するか確認
    print(f"\n--- 種名確認: {target_species} ---")
    exists = verify_species_exists(target_species)
    if exists is False:
        print(f"  ⚠ 警告: '{target_species}' がNCBI Taxonomyに見つかりません。")
        print(f"  → スペルを確認してください。")
        if suggestion:
            print(f"  → もしかして: {suggestion}")
    elif exists is True:
        print(f"  ✓ 種名確認OK")

    # Target
    print(f"\n--- 標的種: {target_species} ({target_gene}) ---")
    seqs = search_and_fetch_sequences(
        target_species,
        target_gene,
        max_results=args.max_seqs,
        min_sequences=args.min_seqs
    )
    if seqs:
        all_sequences["target"] = seqs
        print_sequence_summary(seqs)
        save_sequences(seqs, "target_sequence.fasta")
    else:
        # 修正候補がある場合は自動的に試す
        if suggestion and suggestion.lower() != target_species.lower():
            print(f"\n  → スペル修正版で再検索: {suggestion}")
            seqs = search_and_fetch_sequences(
                suggestion,
                target_gene,
                max_results=args.max_seqs,
                min_sequences=args.min_seqs
            )
            if seqs:
                all_sequences["target"] = seqs
                print_sequence_summary(seqs)
                save_sequences(seqs, "target_sequence.fasta")
                target_species = suggestion  # 修正後の名前を使用

    if "target" not in all_sequences:
        print(f"\n⚠ 警告: 標的種の配列が見つかりませんでした。")
        print("  試したこと:")
        print(f"    - {target_gene} の複数のエイリアスで検索")
        print(f"    - mitochondrion/mitochondrial で検索")
        print("\n  対策:")
        print(f"    1. 種名のスペルを確認: {target_species}")
        print(f"    2. 別の遺伝子を試す: 16S, 18S, ITS, cytb")
        print(f"    3. NCBIで手動検索して配列をダウンロード")
        print(f"    4. sequences/target_sequence.fasta に手動で配置")

    # Related species
    if args.related.strip():
        related_list = [s.strip() for s in args.related.split(",")]
        related_seqs = []
        for sp in related_list:
            if not sp: continue
            print(f"\n--- 近縁種: {sp} ---")

            # スペル修正を試す
            sp_suggestion = suggest_species_name(sp)
            if sp_suggestion and sp_suggestion.lower() != sp.lower():
                print(f"  (スペル修正: {sp} → {sp_suggestion})")
                sp = sp_suggestion

            r_seqs = search_and_fetch_sequences(sp, target_gene, max_results=5, min_sequences=2)
            if r_seqs:
                print_sequence_summary(r_seqs)
            related_seqs.extend(r_seqs)
        if related_seqs:
            all_sequences["related"] = related_seqs
            save_sequences(related_seqs, "related_sequences.fasta")

    # Non-target species
    if args.nontarget.strip():
        nontarget_list = [s.strip() for s in args.nontarget.split(",")]
        nontarget_seqs = []
        for sp in nontarget_list:
            if not sp: continue
            print(f"\n--- 非標的生物: {sp} ---")
            nt_seqs = search_and_fetch_sequences(sp, target_gene, max_results=3, min_sequences=1)
            if nt_seqs:
                print_sequence_summary(nt_seqs)
            nontarget_seqs.extend(nt_seqs)
        if nontarget_seqs:
            all_sequences["nontarget"] = nontarget_seqs
            save_sequences(nontarget_seqs, "nontarget_sequences.fasta")

    # 全配列を1ファイルにまとめる
    all_seqs = []
    for group, seqs in all_sequences.items():
        all_seqs.extend(seqs)

    if all_seqs:
        save_sequences(all_seqs, f"all_{target_gene}_sequences.fasta")

    print("\n" + "=" * 60)
    print("配列取得完了")
    print("=" * 60)

    # サマリー
    print("\n取得配列サマリー:")
    for group, seqs in all_sequences.items():
        accessions = [s.id for s in seqs]
        print(f"  {group}: {len(seqs)}配列")
        print(f"    アクセッション: {', '.join(accessions)}")

    print(f"\n合計: {len(all_seqs)}配列")
    print(f"保存先: {OUTPUT_DIR}")

    # 標的配列が見つからなかった場合は警告を強調
    if "target" not in all_sequences:
        print("\n" + "!" * 60)
        print("重要: 標的種の配列が見つかりませんでした！")
        print("後続のステップは実行できません。")
        print("!" * 60)
        sys.exit(1)

    print("\n✓ 配列取得成功！次のステップに進めます。")


if __name__ == "__main__":
    main()
