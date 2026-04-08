#!/usr/bin/env python3
"""
auto_visualize.py — 種名 (またはプロジェクトディレクトリ名) を指定するだけで、
最も推奨されるプライマーに対して target / nontarget / all / excluded の各 FASTA を
DNA Dynamo .cow ファイルへ自動アラインメント可視化する。

前提:
    既に run_full_workflow.py で reports/<Species>_<Gene>/ が生成済みであること
    （sequences/, primer_candidates/primer_candidates.csv が存在する）。

使い方:
    python3 scripts/auto_visualize.py --project Helicoverpa_armigera_COI
    python3 scripts/auto_visualize.py --project Helicoverpa_armigera_COI --top-n 3
    python3 scripts/auto_visualize.py --project Helicoverpa_armigera_COI --no-automation

オプション:
    --project NAME         reports/<NAME>/ ディレクトリを処理
    --top-n N              上位 N プライマー (Quality_Score 降順) を処理 (default: 1)
    --categories CAT...    target / nontarget / all / excluded から選ぶ (default: 全部)
    --no-automation        cow 生成のみ、DNA Dynamo は開かない
    --skip-existing        既に存在する .cow はスキップ
"""

import argparse
import csv
import subprocess
import sys
import time
from pathlib import Path


PKG_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = PKG_ROOT / "scripts"
REPORTS_DIR = PKG_ROOT / "reports"
WORKFLOW_SCRIPT = SCRIPTS_DIR / "generate_primer_cow_workflow.py"

# 出力フォルダ名（reports/<project>/<OUTPUT_FOLDER>/...）
OUTPUT_FOLDER = "プライマー結合"

# カテゴリキー → 日本語ファイル名
JAPANESE_NAMES = {
    "target": "ターゲット配列のみ",
    "nontarget": "非標的配列のみ",
    "nontarget_strict": "非標的配列のみ_strict",
    "all": "すべての配列",
}


def japanese_filename(category_key: str) -> str:
    """カテゴリキーを日本語ファイル名（拡張子なし）に変換"""
    if category_key in JAPANESE_NAMES:
        return JAPANESE_NAMES[category_key]
    if category_key.startswith("excluded-"):
        return "除外した種_" + category_key[len("excluded-"):]
    return category_key


def species_prefix_from_project(project_name: str) -> str:
    """プロジェクト名から短いプレフィックスを生成 (例: Helicoverpa_armigera_COI → Ha-COI)"""
    parts = project_name.split("_")
    if len(parts) >= 3:
        return parts[0][0].upper() + parts[1][0].lower() + "-" + parts[2]
    return project_name


def pick_top_primers(csv_path: Path, top_n: int) -> list[dict]:
    """primer_candidates.csv を Quality_Score 降順で並べて上位 N 件を返す"""
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    rows.sort(key=lambda r: float(r.get("Quality_Score", 0) or 0), reverse=True)
    return rows[:top_n]


def discover_fasta_categories(sequences_dir: Path, requested: list[str]) -> dict[str, Path]:
    """sequences ディレクトリから処理対象の FASTA を見つける

    Returns: { "target": Path, "nontarget": Path, "all": Path, "excluded_<name>": Path, ... }
    """
    found: dict[str, Path] = {}

    if "target" in requested:
        p = sequences_dir / "target_sequence.fasta"
        if p.exists():
            found["target"] = p

    if "nontarget" in requested:
        # related_sequences.fasta が設計時の "nontarget" 候補
        p = sequences_dir / "related_sequences.fasta"
        if p.exists():
            found["nontarget"] = p
        # 旧形式の nontarget_sequences.fasta もあれば併用
        p2 = sequences_dir / "nontarget_sequences.fasta"
        if p2.exists():
            found["nontarget_strict"] = p2

    if "all" in requested:
        # all_<gene>_sequences.fasta などを探す
        candidates = sorted(sequences_dir.glob("all_*sequences.fasta"))
        if candidates:
            found["all"] = candidates[0]

    if "excluded" in requested:
        for p in sorted(sequences_dir.glob("excluded_*.fasta")):
            # excluded_Spodoptera_litura_COI.fasta → key: excluded_Spodoptera_litura
            key = p.stem
            # 短くする: excluded_Spodoptera_litura_COI → excluded-Spodoptera_litura
            key = key.replace("excluded_", "excluded-")
            found[key] = p

    return found


def category_label(category_key: str) -> str:
    """カテゴリキーから --input-label の値を決定"""
    if category_key == "target":
        return "TARGET"
    if category_key.startswith("nontarget"):
        return "NONTARGET"
    if category_key == "all":
        return "ALL"
    if category_key.startswith("excluded"):
        return "EXCLUDED"
    return "TARGET"


def run_one(
    primer_csv: Path,
    primer_id: str,
    primer_name: str,
    input_fasta: Path,
    input_label: str,
    output_cow: Path,
    no_automation: bool,
    auto_quit: bool,
) -> int:
    """generate_primer_cow_workflow.py を 1 回呼び出す"""
    cmd = [
        sys.executable,
        str(WORKFLOW_SCRIPT),
        "--primers-csv", str(primer_csv),
        "--primer-id", primer_id,
        "--input-fasta", str(input_fasta),
        "--input-label", input_label,
        "--output", str(output_cow),
        "--primer-name", primer_name,
    ]
    if no_automation:
        cmd.append("--no-automation")
    if auto_quit:
        cmd.append("--auto-quit")
    print("\n" + "=" * 70)
    print(f">>> {primer_name}  ({input_label}: {input_fasta.name})")
    print("=" * 70)
    return subprocess.run(cmd).returncode


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--project", "-P", required=True,
                        help="reports/ 配下のプロジェクト名 (例: Helicoverpa_armigera_COI)")
    parser.add_argument("--top-n", type=int, default=1,
                        help="上位 N プライマーを処理 (default: 1)")
    parser.add_argument(
        "--categories", nargs="+",
        default=["target", "nontarget", "all", "excluded"],
        choices=["target", "nontarget", "all", "excluded"],
        help="処理する FASTA カテゴリ (default: 全部)",
    )
    parser.add_argument("--no-automation", action="store_true",
                        help="cow 生成のみで DNA Dynamo は起動しない")
    parser.add_argument("--skip-existing", action="store_true",
                        help="既に存在する .cow はスキップ")
    args = parser.parse_args()

    # 入力サニタイズ: --project は許可文字のみ
    import re
    if not re.match(r"^[A-Za-z0-9._-]+$", args.project) or ".." in args.project:
        print(f"エラー: --project には英数字 / アンダースコア / ハイフン / ドットのみ使用可能 (..は不可): {args.project}")
        sys.exit(1)

    project_dir = (REPORTS_DIR / args.project).resolve()
    # 防御深度: 解決後のパスが必ず REPORTS_DIR 配下にあることを確認
    try:
        project_dir.relative_to(REPORTS_DIR.resolve())
    except ValueError:
        print(f"エラー: 解決されたプロジェクトパスが reports/ の外にあります: {project_dir}")
        sys.exit(1)
    if not project_dir.exists():
        print(f"エラー: プロジェクトが見つかりません: {project_dir}")
        print(f"  → 先に auto_full.py で {args.project} を作成してください")
        sys.exit(1)

    csv_path = project_dir / "primer_candidates" / "primer_candidates.csv"
    if not csv_path.exists():
        print(f"エラー: primer_candidates.csv が見つかりません: {csv_path}")
        sys.exit(1)

    sequences_dir = project_dir / "sequences"
    if not sequences_dir.exists():
        print(f"エラー: sequences ディレクトリが見つかりません: {sequences_dir}")
        sys.exit(1)

    output_dir = project_dir / OUTPUT_FOLDER
    output_dir.mkdir(parents=True, exist_ok=True)

    top_primers = pick_top_primers(csv_path, args.top_n)
    if not top_primers:
        print("エラー: プライマー候補が空です")
        sys.exit(1)

    fasta_map = discover_fasta_categories(sequences_dir, args.categories)
    if not fasta_map:
        print(f"エラー: 処理対象の FASTA が見つかりません ({args.categories})")
        sys.exit(1)

    prefix = species_prefix_from_project(args.project)

    print("=" * 70)
    print("auto_visualize: マルチカテゴリ DNA Dynamo 自動可視化")
    print("=" * 70)
    print(f"  Project       : {args.project}")
    print(f"  Prefix        : {prefix}")
    print(f"  Top primers   : {args.top_n}  (Quality_Score 降順)")
    for p in top_primers:
        print(f"    - ID={p['ID']:3} F={p['Forward']} R={p['Reverse']} Q={p.get('Quality_Score','?')}")
    print(f"  FASTA targets :")
    for k, v in fasta_map.items():
        print(f"    - {k:25} → {v.name}")
    print(f"  Automation    : {'OFF (no-automation)' if args.no_automation else 'ON (DNA Dynamo)'}")
    print()

    total = len(top_primers) * len(fasta_map)
    done = 0
    failed = []
    skipped = []
    succeeded = []

    multi_primer = len(top_primers) > 1

    for primer in top_primers:
        pid = primer["ID"]
        # top-N (N>1) の場合は primer ごとにサブフォルダを作って衝突回避
        if multi_primer:
            primer_subdir = output_dir / f"プライマー候補{pid}"
            primer_subdir.mkdir(parents=True, exist_ok=True)
            current_dir = primer_subdir
        else:
            current_dir = output_dir

        for category, fasta_path in fasta_map.items():
            done += 1
            jp_name = japanese_filename(category)
            # AppleScript automation 用の primer-name (DNA Dynamo ウィンドウタイトルに使われる)
            display_name = f"{prefix}-{pid}-{category}" if multi_primer else f"{prefix}-{pid}"
            output_cow = current_dir / f"{jp_name}.cow"
            print(f"\n[{done}/{total}]")

            if args.skip_existing and output_cow.exists():
                print(f"  skip (既存): {output_cow}")
                skipped.append(jp_name)
                continue

            label = category_label(category)
            rc = run_one(
                primer_csv=csv_path,
                primer_id=pid,
                primer_name=display_name,
                input_fasta=fasta_path,
                input_label=label,
                output_cow=output_cow,
                no_automation=args.no_automation,
                auto_quit=not args.no_automation,
            )
            if rc != 0:
                print(f"  ✗ 失敗: {jp_name}")
                failed.append(jp_name)
            else:
                print(f"  ✓ 成功: {jp_name}")
                succeeded.append(jp_name)
                # cow ファイル以外の sidecar (.fasta, .gb) を削除
                for ext in (".fasta", ".gb"):
                    sidecar = output_cow.with_suffix(ext)
                    if sidecar.exists():
                        try:
                            sidecar.unlink()
                        except OSError:
                            pass
                # automation 後は次の cow との衝突回避で少し待つ
                if not args.no_automation:
                    time.sleep(2)

    print("\n" + "=" * 70)
    print("auto_visualize: サマリー")
    print("=" * 70)
    print(f"  成功 : {len(succeeded)}/{total}")
    for n in succeeded:
        print(f"    ✓ {n}")
    if skipped:
        print(f"  スキップ : {len(skipped)}")
        for n in skipped:
            print(f"    - {n}")
    if failed:
        print(f"  失敗 : {len(failed)}")
        for n in failed:
            print(f"    ✗ {n}")
    print()
    print(f"  生成先 : {output_dir}")
    sys.exit(0 if not failed else 1)


if __name__ == "__main__":
    main()
