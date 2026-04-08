#!/usr/bin/env python3
"""
auto_full.py — 種名を入れるだけで NCBI 配列取得 → Primer3 設計 → クロスチェック →
レポート生成 → 最推奨プライマーの DNA Dynamo .cow 可視化までを 1 コマンドで完遂する。

使い方:
    python3 scripts/auto_full.py --target "Helicoverpa armigera" --gene COI
    python3 scripts/auto_full.py --target "Bemisia tabaci" --gene COI --related "Trialeurodes vaporariorum,Aleyrodes proletella"
    python3 scripts/auto_full.py --target "Helicoverpa armigera" --gene COI --at-rich
    python3 scripts/auto_full.py --target "Helicoverpa armigera" --gene COI --skip-design
    python3 scripts/auto_full.py --target "Helicoverpa armigera" --gene COI --no-automation

ステップ:
    1. fetch_sequences.py     NCBI から配列取得
    2. design_primers_primer3.py  Primer3 でプライマー設計
    3. primer_crosscheck.py   標的・近縁種クロスチェック
    4. generate_full_report.py  総合レポート生成
    5. auto_visualize.py      最推奨プライマー (top-1) で 4 カテゴリ cow 生成
"""

import argparse
import subprocess
import sys
from pathlib import Path


PKG_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = PKG_ROOT / "scripts"
REPORTS_DIR = PKG_ROOT / "reports"


def run_step(step_no: int, label: str, cmd: list[str]) -> None:
    print()
    print("=" * 70)
    print(f"STEP {step_no}: {label}")
    print("=" * 70)
    r = subprocess.run(cmd)
    if r.returncode != 0:
        print(f"\n✗ STEP {step_no} 失敗 ({label})")
        sys.exit(r.returncode)


def run_design_with_retry(step_no: int, base_cmd: list[str], at_rich: bool) -> bool:
    """design_primers_primer3.py を実行。失敗時 (AT-rich 未指定なら) 自動で再試行する。

    Returns: 最終的に at_rich モードを使ったかどうか
    """
    label = "Primer3 でプライマー設計 (design_primers_primer3.py)"
    print()
    print("=" * 70)
    print(f"STEP {step_no}: {label}")
    print("=" * 70)

    cmd = list(base_cmd)
    if at_rich:
        cmd.append("--at-rich")

    r = subprocess.run(cmd, capture_output=True, text=True)
    sys.stdout.write(r.stdout)
    if r.stderr:
        sys.stderr.write(r.stderr)

    failed = r.returncode != 0 or "プライマー候補が見つかりませんでした" in r.stdout

    if failed and not at_rich:
        print()
        print("⚠ Primer3 が候補を見つけられませんでした。")
        print("  → AT-rich 配列 (Lepidoptera 等の mtDNA) の可能性があるため、")
        print("    --at-rich モードで自動的に再試行します。")
        print()
        print("=" * 70)
        print(f"STEP {step_no} (retry): Primer3 設計 + --at-rich")
        print("=" * 70)
        retry_cmd = list(base_cmd) + ["--at-rich"]
        r2 = subprocess.run(retry_cmd, capture_output=True, text=True)
        sys.stdout.write(r2.stdout)
        if r2.stderr:
            sys.stderr.write(r2.stderr)
        failed2 = r2.returncode != 0 or "プライマー候補が見つかりませんでした" in r2.stdout
        if failed2:
            print(f"\n✗ STEP {step_no} 失敗 (--at-rich 再試行も不可)")
            sys.exit(r2.returncode if r2.returncode != 0 else 1)
        return True

    if failed:
        print(f"\n✗ STEP {step_no} 失敗 ({label})")
        sys.exit(r.returncode if r.returncode != 0 else 1)

    return at_rich


def safe_name(s: str) -> str:
    """ファイル/ディレクトリ名として安全な文字列に変換する。

    セキュリティ要件:
        - パストラバーサル ('..' / '/' / '\\') を防ぐ
        - 隠しファイル化 (先頭 '.') を防ぐ
        - NUL バイトや制御文字を除去
        - 許可文字: 英数字, アンダースコア, ハイフン, ドット
        - 半角空白は '_' に変換
    """
    import re
    s = (s or "").strip().replace(" ", "_")
    s = re.sub(r"[^A-Za-z0-9._-]", "_", s)
    s = s.lstrip(".")  # 先頭ドット (隠しファイル化) を防ぐ
    s = s.replace("..", "_")  # 念のため二重ドット除去
    return s or "unknown"


def project_folder_name(target: str, gene: str) -> str:
    return f"{safe_name(target)}_{safe_name(gene)}"


def interactive_input(args) -> None:
    """run_full_workflow.py を踏襲したインタラクティブ入力。
    --target が未指定なら対話モードに入る。"""
    print("=" * 70)
    print("auto_full インタラクティブ設定")
    print("=" * 70)
    print()
    print("配列取得とプライマー設計の条件を入力してください。")
    print("(空欄可の項目は何も入力せず Enter で次へ)")
    print()

    if not args.target:
        while True:
            v = input("標的種を入力してください (例: Helicoverpa armigera): ").strip()
            if v:
                args.target = v
                break
            print("  ※ 標的種は必須です")

    if not args.gene:
        while True:
            v = input("対象遺伝子を入力してください (例: COI, matK): ").strip()
            if v:
                args.gene = v
                break
            print("  ※ 遺伝子は必須です")

    if not args.related:
        v = input("近縁種をカンマ区切りで入力 [空欄可]: ").strip()
        if v:
            args.related = v

    if not args.nontarget:
        v = input("非標的生物をカンマ区切りで入力 [空欄可]: ").strip()
        if v:
            args.nontarget = v

    if not args.at_rich:
        v = input("AT-rich モードを使用しますか？ (Lepidoptera 等の mtDNA) [y/N]: ").strip().lower()
        args.at_rich = v in ("y", "yes")

    print()
    print("【入力確認】")
    print(f"  標的種       : {args.target}")
    print(f"  遺伝子       : {args.gene}")
    print(f"  近縁種       : {args.related or '(なし)'}")
    print(f"  非標的生物   : {args.nontarget or '(なし)'}")
    print(f"  AT-rich      : {'ON' if args.at_rich else 'OFF (失敗時は自動フォールバック)'}")
    print()
    input("Enter で開始 (Ctrl+C で中止)...")
    print()


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--target", "-t", default="",
                   help='標的種 (例: "Helicoverpa armigera"、未指定なら対話入力)')
    p.add_argument("--gene", "-g", default="",
                   help="対象遺伝子 (例: COI, matK、未指定なら対話入力)")
    p.add_argument("--related", default="",
                   help="近縁種 (カンマ区切り、空欄可)")
    p.add_argument("--nontarget", default="",
                   help="非標的生物 (カンマ区切り、空欄可)")
    p.add_argument("--at-rich", action="store_true",
                   help="AT-rich モード (Lepidoptera 等の mtDNA 向け)")
    p.add_argument("--max-seqs", type=int, default=10,
                   help="種ごとの最大取得配列数 (default: 10)")
    p.add_argument("--skip-design", action="store_true",
                   help="既に設計済みなら fetch〜report をスキップして可視化のみ")
    p.add_argument("--no-automation", action="store_true",
                   help="cow 生成のみで DNA Dynamo は起動しない")
    p.add_argument("--categories", nargs="+",
                   default=["target", "nontarget", "all", "excluded"],
                   choices=["target", "nontarget", "all", "excluded"],
                   help="可視化する FASTA カテゴリ (default: 全部)")
    p.add_argument("--non-interactive", action="store_true",
                   help="対話入力をスキップ (CLI 引数で全部指定する場合)")
    args = p.parse_args()

    # --target または --gene が未指定なら対話モード
    if not args.non_interactive and (not args.target or not args.gene):
        interactive_input(args)

    project_name = project_folder_name(args.target, args.gene)
    project_dir = (REPORTS_DIR / project_name).resolve()

    # 防御深度: 解決後のパスが必ず REPORTS_DIR 配下にあることを確認
    try:
        project_dir.relative_to(REPORTS_DIR.resolve())
    except ValueError:
        print(f"エラー: 解決されたプロジェクトパスが reports/ の外にあります: {project_dir}")
        print("  → 種名・遺伝子名に '..' などのパストラバーサル文字が含まれていないか確認してください")
        sys.exit(1)

    print("=" * 70)
    print("auto_full: 種名 → DNA Dynamo cow 完全自動化")
    print("=" * 70)
    print(f"  Target       : {args.target}")
    print(f"  Gene         : {args.gene}")
    print(f"  Related      : {args.related or '(なし)'}")
    print(f"  Nontarget    : {args.nontarget or '(なし)'}")
    print(f"  AT-rich      : {'ON' if args.at_rich else 'OFF'}")
    print(f"  Project dir  : {project_dir}")
    print(f"  Skip design  : {'YES' if args.skip_design else 'NO'}")
    print(f"  Automation   : {'OFF' if args.no_automation else 'ON (DNA Dynamo)'}")
    print()

    if not args.skip_design:
        # サブフォルダを作成
        for sub in ("sequences", "primer_candidates", "crosscheck_results", "dnadynamo_files"):
            (project_dir / sub).mkdir(parents=True, exist_ok=True)

        # STEP 1: fetch
        fetch_cmd = [
            sys.executable, str(SCRIPTS_DIR / "fetch_sequences.py"),
            "--target", args.target,
            "--gene", args.gene,
            "--max-seqs", str(args.max_seqs),
            "--output-dir", str(project_dir / "sequences"),
        ]
        if args.related:
            fetch_cmd += ["--related", args.related]
        if args.nontarget:
            fetch_cmd += ["--nontarget", args.nontarget]
        run_step(1, "NCBI 配列取得 (fetch_sequences.py)", fetch_cmd)

        # STEP 2: design (失敗時は --at-rich で自動再試行)
        design_cmd = [
            sys.executable, str(SCRIPTS_DIR / "design_primers_primer3.py"),
            "--input-dir", str(project_dir / "sequences"),
            "--output-dir", str(project_dir / "primer_candidates"),
        ]
        used_at_rich = run_design_with_retry(2, design_cmd, at_rich=args.at_rich)
        if used_at_rich and not args.at_rich:
            print(f"  ℹ 自動 fallback: --at-rich モードで設計が成功しました")

        # STEP 3: crosscheck
        cross_cmd = [
            sys.executable, str(SCRIPTS_DIR / "primer_crosscheck.py"),
            "--sequences-dir", str(project_dir / "sequences"),
            "--primers-dir", str(project_dir / "primer_candidates"),
            "--output-dir", str(project_dir / "crosscheck_results"),
        ]
        run_step(3, "クロスチェック (primer_crosscheck.py)", cross_cmd)

        # STEP 4: report
        report_cmd = [
            sys.executable, str(SCRIPTS_DIR / "generate_full_report.py"),
            "--target", args.target,
            "--gene", args.gene,
            "--project-dir", str(project_dir),
        ]
        run_step(4, "総合レポート生成 (generate_full_report.py)", report_cmd)
    else:
        if not project_dir.exists():
            print(f"エラー: --skip-design 指定だがプロジェクトが存在しません: {project_dir}")
            sys.exit(1)
        print(f"--skip-design 指定: {project_dir} の既存データを使用します")

    # STEP 5: auto_visualize (top-1 のみ)
    visualize_cmd = [
        sys.executable, str(SCRIPTS_DIR / "auto_visualize.py"),
        "--project", project_name,
        "--top-n", "1",
        "--categories", *args.categories,
    ]
    if args.no_automation:
        visualize_cmd.append("--no-automation")
    run_step(5, "最推奨プライマーで DNA Dynamo cow 生成 (auto_visualize.py)", visualize_cmd)

    print()
    print("=" * 70)
    print("✓ auto_full 完全成功")
    print("=" * 70)
    print(f"  プロジェクト : {project_dir}")
    print(f"  cow ファイル : {project_dir / 'dnadynamo_files'}")
    print()


if __name__ == "__main__":
    main()
