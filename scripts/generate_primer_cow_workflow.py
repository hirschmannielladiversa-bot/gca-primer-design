#!/usr/bin/env python3
"""
DNA Dynamo .cow 自動アラインメントワークフロー

generate_primer_cow.py で .cow + .fasta を生成し、続けて DNA Dynamo を起動して
AppleScript で以下を自動実行する:

  1. .fasta を一度開いて recent files に登録（FASTA ビューアダイアログを閉じる）
  2. .cow を開く
  3. Alignments メニュー → "Align this sequence to other open, recent, or filesystem files"
  4. Drag Drop Assembly Window で fasta が選択された状態にして "Assemble Sequences"
  5. "Select Sequences to Align" ダイアログで "Align Selected Sequences"
  6. → コンセンサス + プライマー feature + 24 配列のアラインメントビューが完成

使用例:
    python generate_primer_cow_workflow.py \
        --primers-csv reports/Helicoverpa_armigera_COI/primer_candidates/primer_candidates.csv \
        --primer-id 3 \
        --sequences-dir reports/Helicoverpa_armigera_COI/sequences \
        --include-related \
        --output reports/Helicoverpa_armigera_COI/dnadynamo_files/Ha-COI-3.cow \
        --primer-name Ha-COI-3
"""

import argparse
import subprocess
import time
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.absolute()
GENERATE_SCRIPT = SCRIPT_DIR / "generate_primer_cow.py"

# DNA Dynamo の System Events 上のプロセス名
# DNA Dynamo は Java アプリなので JavaAppLauncher として登録される
DNADYNAMO_PROCESS = "JavaAppLauncher"
# DNA Dynamo の起動コマンド (open -a の引数)
DNADYNAMO_APP = "DNADynamo"


APPLESCRIPT_TEMPLATE = '''
on run argv
    set cowPath to item 1 of argv
    set fastaPath to item 2 of argv
    set fastaFileName to item 3 of argv
    set cowBaseName to item 4 of argv
    set clickHelperPath to item 5 of argv
    set pythonExePath to item 6 of argv
    set processName to "{process}"
    set appName to "{app}"
    set assemblyWinName to "DNADynamo Sequencing - Drag Drop Assembly Window"

    -- ステップ 0: DNA Dynamo を事前起動して完全に立ち上がるまで待つ
    -- 既に起動済みなら no-op、未起動ならスプラッシュを通過するまで待機
    do shell script "open -a " & quoted form of appName
    set t to 0
    repeat until (my procExists(processName)) or t > 60
        delay 0.5
        set t to t + 1
    end repeat
    if not (my procExists(processName)) then
        error "DNA Dynamo process did not start"
    end if
    -- スプラッシュ画面が消えてメイン UI が応答可能になるまでさらに待機
    delay 4

    -- ステップ 1: FASTA を recent files に登録するため open する
    -- 開く前のウィンドウ数を記録 → open → ウィンドウ数増加を待つ → Cmd+W で閉じ → 元の数に戻ったか確認
    set baseWinCount to my winCount(processName)
    do shell script "open -a " & quoted form of appName & " " & quoted form of fastaPath

    set t to 0
    repeat until (my winCount(processName) > baseWinCount) or t > 40
        delay 0.3
        set t to t + 1
    end repeat
    delay 1.0

    -- FASTA ビューアウィンドウを閉じる
    -- このウィンドウは Cmd+W を受け付けない Java JDialog なので、
    -- タイトルバーの close button (AXCloseButton) を直接クリックする
    tell application "System Events"
        tell process processName
            set frontmost to true
            delay 0.5
            -- 戦略 1: Escape
            key code 53
            delay 0.4
            -- 戦略 2: Cmd+W
            keystroke "w" using command down
            delay 0.4
            -- 戦略 3: cow ウィンドウ以外のウィンドウの close button をクリック
            try
                repeat with w in (every window)
                    set wn to name of w
                    if wn does not contain "primer alignment" and wn does not contain cowBaseName then
                        try
                            -- AXCloseButton をクリック
                            click (first button of w whose subrole is "AXCloseButton")
                        on error
                            try
                                -- フォールバック: button 1
                                click button 1 of w
                            end try
                        end try
                    end if
                end repeat
            end try
            delay 0.5
        end tell
    end tell

    -- ダイアログが実際に閉じるまで待つ
    set t to 0
    repeat until (my winCount(processName) ≤ baseWinCount) or t > 20
        delay 0.3
        set t to t + 1
    end repeat
    delay 0.5

    -- ステップ 2: .cow を開く
    set baseWinCount2 to my winCount(processName)
    do shell script "open -a " & quoted form of appName & " " & quoted form of cowPath

    -- cow ウィンドウ（名前に primer alignment or cowBaseName が含まれる）が出現するまで待機
    set t to 0
    set cowReady to false
    repeat until cowReady or t > 60
        delay 0.5
        set t to t + 1
        try
            tell application "System Events"
                tell process processName
                    repeat with w in (every window)
                        set wn to name of w
                        if wn contains "primer alignment" or wn contains cowBaseName then
                            set cowReady to true
                            exit repeat
                        end if
                    end repeat
                end tell
            end tell
        end try
    end repeat
    if not cowReady then
        error "cow window did not appear within timeout"
    end if
    delay 1.0

    -- ステップ 3 前準備: DNA Dynamo を確実にアクティブ化する
    -- Java Swing は activate イベントで menu bar を完成させるため、
    -- "tell application ... activate" を使う必要がある
    try
        tell application appName to activate
    end try
    delay 1.0

    -- ステップ 3-7: メニュー操作 → アラインメント実行 → 保存
    tell application "System Events"
        tell process processName
            set frontmost to true
            delay 0.8

            -- cow ウィンドウを明示的に最前面に持ってくる
            set raised to false
            try
                repeat with w in (every window)
                    set wn to name of w
                    if wn contains "primer alignment" or wn contains cowBaseName then
                        perform action "AXRaise" of w
                        set raised to true
                        exit repeat
                    end if
                end repeat
            end try
            if not raised then
                error "could not raise cow window"
            end if
            delay 1.0

            -- もう一度アクティブ化（AXRaise でフォーカスが戻る可能性に対する保険）
            set frontmost to true
            delay 0.5

            -- ステップ 3: cow ウィンドウ内のツールバーから "Align To" ボタンをクリック
            -- 理由: Java accessibility は menu items を公開しないが、
            -- Java Swing の Window 内の button は AX tree に出るためクリック可能。
            -- ツールバーの "Align To" ボタンは Alignments メニューと同等の機能を持つ。
            set cowWin to missing value
            try
                repeat with w in (every window)
                    set wn to name of w
                    if wn contains "primer alignment" or wn contains cowBaseName then
                        set cowWin to w
                        exit repeat
                    end if
                end repeat
            end try
            if cowWin is missing value then
                error "cow window not found for button click"
            end if

            -- entire contents of で再帰的に button を探す
            set targetBtn to missing value
            set foundBtnNames to ""
            try
                set allElems to entire contents of cowWin
                repeat with elem in allElems
                    try
                        set elemClass to class of elem
                        if elemClass is button then
                            set bn to ""
                            try
                                set bn to name of elem
                            end try
                            if bn is not "" then
                                set foundBtnNames to foundBtnNames & bn & ", "
                                if bn is "Align To" then
                                    set targetBtn to elem
                                end if
                            end if
                        end if
                    end try
                end repeat
            on error errMsg
                error "entire contents enumeration failed: " & errMsg
            end try

            if targetBtn is missing value then
                error "Align To button not found. Available buttons: [" & foundBtnNames & "]"
            end if

            click targetBtn
            delay 1.5

            -- ステップ 4: Drag Drop Assembly Window 出現待ち
            set t to 0
            repeat until (exists window assemblyWinName) or t > 50
                delay 0.3
                set t to t + 1
            end repeat
            if not (exists window assemblyWinName) then
                error "Drag Drop Assembly Window did not appear"
            end if
            delay 0.5

            -- ステップ 5: ファイルリストから fasta 行を見つけてチェックボックスを ON
            -- 階層パスは Java Swing の更新で変わる可能性があるため、
            -- entire contents of で全要素を列挙し、static text "Ha-COI-3.fasta" の
            -- Y 座標と最も近い checkbox を選ぶ（位置ベースマッチング）。
            set assemblyWin to window assemblyWinName
            set allElems to entire contents of assemblyWin

            -- 1) Ha-COI-3.fasta の static text を探して Y 座標取得
            set fastaTextY to -1
            repeat with elem in allElems
                try
                    if class of elem is static text then
                        set elemName to ""
                        try
                            set elemName to name of elem
                        end try
                        if elemName is fastaFileName then
                            set tp to position of elem
                            set fastaTextY to (item 2 of tp)
                            exit repeat
                        end if
                    end if
                end try
            end repeat
            if fastaTextY is -1 then
                error "static text '" & fastaFileName & "' not found in assembly window"
            end if

            -- 2) Y 座標が最も近い checkbox を探す
            -- ただし上部のフィルタ checkbox (All, DNA, abi/scf trace, other, folders) は除外。
            -- これらはウィンドウ上部にあり、Y 座標が小さい。
            -- ファイルリスト領域はウィンドウ中央。fasta の Y ± 10 px 以内の checkbox を選ぶ。
            set bestCheckbox to missing value
            set bestDist to 100000
            repeat with elem in allElems
                try
                    if class of elem is checkbox then
                        set cbPos to position of elem
                        set cbY to (item 2 of cbPos)
                        set d to cbY - fastaTextY
                        if d < 0 then set d to -d
                        if d < bestDist then
                            set bestDist to d
                            set bestCheckbox to elem
                        end if
                    end if
                end try
            end repeat
            if bestCheckbox is missing value or bestDist > 30 then
                -- 診断: すべての checkbox の位置と value をダンプ
                set cbDump to ""
                repeat with elem in allElems
                    try
                        if class of elem is checkbox then
                            set cbP to position of elem
                            set cbV to 0
                            try
                                set cbV to value of elem
                            end try
                            set cbDump to cbDump & "Y=" & (item 2 of cbP) & " val=" & cbV & "; "
                        end if
                    end try
                end repeat
                error "no checkbox close to fasta text Y=" & fastaTextY & " (best dist=" & bestDist & ") | checkboxes: " & cbDump
            end if
            -- Java Swing 対策: size=347x16 の行全体が checkbox として公開されているため、
            -- click の中心座標はテキスト領域に落ちてチェックボックス本体をヒットしない。
            -- Space キーも意図しないウィンドウにフォーカスが飛ぶリスクがある。
            -- 唯一確実なのは「checkbox graphic が実際にある右端座標を Quartz で直接クリック」する方法。
            set cbPos to position of bestCheckbox
            set cbSize to size of bestCheckbox
            set clickX to (item 1 of cbPos) + (item 1 of cbSize) - 8
            set clickY to (item 2 of cbPos) + ((item 2 of cbSize) / 2)

            -- Quartz (pyobjc) による座標クリック
            -- セキュリティ: Python パスは Python 側で事前解決され argv 経由で渡される。
            -- これにより AppleScript 内でのパスハードコードを排除し、PATH 操作攻撃を防ぐ。
            -- 全引数 (path / 座標) を quoted form of で囲んでシェル注入を防ぐ。
            set pyQ to quoted form of pythonExePath
            set helperQ to quoted form of clickHelperPath
            set xQ to quoted form of (clickX as string)
            set yQ to quoted form of (clickY as string)
            try
                do shell script pyQ & " " & helperQ & " " & xQ & " " & yQ
            on error errMsg
                error "Quartz click failed: " & errMsg
            end try
            delay 0.8

            -- クリック後に value が 1 になっているか検証
            set cbValue to 0
            try
                set cbValue to value of bestCheckbox
            end try
            if cbValue is not 1 then
                error "checkbox click did not toggle. Quartz click at (" & clickX & "," & clickY & ") pos=" & (item 1 of cbPos) & "," & (item 2 of cbPos) & " size=" & (item 1 of cbSize) & "," & (item 2 of cbSize) & " final value=" & cbValue
            end if
            delay 0.3

            -- ステップ 6: "Assemble Sequences" ボタンを entire contents から探してクリック
            set assembleBtn to missing value
            repeat with elem in allElems
                try
                    if class of elem is button then
                        set bn to ""
                        try
                            set bn to name of elem
                        end try
                        if bn is "Assemble Sequences" then
                            set assembleBtn to elem
                            exit repeat
                        end if
                    end if
                end try
            end repeat
            if assembleBtn is missing value then
                error "Assemble Sequences button not found"
            end if
            click assembleBtn
            delay 1.5

            -- ステップ 7: "Select Sequences to Align" ダイアログ → Align Selected Sequences
            set t to 0
            repeat until (exists window "Select Sequences to Align") or t > 30
                delay 0.3
                set t to t + 1
            end repeat
            if exists window "Select Sequences to Align" then
                click button "Align Selected Sequences" of window "Select Sequences to Align"
                delay 1
            end if

            -- ステップ 8: アラインメント完了待ち（Drag Drop Assembly Window が消えるまで）
            set t to 0
            repeat while (exists window assemblyWinName) and t < 60
                delay 0.5
                set t to t + 1
            end repeat
            delay 1.5

            -- ステップ 9: Cmd+S で保存
            keystroke "s" using command down
            delay 1.5

            -- 保存ダイアログが出た場合は Enter で確定
            try
                if (count of windows whose subrole is "AXDialog") > 0 then
                    keystroke return
                    delay 1
                end if
            end try
        end tell
    end tell

    return "OK"
end run

on procExists(pn)
    tell application "System Events"
        return (exists process pn)
    end tell
end procExists

on winCount(pn)
    try
        tell application "System Events"
            tell process pn
                return count of windows
            end tell
        end tell
    on error
        return 0
    end try
end winCount
'''


def run_generate(args):
    """generate_primer_cow.py を呼び出す"""
    cmd = [sys.executable, str(GENERATE_SCRIPT)]
    if args.primers_csv:
        cmd += ["--primers-csv", args.primers_csv]
    if args.primer_id:
        cmd += ["--primer-id", args.primer_id]
    if args.forward:
        cmd += ["--forward", args.forward]
    if args.reverse:
        cmd += ["--reverse", args.reverse]
    if args.include_related:
        cmd += ["--include-related"]
    if args.input_fasta:
        cmd += ["--input-fasta", args.input_fasta]
        cmd += ["--input-label", args.input_label]
    elif args.sequences_dir:
        cmd += ["--sequences-dir", args.sequences_dir]
    cmd += [
        "--output", args.output,
        "--primer-name", args.primer_name,
    ]
    if args.template:
        cmd += ["--template", args.template]
    if args.flank is not None:
        cmd += ["--flank", str(args.flank)]
    print("===== ファイル生成 =====")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print("エラー: ファイル生成に失敗しました")
        sys.exit(1)


def resolve_python_with_quartz() -> str | None:
    """Quartz (pyobjc) が import できる python3 の絶対パスを返す。

    セキュリティ: Python 側で事前解決することで、AppleScript 内の
    `do shell script "/opt/homebrew/bin/python3 ..."` のようなパスハードコードを避け、
    PATH 操作攻撃の余地をなくす。
    """
    candidates = [
        sys.executable,
        "/opt/homebrew/bin/python3",
        "/usr/local/bin/python3",
        "/Library/Frameworks/Python.framework/Versions/3.13/bin/python3",
        "/Library/Frameworks/Python.framework/Versions/3.12/bin/python3",
        "/Library/Frameworks/Python.framework/Versions/3.11/bin/python3",
        "/usr/bin/python3",
    ]
    seen = set()
    for p in candidates:
        if not p or p in seen:
            continue
        seen.add(p)
        if not Path(p).exists():
            continue
        try:
            r = subprocess.run(
                [p, "-c", "import Quartz"],
                capture_output=True, timeout=5,
            )
            if r.returncode == 0:
                return p
        except (subprocess.TimeoutExpired, OSError):
            continue
    return None


def quit_dnadynamo():
    """DNA Dynamo を完全終了して次のワークフロー実行のためにクリーンな状態にする"""
    print("  DNA Dynamo を終了中...")
    subprocess.run(
        ["osascript", "-e", 'tell application "DNADynamo" to quit'],
        capture_output=True,
    )
    # プロセスが消えるまで最大 8 秒待つ
    for _ in range(16):
        r = subprocess.run(["pgrep", "-f", "DNADynamo"], capture_output=True)
        if r.returncode != 0:
            break
        time.sleep(0.5)
    time.sleep(1.0)


def run_applescript(cow_path: Path, fasta_path: Path):
    """AppleScript を実行して DNA Dynamo を自動操作"""
    script = APPLESCRIPT_TEMPLATE.format(
        process=DNADYNAMO_PROCESS,
        app=DNADYNAMO_APP,
    )
    print("\n===== DNA Dynamo 自動操作 =====")
    print(f"  cow:   {cow_path.name}")
    print(f"  fasta: {fasta_path.name}")
    print("  AppleScript 実行中（DNA Dynamo を操作するため、画面に触れないでください）...")
    click_helper = Path(__file__).parent / "_click_at.py"

    # セキュリティ: AppleScript 内で /opt/homebrew/bin/python3 をハードコードすると
    # PATH 操作攻撃の余地が生じる。Python 側で Quartz 利用可能な python3 を解決して
    # 絶対パスを AppleScript に argv 経由で渡す。
    python_path = resolve_python_with_quartz()
    if python_path is None:
        print("\n  ✗ Quartz (pyobjc) がインストールされた python3 が見つかりません")
        print("  対処方法: pip3 install pyobjc-framework-Quartz")
        sys.exit(1)

    result = subprocess.run(
        ["osascript", "-e", script,
         str(cow_path), str(fasta_path), fasta_path.name, cow_path.stem,
         str(click_helper), python_path],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print("\n  ✗ AppleScript エラー:")
        print("  " + result.stderr.replace("\n", "\n  "))
        print("\n  対処方法:")
        print("  - DNADynamo がインストールされているか確認")
        print("  - 「システム設定 > プライバシーとセキュリティ > アクセシビリティ」で")
        print("    Terminal (または python) に許可を与える")
        print("  - DNADynamo のメニュー名・ボタン名が変わっていないか確認")
        sys.exit(1)
    print(f"  ✓ 完了: {result.stdout.strip()}")


def main():
    parser = argparse.ArgumentParser(
        description="DNA Dynamo .cow 完全自動アラインメントワークフロー",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--sequences-dir", "-s",
                        help="target_sequence.fasta などを含むディレクトリ")
    parser.add_argument("--input-fasta",
                        help="単一 FASTA を入力にするモード（--sequences-dir の代替）")
    parser.add_argument("--input-label", default="TARGET",
                        choices=["TARGET", "NONTARGET", "RELATED", "ALL", "EXCLUDED"])
    parser.add_argument("--output", "-o", required=True, help="出力 .cow パス")
    parser.add_argument("--primer-name", "-n", default="Primer")
    parser.add_argument("--primers-csv", "-p")
    parser.add_argument("--primer-id")
    parser.add_argument("--forward")
    parser.add_argument("--reverse")
    parser.add_argument("--template")
    parser.add_argument("--flank", type=int)
    parser.add_argument("--include-related", action="store_true")
    parser.add_argument("--no-automation", action="store_true",
                        help="ファイル生成のみ実行（DNA Dynamo は起動しない）")
    parser.add_argument("--auto-quit", action="store_true",
                        help="開始前に DNA Dynamo を quit してクリーンな状態にする")

    args = parser.parse_args()

    if not args.input_fasta and not args.sequences_dir:
        print("エラー: --sequences-dir または --input-fasta のいずれかを指定してください")
        sys.exit(1)

    # ステップ 0: クリーンな状態にする
    if args.auto_quit and not args.no_automation:
        quit_dnadynamo()

    # ステップ 1: ファイル生成
    run_generate(args)

    if args.no_automation:
        print("\n--no-automation 指定のためここで終了します。")
        return

    # ステップ 2: AppleScript で DNA Dynamo を自動操作
    cow_path = Path(args.output)
    fasta_path = cow_path.with_suffix(".fasta")
    if not cow_path.exists() or not fasta_path.exists():
        print(f"エラー: 生成されたファイルが見つかりません ({cow_path}, {fasta_path})")
        sys.exit(1)

    run_applescript(cow_path, fasta_path)

    print("\n===== 完了 =====")
    print("DNA Dynamo にアラインメントウィンドウが表示されているはずです。")


if __name__ == "__main__":
    main()
