#!/bin/zsh
# verify_slw_count.sh — 指定プロジェクトの「プライマー結合/」配下の .cow ファイルの
# SLW count を一覧表示
#
# 使い方: ./verify_slw_count.sh Helicoverpa_armigera_COI
#
# セキュリティ注意:
#   - 引数 PROJECT_NAME は許可文字 [A-Za-z0-9._-] のみ受け付ける
#   - Python ソース内に変数を直接展開しない (環境変数経由で渡す)

set -e
cd "$(dirname "$0")"

PROJECT="${1:-}"
if [ -z "$PROJECT" ]; then
    echo "使い方: $0 <PROJECT_NAME>"
    echo "例:    $0 Helicoverpa_armigera_COI"
    exit 1
fi

# 入力サニタイズ: 許可文字のみ
case "$PROJECT" in
    *[!A-Za-z0-9._-]*|.*|*..*)
        echo "エラー: PROJECT_NAME には英数字 / アンダースコア / ドット / ハイフンのみ使用可能 (..は不可)"
        exit 1
        ;;
esac

DIR="reports/$PROJECT/プライマー結合"

if [ ! -d "$DIR" ]; then
    echo "エラー: $DIR が見つかりません"
    exit 1
fi

PYTHON=/opt/homebrew/bin/python3
[ ! -x "$PYTHON" ] && PYTHON=/Library/Frameworks/Python.framework/Versions/3.13/bin/python3
[ ! -x "$PYTHON" ] && PYTHON=python3

# 環境変数で値を渡し、ヒアドキュメント内での変数展開を完全に抑止 ('PYEOF')
PROJECT_NAME="$PROJECT" DIR_PATH="$DIR" "$PYTHON" - <<'PYEOF'
import os
import glob
import zipfile

try:
    import javaobj
except ImportError:
    print("Error: javaobj-py3 が必要です: pip install javaobj-py3")
    raise SystemExit(1)

project = os.environ["PROJECT_NAME"]
dir_path = os.environ["DIR_PATH"]
files = sorted(glob.glob(os.path.join(dir_path, "*.cow")))
if not files:
    print("(.cow ファイルなし)")
else:
    print(f"プロジェクト: {project}")
    print("=" * 70)
    for p in files:
        try:
            spw = javaobj.loads(zipfile.ZipFile(p).read("SPW"))
            n = len(spw.slwArray)
            size = os.path.getsize(p)
            print(f"  {n:4} SLW  {size:8} bytes  {os.path.basename(p)}")
        except Exception as e:
            print(f"  ERROR  {os.path.basename(p)}: {e}")
PYEOF
