#!/bin/zsh
# auto_full.sh — 種名を入れるだけで配列取得→設計→cow ファイル生成まで一発で実行
#
# インタラクティブ実行 (引数なし、または引数不足):
#   ./auto_full.sh
#       → 標的種 / 遺伝子 / 近縁種 / 非標的 / AT-rich を順次入力
#
# 直接実行 (CLI 引数指定):
#   ./auto_full.sh "Helicoverpa armigera" COI
#   ./auto_full.sh "Helicoverpa armigera" COI --at-rich
#   ./auto_full.sh "Bemisia tabaci" COI --related "Trialeurodes vaporariorum,Aleyrodes proletella"
#   ./auto_full.sh "Helicoverpa armigera" COI --skip-design
#   ./auto_full.sh "Helicoverpa armigera" COI --no-automation

set -e
cd "$(dirname "$0")"

if [ -z "$1" ]; then
    # 完全インタラクティブ
    python3 scripts/auto_full.py
elif [ -z "$2" ]; then
    # 種名のみ → 遺伝子だけ追加で聞く
    python3 scripts/auto_full.py --target "$1"
else
    TARGET="$1"
    GENE="$2"
    shift 2
    python3 scripts/auto_full.py --target "$TARGET" --gene "$GENE" "$@"
fi
