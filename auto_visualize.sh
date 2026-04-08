#!/bin/zsh
# auto_visualize.sh — プロジェクト名を引数にして DNA Dynamo マルチカテゴリ可視化を実行
#
# 使い方:
#   ./auto_visualize.sh Helicoverpa_armigera_COI                    # top1, 全カテゴリ, 自動操作 ON
#   ./auto_visualize.sh Helicoverpa_armigera_COI 3                  # top3
#   ./auto_visualize.sh Helicoverpa_armigera_COI 1 --no-automation  # cow 生成のみ
#
# 前提: reports/<NAME>/sequences/, primer_candidates/primer_candidates.csv が存在
set -e
cd "$(dirname "$0")"

if [ -z "$1" ]; then
    echo "使い方: $0 <PROJECT_NAME> [TOP_N] [追加オプション...]"
    echo "例:    $0 Helicoverpa_armigera_COI"
    echo "       $0 Helicoverpa_armigera_COI 3"
    echo "       $0 Helicoverpa_armigera_COI 1 --no-automation"
    exit 1
fi

PROJECT="$1"
shift

TOP_N=1
if [ -n "$1" ] && [[ "$1" =~ ^[0-9]+$ ]]; then
    TOP_N="$1"
    shift
fi

python3 scripts/auto_visualize.py --project "$PROJECT" --top-n "$TOP_N" "$@"
