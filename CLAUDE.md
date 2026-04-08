# CLAUDE.md - GCA Primer Auto Design Package

**プロジェクト**: GCA (Gut Content Analysis) プライマー自動設計パッケージ
**最終更新**: 2026-04-08
**バージョン**: v1.5

---

## プロジェクト概要

NCBI から標的種の配列を取得し、Primer3 で GCA 用プライマーを自動設計し、結果を DNA Dynamo `.cow` ファイルとして可視化するまでをワンコマンドで完結するワークフロー。

### 主要機能

1. **配列取得** (`fetch_sequences.py`)
   - NCBI から標的種・近縁種・非標的生物の配列を取得
   - 複数クエリパターン (`COI`, `CO1`, `cox1`, ...) で網羅検索
   - k-mer 類似度による異常配列フィルタリング

2. **プライマー設計** (`design_primers_primer3.py`)
   - Primer3 を使用した GCA 最適化設計
   - 産物サイズ 100–160 bp (分解 DNA 対応)
   - `--at-rich` で Lepidoptera 等の AT-rich mtDNA 対応

3. **クロスチェック** (`primer_crosscheck.py`)
   - 近縁種・非標的生物との交差反応評価
   - ミスマッチ数によるリスク判定

4. **レポート生成** (`generate_full_report.py`)
   - 包括的 Markdown レポート

5. **DNA Dynamo 可視化** (`generate_primer_cow.py` + `generate_primer_cow_workflow.py`)
   - プライマー × 配列のアラインメントを `.cow` 形式で生成
   - 4 カテゴリ (target / nontarget / all / excluded) を一括出力
   - macOS 上で DNA Dynamo を AppleScript で自動操作 (オプション)

6. **完全自動化** (`auto_full.py` / `auto_full.sh`)
   - 種名を入力するだけで配列取得から `.cow` 生成まで一発実行
   - 対話モード / CLI モード両対応
   - Primer3 失敗時に `--at-rich` で自動フォールバック

---

## フォルダ構造

```
gca_auto_design_package/
├── CLAUDE.md                       # このファイル
├── README.md                       # ユーザー向け概要
├── USAGE.txt                       # 使い方早見表
├── requirements.txt
├── auto_full.sh                    # 完全自動化エントリポイント (推奨)
├── auto_visualize.sh               # 可視化のみエントリポイント
├── verify_slw_count.sh             # SLW 検証ヘルパー
├── run_full_workflow.py            # 旧インタラクティブワークフロー
├── scripts/
│   ├── auto_full.py                # 完全自動化オーケストレータ
│   ├── auto_visualize.py           # マルチカテゴリ可視化
│   ├── fetch_sequences.py
│   ├── design_primers_primer3.py
│   ├── primer_crosscheck.py
│   ├── generate_full_report.py
│   ├── generate_dnadynamo_batch.py
│   ├── generate_primer_cow.py      # cow ファイルビルダー
│   ├── generate_primer_cow_workflow.py  # AppleScript 自動操作付き
│   └── _click_at.py                # Quartz 座標クリック ヘルパー
├── templates/                      # cow テンプレート
└── reports/
    └── {Species}_{Gene}/
        ├── sequences/              # FASTA (target/related/all/excluded)
        ├── primer_candidates/      # Primer3 出力
        ├── crosscheck_results/     # クロスチェック
        ├── report.md
        └── プライマー結合/         # cow ファイル
            ├── ターゲット配列のみ.cow
            ├── 非標的配列のみ.cow
            ├── すべての配列.cow
            └── 除外した種_*.cow
```

---

## 使用方法

### 1. 完全自動化 (最推奨)

```bash
./auto_full.sh
```

→ 標的種・遺伝子・近縁種・非標的・AT-rich を順次対話入力 → NCBI 取得 → Primer3 設計 → クロスチェック → レポート → 4 カテゴリ `.cow` 生成 → DNA Dynamo 自動操作 まで一発完遂。

または引数指定:

```bash
./auto_full.sh "Helicoverpa armigera" COI
./auto_full.sh "Helicoverpa armigera" COI --at-rich
./auto_full.sh "Bemisia tabaci" COI --related "Trialeurodes vaporariorum,Aleyrodes proletella"
./auto_full.sh "Helicoverpa armigera" COI --skip-design     # 既存設計の cow だけ更新
./auto_full.sh "Helicoverpa armigera" COI --no-automation   # cow 生成のみ DNA Dynamo 起動なし
```

### 2. 可視化のみ (設計済みのプロジェクトに対して)

```bash
./auto_visualize.sh Helicoverpa_armigera_COI
```

### 3. SLW 数の検証

```bash
./verify_slw_count.sh Helicoverpa_armigera_COI
```

---

## AT-rich モードについて

Lepidoptera (鱗翅目)、Diptera (双翅目) などの**ミトコンドリア DNA は AT 含有率が 70% を超える**ことが多く、Primer3 のデフォルトパラメータでは候補が見つからないことがあります。

`--at-rich` を指定すると以下が緩和されます:

- `PRIMER_MIN_TM`: 50°C
- `PRIMER_MAX_TM`: 60°C
- `PRIMER_MIN_GC`: 20%
- `PRIMER_MAX_GC`: 60%

`auto_full.py` は **Primer3 失敗時に `--at-rich` で自動フォールバック**します。判断不要です。

---

## DNA Dynamo 自動操作について

`generate_primer_cow_workflow.py` は AppleScript + Quartz で DNA Dynamo を自動操作し、`.cow` ファイルにアラインメント結果を埋め込みます。

### 環境要件

- macOS
- DNA Dynamo がインストール済み
- システム設定 → プライバシーとセキュリティ → アクセシビリティ で Terminal (または Python) を許可
- `/opt/homebrew/bin/python3` または `/Library/Frameworks/Python.framework/Versions/*/bin/python3` に Quartz (pyobjc) が入っていること

### 仕組み

1. DNA Dynamo を起動 → cow ファイルを開く
2. Java accessibility では menu items が公開されないため、ツールバーの `Align To` ボタンを AX tree から探してクリック
3. Drag Drop Assembly Window で `Ha-COI-3.fasta` の static text を Y 座標一致で特定
4. **Quartz (pyobjc) で行右端をピクセル座標クリック** ← ここが決定打
5. `Assemble Sequences` → `Align Selected Sequences` → 完了待機 → Cmd+S 保存

詳細な技術的経緯は git ログを参照。

---

## 配列フィルタリング

### 類似度閾値

- `MIN_SIMILARITY_THRESHOLD = 0.15` : 絶対最小
- `OUTLIER_THRESHOLD = 0.20` : グローバル平均からの許容下限
- `MIN_SEQUENCES_TO_KEEP = 3` : 最低保持配列数

### 出力ファイル

- `target_sequence.fasta` : 標的種フィルタリング後
- `related_sequences.fasta` : 近縁種
- `all_{Gene}_sequences.fasta` : 全配列
- `excluded_{Species}_{Gene}.fasta` : 除外配列

---

## 依存パッケージ

```bash
pip install -r requirements.txt
# = biopython, primer3-py, javaobj-py3
```

DNA Dynamo 自動操作のみ追加で:
```bash
# pyobjc がない場合
pip install pyobjc-framework-Quartz
```

---

## 注意事項

1. **NCBI Entrez API**: 利用前に `Entrez.email` を自分のメールアドレスに設定すること
2. **プライマー検証**: in-silico 設計のため、実験的検証が必要
3. **DNA Dynamo 自動操作**: 実行中は画面に触れない (キーボード入力やウィンドウ切替で失敗する)
4. **AT-rich 自動フォールバック**: 鱗翅目などで Primer3 が候補を見つけられなくても自動再試行されるが、それでもダメな場合は `--max-seqs` を増やすか配列を見直す

---

## 関連ドキュメント

- `manual.html` : ユーザーマニュアル (Q&A 付き)
- `README.md` : 概要
- `USAGE.txt` : 使い方早見表
- `~/.claude/rules/primer-design-validation.md` : プライマー検証プロトコル
- `~/.claude/rules/work-integrity-standards.md` : 作業品質基準

---

## 変更履歴

### v1.5 (2026-04-08)
- `auto_full.py` / `auto_full.sh` 追加: 種名を入力するだけで設計から DNA Dynamo cow まで一発実行
- インタラクティブ入力モード (run_full_workflow 互換)
- Primer3 失敗時に `--at-rich` で自動フォールバック
- 出力フォルダを `dnadynamo_files/` から `プライマー結合/` に変更
- 日本語ファイル名: `ターゲット配列のみ.cow` `非標的配列のみ.cow` `すべての配列.cow` `除外した種_*.cow`
- `.fasta`, `.gb` sidecar ファイルは自動削除 (`.cow` のみ保持)

### v1.4 (2026-04-07)
- `generate_primer_cow.py` の単一 FASTA モード (`--input-fasta`) 追加
- `auto_visualize.py` 追加: 4 カテゴリ FASTA を一括処理
- `generate_primer_cow_workflow.py` の AppleScript 自動操作完成
  - DNA Dynamo の Java accessibility 制限を Quartz 座標クリックで回避
  - Drag Drop Assembly Window のチェックボックスを Y 座標マッチング → 行右端 Quartz クリック
  - DNA Dynamo 自動 quit/relaunch (`--auto-quit`)

### v1.3 (2026-03-24)
- DNAdynamo 用 GenBank ファイル一括生成 (`generate_dnadynamo_batch.py`)
- 統合ワークフローに DNAdynamo 生成ステップ追加

### v1.2 (2026-03-24)
- AT-rich 配列用パラメータ (`--at-rich`) を追加 (Lepidoptera 等の mtDNA 対応)

### v1.1 (2026-03-24)
- DNA Dynamo ファイル生成機能を追加、フォルダ構造を更新

---

*Maintained by TKG.M*
