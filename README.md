# GCA Primer Auto Design

**種名を入れるだけで、NCBI 配列取得 → Primer3 設計 → クロスチェック → DNA Dynamo 可視化までを完全自動化する macOS 向けツールキット**

GCA (Gut Content Analysis: 消化管内容物 DNA 解析) 用のプライマーを設計する研究者のために作られました。

---

## 特徴

- **ワンコマンド** — `./auto_full.sh` で対話モードが起動し、種名・遺伝子名を入力するだけ
- **AT-rich 自動対応** — 鱗翅目・双翅目の AT-rich mtDNA でも Primer3 失敗時に `--at-rich` で自動再試行
- **4 カテゴリ可視化** — ターゲット / 非標的 / 全配列 / 除外種 を別々の DNA Dynamo `.cow` ファイルに自動アラインメント
- **macOS 完全自動操作** — Quartz API 経由で DNA Dynamo の Drag Drop Assembly Window を自動操作 (オプション)
- **クロス汚染防止** — 各種が `reports/<Species>_<Gene>/` の独立フォルダに収まる

---

## クイックスタート

```bash
git clone https://github.com/TKG-M/gca-primer-design.git
cd gca-primer-design

# 依存パッケージ
pip3 install -r requirements.txt

# (可視化を使う場合のみ)
pip3 install pyobjc-framework-Quartz

# 実行権限
chmod +x auto_full.sh auto_visualize.sh verify_slw_count.sh

# 一発実行
./auto_full.sh
```

対話入力:

```
標的種を入力してください (例: Helicoverpa armigera): Helicoverpa armigera
対象遺伝子を入力してください (例: COI, matK): COI
近縁種をカンマ区切りで入力 [空欄可]: Helicoverpa assulta,Spodoptera exigua
非標的生物をカンマ区切りで入力 [空欄可]: Spodoptera litura
AT-rich モードを使用しますか？ [y/N]: y
```

完了後:

```
reports/Helicoverpa_armigera_COI/プライマー結合/
├── ターゲット配列のみ.cow
├── 非標的配列のみ.cow
├── すべての配列.cow
└── 除外した種_Spodoptera_litura.cow
```

---

## ドキュメント

詳細な使い方・トラブルシューティング・Q&A は **[manual.html](manual.html)** をブラウザで開いてください。

---

## 動作環境

- **macOS** 10.15 以降 (DNA Dynamo 自動操作は macOS 専用)
- **Python** 3.10 以降
- **DNA Dynamo** (可視化部分のみ。ファイル生成だけなら不要)

Linux / Windows でも配列取得・設計・cow ファイル生成までは動きますが、cow を DNA Dynamo で表示する自動操作部分は使えません (`--no-automation` を付けて実行してください)。

---

## ライセンス

MIT License — 詳細は [LICENSE](LICENSE) を参照。

外部ライブラリ:
- [Biopython](https://biopython.org/) — Biopython License
- [primer3-py](https://github.com/libnano/primer3-py) — GPL v2
- [javaobj-py3](https://github.com/tcalmant/python-javaobj) — Apache 2.0
- [pyobjc](https://pyobjc.readthedocs.io/) — MIT

---

## 著者

**TKG.M**

バグ報告・機能要望は GitHub Issues にお願いします。
