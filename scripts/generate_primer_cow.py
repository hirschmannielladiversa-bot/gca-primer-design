#!/usr/bin/env python3
"""
DNA Dynamo .cow ファイル生成スクリプト

NCBIから取得した標的種・非標的種のアクセッション配列をプライマー結合位置で
アラインメントし、DNA Dynamo の Clustal mode で開ける .cow ファイルを生成する。

特徴:
- 各配列はプライマーの 5' 端で位置揃えされる（forward 開始位置で縦に揃う）
- 視覚的にミスマッチが分かる（標的種は一致、非標的種は不一致箇所が色分けで目立つ）
- 上部にコンセンサス配列、Forward / Reverse プライマーは feature として表示

依存:
    pip install biopython javaobj-py3

使用例:
    python generate_primer_cow.py \
        --primers-csv reports/Helicoverpa_armigera_COI/primer_candidates/primer_candidates.csv \
        --primer-id 3 \
        --sequences-dir reports/Helicoverpa_armigera_COI/sequences \
        --output reports/Helicoverpa_armigera_COI/dnadynamo_files/Ha-COI-3.cow \
        --primer-name Ha-COI-3

テンプレート:
    既存の .cow ファイル（デフォルト: モモアカアブラムシ.cow）を構造テンプレートとして
    読み込み、内部の SPW (Java serialized SequencingProjectWrapper)、xml、history を
    新しいデータで差し替える。
"""

import argparse
import copy
import csv
import hashlib
import io
import re
import sys
import zipfile
from collections import Counter
from datetime import datetime
from pathlib import Path

try:
    import javaobj
except ImportError:
    print("Error: javaobj-py3 が必要です: pip install javaobj-py3")
    sys.exit(1)

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    print("Error: Biopython が必要です: pip install biopython")
    sys.exit(1)


# ===== デフォルトテンプレート =====
# パッケージ同梱の templates/template_consensus.cow を使う
DEFAULT_TEMPLATE = str(
    Path(__file__).resolve().parent.parent / "templates" / "template_consensus.cow"
)


# ===== IUPAC ミスマッチ許容マッチング =====
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


def base_match(primer_base: str, template_base: str) -> bool:
    if template_base not in "ACGT":
        return False
    return template_base in IUPAC.get(primer_base, primer_base)


def find_best_binding(template: str, primer: str):
    """テンプレート上で primer の最良結合位置を全長スライディングで探索"""
    template = template.upper()
    primer = primer.upper()
    n, m = len(template), len(primer)
    if m == 0 or n < m:
        return None
    best = None
    for i in range(n - m + 1):
        mismatches = 0
        for j in range(m):
            if not base_match(primer[j], template[i + j]):
                mismatches += 1
                if best is not None and mismatches > best[0]:
                    break
        if best is None or mismatches < best[0]:
            best = (mismatches, i)
            if best[0] == 0:
                break
    if best is None:
        return None
    mc, start = best
    return (start, start + m, mc)


# ===== 配列ウィンドウ抽出とアラインメント =====
def extract_window(record, forward, reverse, flank=30):
    """
    配列から forward...reverse のアンプリコン領域（+ flanking）を切り出す。
    Returns dict or None
    """
    template = str(record.seq).upper()
    f_hit = find_best_binding(template, forward)
    if f_hit is None:
        return None
    f_start, f_end, f_mc = f_hit

    rev_rc = str(Seq(reverse.upper()).reverse_complement())
    r_hit = find_best_binding(template, rev_rc)

    if r_hit is not None:
        r_start, r_end, r_mc = r_hit
    else:
        # reverse が見つからない場合はおおよそ 200bp 後ろまで取る
        r_end = min(len(template), f_start + 250)
        r_start = max(f_start, r_end - len(reverse))
        r_mc = -1

    win_start = max(0, f_start - flank)
    win_end = min(len(template), r_end + flank)
    window = template[win_start:win_end]

    return {
        "name": record.id,
        "description": record.description,
        "full_sequence": template,           # 元 FASTA のフル配列
        "window": window,
        "fw_offset_in_window": f_start - win_start,
        "rv_end_in_window": r_end - win_start,
        "fw_mismatches": f_mc,
        "rv_mismatches": r_mc,
        "amplicon_size": r_end - f_start,
    }


def build_alignment(records_with_label, forward, reverse, flank=30):
    """
    複数配列を Forward プライマー開始位置で揃えた MSA を構築する。

    全配列の gappedSequence を同じ幅にする。
    """
    extracted = []
    skipped = []
    for rec, label in records_with_label:
        info = extract_window(rec, forward, reverse, flank)
        if info is None:
            skipped.append((rec.id, label))
            continue
        info["label"] = label
        extracted.append(info)

    if not extracted:
        return None, skipped

    max_before = max(e["fw_offset_in_window"] for e in extracted)
    max_after = max(len(e["window"]) - e["fw_offset_in_window"] for e in extracted)
    total_width = max_before + max_after

    for e in extracted:
        before_pad = max_before - e["fw_offset_in_window"]
        after_pad = total_width - before_pad - len(e["window"])
        e["gapped"] = "." * before_pad + e["window"] + "." * after_pad
        e["aligned_fw_start"] = max_before  # 0-based
        e["aligned_fw_end"] = max_before + len(forward)
        e["aligned_rv_end"] = before_pad + e["rv_end_in_window"]
        e["aligned_rv_start"] = e["aligned_rv_end"] - len(reverse)

    # コンセンサスは TARGET 行のみから計算（標的種のシグナルを保つ）
    target_rows = [e["gapped"] for e in extracted if e["label"] == "TARGET"]
    if not target_rows:
        target_rows = [e["gapped"] for e in extracted]
    consensus = compute_consensus(target_rows)

    return {
        "rows": extracted,
        "consensus": consensus,
        "width": total_width,
        "fw_start": max_before,                       # 0-based
        "fw_end": max_before + len(forward),          # 0-based exclusive
    }, skipped


def compute_consensus(rows):
    if not rows:
        return ""
    width = len(rows[0])
    out = []
    for i in range(width):
        col = [r[i] for r in rows if r[i] not in (".", "-", "N")]
        if not col:
            out.append("N")
        else:
            base, _ = Counter(col).most_common(1)[0]
            out.append(base)
    return "".join(out)


def sanitize_name(name: str) -> str:
    """DNA Dynamo の表示用に空白などを置換"""
    return re.sub(r"[\s/]+", "_", name)


def md5_hex(s: str) -> str:
    return hashlib.md5(s.encode("utf-8")).hexdigest()


# ===== SPW (Java シリアライズ) 構築 =====
def build_spw(template_spw, alignment, primer_name: str):
    """
    テンプレート SPW を「コンセンサスのみ」に変更する。

    DNA Dynamo の Java 直列化のバックリファレンス制約により、slwArray の一部を
    変更すると Sequencing Project として認識されなくなる。そのためデータ行は
    SPW に入れず、代わりに同名の .fasta を別途出力してドラッグ＆ドロップで
    取り込んでもらう運用とする。
    """
    if len(template_spw.slwArray) < 1:
        raise ValueError("テンプレート .cow が空です")
    template_consensus_slw = copy.deepcopy(template_spw.slwArray[0])
    template_spw.slwArray.clear()

    cs = copy.deepcopy(template_consensus_slw)
    cs.gappedSequence = alignment["consensus"]
    cs.name = None
    cs.sequence = None
    cs.accession = None
    cs.annotationXML = ""
    template_spw.slwArray.append(cs)

    template_spw.consensusMode = 3
    template_spw.showConsensusInClustalMode = False
    template_spw.lastEdited = int(datetime.now().timestamp() * 1000)
    return template_spw


def write_alignment_fasta(fasta_path: Path, alignment, primer_name: str):
    """
    DNA Dynamo にドラッグ＆ドロップ用の FASTA を書き出す。
    各配列はプライマー結合付近の窓を抽出したもの（gapped 表記の '.' は除去）。
    """
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    with open(fasta_path, "w") as f:
        for row in alignment["rows"]:
            ungapped = row["gapped"].replace(".", "")
            label = row["label"]
            f.write(f">[{label}]_{row['name']} mm_F={row['fw_mismatches']} mm_R={row['rv_mismatches']}\n")
            # 60bp 折返し
            for i in range(0, len(ungapped), 60):
                f.write(ungapped[i:i + 60] + "\n")


def write_alignment_genbank(gb_path: Path, alignment, primer_name: str,
                            forward: str, reverse: str, target_organism: str):
    """
    DNA Dynamo の "Open Sequences In Alignment Editor" 用マルチレコード GenBank。

    レコード1: コンセンサス（プライマー primer_bind feature 付き）
    レコード2+: 標的種・非標的種の各配列（feature なし）
    """
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    records = []
    consensus = alignment["consensus"]
    fw_start = alignment["fw_start"]                # 0-based
    fw_end = alignment["fw_end"]                    # 0-based exclusive
    # Reverse は中央値の位置を採る
    rv_starts = [r["aligned_rv_start"] for r in alignment["rows"]]
    rv_ends = [r["aligned_rv_end"] for r in alignment["rows"]]
    rv_start = sorted(rv_starts)[len(rv_starts) // 2]
    rv_end = sorted(rv_ends)[len(rv_ends) // 2]

    # === レコード 1: コンセンサス + プライマー feature ===
    cons_rec = SeqRecord(
        Seq(consensus),
        id=f"{primer_name}_consensus"[:16],
        name=f"{primer_name}_consensus"[:16],
        description=f"{target_organism} consensus for primer {primer_name}",
    )
    cons_rec.annotations["molecule_type"] = "DNA"
    cons_rec.annotations["topology"] = "linear"
    cons_rec.annotations["organism"] = target_organism
    cons_rec.annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()

    cons_rec.features = [
        SeqFeature(
            FeatureLocation(fw_start, fw_end, strand=1),
            type="primer_bind",
            qualifiers={
                "label": ["Forward"],
                "note": [forward],
                "PCR_conditions": [f"primer for {primer_name}"],
            },
        ),
        SeqFeature(
            FeatureLocation(rv_start, rv_end, strand=-1),
            type="primer_bind",
            qualifiers={
                "label": ["Reverse"],
                "note": [reverse],
                "PCR_conditions": [f"primer for {primer_name}"],
            },
        ),
        SeqFeature(
            FeatureLocation(fw_start, rv_end, strand=1),
            type="misc_feature",
            qualifiers={
                "label": ["PCR_product"],
                "note": [f"{rv_end - fw_start} bp amplicon"],
            },
        ),
    ]
    records.append(cons_rec)

    # === レコード 2+: 各配列（gapped から '.' を除去したもの） ===
    for row in alignment["rows"]:
        ungapped = row["gapped"].replace(".", "")
        rec_id = sanitize_name(f"{row['label']}_{row['name']}")[:16]
        rec = SeqRecord(
            Seq(ungapped),
            id=rec_id,
            name=rec_id,
            description=f"[{row['label']}] {row['name']} mm_F={row['fw_mismatches']} mm_R={row['rv_mismatches']}",
        )
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["topology"] = "linear"
        rec.annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()
        records.append(rec)

    gb_path.parent.mkdir(parents=True, exist_ok=True)
    with open(gb_path, "w") as f:
        SeqIO.write(records, f, "genbank")


# ===== xml 構築 =====
def build_annotation_entry(
    template_entry: str,
    label: str,
    start: int,         # 1-based
    end: int,           # 1-based inclusive
    feature_name: str,
    display_text: str,
    comments: str,
    color: int,
    is_complement: bool = False,
) -> str:
    """テンプレートの primer_bind entry をベースに新しい entry を作る"""
    e = template_entry
    e = re.sub(r"<label>.*?</label>", f"<label>{escape_xml(label)}</label>", e, count=1)
    e = re.sub(r"<startOfAnnotation>\d+</startOfAnnotation>",
               f"<startOfAnnotation>{start}</startOfAnnotation>", e)
    e = re.sub(r"<endOfAnnotation>\d+</endOfAnnotation>",
               f"<endOfAnnotation>{end}</endOfAnnotation>", e)
    e = re.sub(r"<comments>.*?</comments>",
               f"<comments>{escape_xml(comments)}</comments>", e, count=1, flags=re.DOTALL)
    e = re.sub(r"<annoColor>-?\d+</annoColor>",
               f"<annoColor>{color}</annoColor>", e)
    e = re.sub(r"<displayText>.*?</displayText>",
               f"<displayText>{escape_xml(display_text)}</displayText>", e, count=1)
    e = re.sub(r"<featureName>.*?</featureName>",
               f"<featureName>{escape_xml(feature_name)}</featureName>", e, count=1)
    # isComp タグは存在する場合と存在しない場合がある
    if is_complement:
        if "<isComp>" in e:
            e = re.sub(r"<isComp>(true|false)</isComp>", "<isComp>true</isComp>", e)
        else:
            e = e.replace("<direction>", "<isComp>true</isComp><direction>", 1)
    else:
        e = re.sub(r"<isComp>(true|false)</isComp>", "", e)
    return e


def escape_xml(s: str) -> str:
    return (s.replace("&", "&amp;")
             .replace("<", "&lt;")
             .replace(">", "&gt;"))


def build_xml(template_xml: str, alignment, primer_name: str, forward: str, reverse: str,
              accessions_used: list, target_organism: str) -> str:
    """テンプレート xml を流用して、配列・アノテーション・notes を新しいものに置換"""
    consensus = alignment["consensus"]
    consensus_disp = consensus.replace("N", "n")  # 表示用
    width = len(consensus)
    fw_start_1 = alignment["fw_start"] + 1                    # 1-based
    fw_end_1 = alignment["fw_end"]                            # 1-based inclusive
    # Reverse の表示位置：rows の中央値を採る（標的種で見える位置）
    rv_starts = [r["aligned_rv_start"] for r in alignment["rows"] if r["aligned_rv_start"] >= 0]
    rv_ends = [r["aligned_rv_end"] for r in alignment["rows"] if r["aligned_rv_end"] >= 0]
    if rv_starts and rv_ends:
        rv_start_1 = sorted(rv_starts)[len(rv_starts) // 2] + 1
        rv_end_1 = sorted(rv_ends)[len(rv_ends) // 2]
    else:
        rv_start_1 = max(1, width - len(reverse) + 1)
        rv_end_1 = width

    seq_md5 = md5_hex(consensus)
    accession_id = f"{primer_name}_alignment"

    # === <sequence1> 内の置換 ===
    new_xml = template_xml
    new_xml = re.sub(r"<seq>.*?</seq>",
                     f"<seq>{consensus}</seq>", new_xml, count=1)
    new_xml = re.sub(
        r"<definitionGB>.*?</definitionGB>",
        f"<definitionGB>DEFINITION  {escape_xml(target_organism)} primer alignment "
        f"({primer_name}) from {len(alignment['rows'])} accessions.</definitionGB>",
        new_xml, count=1, flags=re.DOTALL,
    )
    new_xml = re.sub(
        r"<sourceGB>.*?</sourceGB>",
        f"<sourceGB>SOURCE      {escape_xml(target_organism)}\n"
        f"  ORGANISM  {escape_xml(target_organism)}\n"
        f"            Unclassified.</sourceGB>",
        new_xml, count=1, flags=re.DOTALL,
    )
    new_xml = re.sub(
        r"<organismGB>.*?</organismGB>",
        f"<organismGB>  ORGANISM  {escape_xml(target_organism)}\n"
        f"            Unclassified.</organismGB>",
        new_xml, count=1, flags=re.DOTALL,
    )
    new_xml = re.sub(
        r"<accessionGB>.*?</accessionGB>",
        f"<accessionGB>ACCESSION   {accession_id}</accessionGB>",
        new_xml, count=1,
    )
    new_xml = re.sub(
        r"<seqHashGB>.*?</seqHashGB>",
        f"<seqHashGB>{seq_md5}</seqHashGB>",
        new_xml, count=1,
    )

    # === annotationBox 内のエントリを差し替え ===
    entries = re.findall(r"<entry>.*?</entry>", new_xml, re.DOTALL)
    # 既存の最初の entry (source) と primer_bind entry をテンプレートとして使う
    template_source = entries[0] if entries else None
    template_primer = None
    for e in entries:
        m = re.search(r"<featureName>(.*?)</featureName>", e)
        if m and m.group(1) == "primer_bind":
            template_primer = e
            break
    if template_primer is None:
        template_primer = entries[2] if len(entries) > 2 else template_source

    accs_str = ", ".join(accessions_used) if accessions_used else "(none)"
    new_entries = []

    # source feature（全長）
    if template_source:
        src_comments = (
            "Original Genbank Feature Record....\n"
            f"     source          1..{width}\n"
            f'                     /organism="{target_organism}"\n'
            f'                     /mol_type="genomic DNA"\n'
            f'                     /note="Primer alignment ({primer_name}) from '
            f'{len(alignment["rows"])} accessions"\n'
        )
        new_entries.append(build_annotation_entry(
            template_source,
            label=f"Imported from genbank file {accession_id}",
            start=1, end=width,
            feature_name="source",
            display_text=f"{primer_name} alignment ({len(alignment['rows'])} seqs)",
            comments=src_comments,
            color=-3618561,
        ))

    # PCR_product feature
    new_entries.append(build_annotation_entry(
        template_primer,
        label=f"Imported from genbank file {accession_id}",
        start=fw_start_1, end=rv_end_1,
        feature_name="misc_feature",
        display_text="PCR_product",
        comments=(
            "Original Genbank Feature Record....\n"
            f"     misc_feature    {fw_start_1}..{rv_end_1}\n"
            f'                     /label="PCR_product"\n'
            f'                     /note="Accessions used: {accs_str}"\n'
        ),
        color=-3618561,
    ))

    # Forward primer_bind
    new_entries.append(build_annotation_entry(
        template_primer,
        label=f"Imported from genbank file {accession_id}",
        start=fw_start_1, end=fw_end_1,
        feature_name="primer_bind",
        display_text="Forward",
        comments=(
            "Original Genbank Feature Record....\n"
            f"     primer_bind     {fw_start_1}..{fw_end_1}\n"
            f'                     /label="Forward"\n'
            f'                     /note="{forward}"\n'
            f'                     /PCR_conditions="primer for {primer_name}"\n'
        ),
        color=-52480,
    ))

    # Reverse primer_bind (complement)
    new_entries.append(build_annotation_entry(
        template_primer,
        label=f"Imported from genbank file {accession_id}",
        start=rv_start_1, end=rv_end_1,
        feature_name="primer_bind",
        display_text="Reverse",
        comments=(
            "Original Genbank Feature Record....\n"
            f"     primer_bind     complement({rv_start_1}..{rv_end_1})\n"
            f'                     /label="Reverse"\n'
            f'                     /note="{reverse}"\n'
            f'                     /PCR_conditions="primer for {primer_name}"\n'
        ),
        color=-52480,
        is_complement=True,
    ))

    new_annotation_box = "<annotationBox>" + "".join(new_entries) + "</annotationBox>"
    new_xml = re.sub(r"<annotationBox>.*?</annotationBox>",
                     new_annotation_box, new_xml, count=1, flags=re.DOTALL)

    # === notes (GenBank text 表記) ===
    notes_text = build_genbank_notes(
        accession_id, target_organism, consensus, primer_name,
        forward, reverse, fw_start_1, fw_end_1, rv_start_1, rv_end_1, accessions_used,
    )
    new_xml = re.sub(r"<notes>.*?</notes>",
                     f"<notes>{escape_xml(notes_text)}</notes>",
                     new_xml, count=1, flags=re.DOTALL)

    return new_xml


def build_genbank_notes(accession_id, organism, consensus, primer_name,
                        forward, reverse, fw_s, fw_e, rv_s, rv_e, accessions):
    """LOCUS ... // 形式の GenBank テキストを構築"""
    width = len(consensus)
    today = datetime.now().strftime("%d-%b-%Y").upper()
    accs_lines = "\n".join(f"            - {a}" for a in accessions)
    lines = [
        f"LOCUS       {accession_id:<22}{width} bp    DNA     linear   INV {today}",
        f"DEFINITION  {organism} primer alignment ({primer_name}) from {len(accessions)} accessions.",
        f"ACCESSION   {accession_id}",
        f"VERSION     {accession_id}.1",
        f"KEYWORDS    GCA primer; primer alignment.",
        f"SOURCE      {organism}",
        f"  ORGANISM  {organism}",
        f"            Unclassified.",
        f"COMMENT     Primer alignment generated from the following accessions:",
        accs_lines,
        f"",
        f"            Forward primer: {forward}",
        f"            Reverse primer: {reverse}",
        f"FEATURES             Location/Qualifiers",
        f"     source          1..{width}",
        f'                     /organism="{organism}"',
        f'                     /mol_type="genomic DNA"',
        f"     misc_feature    {fw_s}..{rv_e}",
        f'                     /label="PCR_product"',
        f"     primer_bind     {fw_s}..{fw_e}",
        f'                     /label="Forward"',
        f'                     /note="{forward}"',
        f"     primer_bind     complement({rv_s}..{rv_e})",
        f'                     /label="Reverse"',
        f'                     /note="{reverse}"',
        f"ORIGIN",
    ]
    # 60bp / line, 1-based 数字
    for i in range(0, width, 60):
        chunk = consensus[i:i + 60].lower()
        groups = " ".join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        lines.append(f"{i + 1:>9} {groups}")
    lines.append("//")
    return "\n".join(lines) + "\n"


# ===== history 構築 =====
def build_history(template_history: str, accession_id: str, consensus: str,
                  new_annotations_xml: str) -> str:
    """テンプレート history を新しい seq/annotation MD5 で更新"""
    seq_md5 = md5_hex(consensus)
    anno_md5 = md5_hex(new_annotations_xml)
    h = template_history
    h = re.sub(r"<seqh123>.*?</seqh123>",
               f"<seqh123>{consensus}</seqh123>", h, count=1, flags=re.DOTALL)
    h = re.sub(r"<seqMD5h123>.*?</seqMD5h123>",
               f"<seqMD5h123>{seq_md5}</seqMD5h123>", h, count=1)
    h = re.sub(r"<annoMD5h123>.*?</annoMD5h123>",
               f"<annoMD5h123>{anno_md5}</annoMD5h123>", h, count=1)
    h = re.sub(r"<fileNameh123>.*?</fileNameh123>",
               f"<fileNameh123>{escape_xml(accession_id)}</fileNameh123>", h, count=1)
    h = re.sub(r"<annotationsXMLh123>.*?</annotationsXMLh123>",
               f"<annotationsXMLh123>{new_annotations_xml}</annotationsXMLh123>",
               h, count=1, flags=re.DOTALL)
    return h


# セキュリティ: zip-bomb / zip-slip 防御の上限値
MAX_TEMPLATE_MEMBER_SIZE = 50 * 1024 * 1024   # 50 MB / member
MAX_TEMPLATE_TOTAL_SIZE = 200 * 1024 * 1024   # 200 MB / total


# ===== .cow パッキング =====
def pack_cow(output_path: Path, template_path: Path, new_xml: str,
             new_history: str, new_spw_bytes: bytes):
    """テンプレートを基に新しい .cow ファイルを作成

    セキュリティ:
        - zip member サイズ上限を超えたら拒否 (zip bomb 防御)
        - メンバー名に '..' / 絶対パスが含まれる場合は拒否 (zip slip 防御)
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(template_path, "r") as zin:
        # 事前検査
        total_size = 0
        for info in zin.infolist():
            name = info.filename
            if name.startswith("/") or name.startswith("\\") or ".." in name.replace("\\", "/").split("/"):
                raise ValueError(f"テンプレート .cow に不正なパスが含まれています: {name}")
            if info.file_size > MAX_TEMPLATE_MEMBER_SIZE:
                raise ValueError(
                    f"テンプレート .cow のメンバー '{name}' が上限 ({MAX_TEMPLATE_MEMBER_SIZE} bytes) を超えています"
                )
            total_size += info.file_size
            if total_size > MAX_TEMPLATE_TOTAL_SIZE:
                raise ValueError(
                    f"テンプレート .cow の合計サイズが上限 ({MAX_TEMPLATE_TOTAL_SIZE} bytes) を超えています"
                )

        with zipfile.ZipFile(output_path, "w", zipfile.ZIP_DEFLATED) as zout:
            for name in zin.namelist():
                if name in ("xml", "xml2", "xml3"):
                    zout.writestr(name, new_xml.encode("utf-8"))
                elif name == "history":
                    zout.writestr(name, new_history.encode("utf-8"))
                elif name == "SPW":
                    zout.writestr(name, new_spw_bytes)
                else:
                    # NIW, xls, blastSaves はそのまま
                    zout.writestr(name, zin.read(name))


# ===== FASTA 読み込み =====
def load_records(path: Path):
    if not path.exists():
        return []
    return list(SeqIO.parse(path, "fasta"))


def detect_organism(records, fallback="Unknown organism"):
    if not records:
        return fallback
    desc = records[0].description
    parts = desc.split(maxsplit=3)
    if len(parts) >= 3:
        return f"{parts[1]} {parts[2]}"
    return fallback


def read_primer_from_csv(csv_path: Path, primer_id=None):
    with open(csv_path, "r") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"プライマー候補が空です: {csv_path}")
    if primer_id is not None:
        for r in rows:
            if str(r.get("ID", "")).strip() == str(primer_id).strip():
                return r
        raise ValueError(f"primer_id={primer_id} が見つかりません")
    return rows[0]


# ===== メイン =====
def main():
    parser = argparse.ArgumentParser(
        description="DNA Dynamo .cow ファイル生成（プライマー × NCBI配列アラインメント）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--sequences-dir", "-s",
                        help="target_sequence.fasta などを含むディレクトリ（--input-fasta 未指定時に必須）")
    parser.add_argument("--input-fasta",
                        help="単一 FASTA を入力にするモード（--sequences-dir の代替）")
    parser.add_argument("--input-label", default="TARGET",
                        choices=["TARGET", "NONTARGET", "RELATED", "ALL", "EXCLUDED"],
                        help="--input-fasta 使用時の全レコードのラベル (default: TARGET)")
    parser.add_argument("--output", "-o", required=True,
                        help="出力 .cow ファイルパス")
    parser.add_argument("--primer-name", "-n", default="Primer",
                        help="プライマー名（表示用）")
    parser.add_argument("--primers-csv", "-p", help="プライマー候補CSV")
    parser.add_argument("--primer-id", help="CSV の ID 列で選択（未指定なら先頭）")
    parser.add_argument("--forward", help="Forward プライマー配列（CSV不使用時）")
    parser.add_argument("--reverse", help="Reverse プライマー配列（CSV不使用時）")
    parser.add_argument("--template", default=DEFAULT_TEMPLATE,
                        help=(f"テンプレート .cow ファイル (default: {Path(DEFAULT_TEMPLATE).name}). "
                              "セキュリティ上の理由から、必ず信頼できる .cow ファイル "
                              "(自分で生成したもの、または同梱のテンプレート) のみ指定してください。"))
    parser.add_argument("--flank", type=int, default=30,
                        help="プライマー外側に表示する flanking 塩基数 (default: 30)")
    parser.add_argument("--include-related", action="store_true",
                        help="related_sequences.fasta も含める（非標的扱い）")

    args = parser.parse_args()

    # プライマー
    if args.primers_csv:
        row = read_primer_from_csv(Path(args.primers_csv), args.primer_id)
        forward = row["Forward"].strip()
        reverse = row["Reverse"].strip()
        if args.primer_name == "Primer":
            primer_name = f"Primer-{row.get('ID', 'X')}"
        else:
            primer_name = args.primer_name
    else:
        if not (args.forward and args.reverse):
            print("エラー: --primers-csv または --forward/--reverse を指定してください")
            sys.exit(1)
        forward = args.forward.strip()
        reverse = args.reverse.strip()
        primer_name = args.primer_name

    output_path = Path(args.output)
    template_path = Path(args.template)

    if not template_path.exists():
        print(f"エラー: テンプレートが見つかりません: {template_path}")
        sys.exit(1)

    # 配列読み込み
    if args.input_fasta:
        # 単一 FASTA モード
        input_path = Path(args.input_fasta)
        if not input_path.exists():
            print(f"エラー: --input-fasta が見つかりません: {input_path}")
            sys.exit(1)
        all_records = load_records(input_path)
        label = args.input_label
        records_with_label = [(r, label) for r in all_records]
        target = all_records  # organism 検出用
        nontarget = []
        related = []
        print("=" * 60)
        print("DNA Dynamo .cow 生成 (single FASTA mode)")
        print("=" * 60)
        print(f"  Primer name : {primer_name}")
        print(f"  Forward     : {forward}")
        print(f"  Reverse     : {reverse}")
        print(f"  Input FASTA : {input_path.name}")
        print(f"  Label       : {label}")
        print(f"  Records     : {len(all_records)}")
        print()
        if not all_records:
            print("エラー: 配列が1つもありません")
            sys.exit(1)
    else:
        if not args.sequences_dir:
            print("エラー: --sequences-dir または --input-fasta のいずれかを指定してください")
            sys.exit(1)
        sequences_dir = Path(args.sequences_dir)
        target = load_records(sequences_dir / "target_sequence.fasta")
        nontarget = load_records(sequences_dir / "nontarget_sequences.fasta")
        related = load_records(sequences_dir / "related_sequences.fasta") if args.include_related else []

        print("=" * 60)
        print("DNA Dynamo .cow 生成")
        print("=" * 60)
        print(f"  Primer name : {primer_name}")
        print(f"  Forward     : {forward}")
        print(f"  Reverse     : {reverse}")
        print(f"  Template    : {template_path.name}")
        print(f"  Sequences   :")
        print(f"    target    : {len(target)}")
        print(f"    nontarget : {len(nontarget)}")
        if args.include_related:
            print(f"    related   : {len(related)}")
        print()

        if not target and not nontarget and not related:
            print("エラー: 配列が1つもありません")
            sys.exit(1)

        # ラベル付き配列リスト
        records_with_label = (
            [(r, "TARGET") for r in target]
            + [(r, "RELATED") for r in related]
            + [(r, "NONTARGET") for r in nontarget]
        )

    # アラインメント構築
    alignment, skipped = build_alignment(records_with_label, forward, reverse, flank=args.flank)
    if alignment is None:
        print("エラー: どの配列にも Forward プライマーが見つかりませんでした")
        sys.exit(1)

    print(f"  Alignment width : {alignment['width']} bp")
    print(f"  Forward column  : {alignment['fw_start'] + 1}..{alignment['fw_end']}")
    print(f"  Aligned rows    : {len(alignment['rows'])}")
    if skipped:
        print(f"  Skipped         : {len(skipped)} 配列 (forward 結合無し)")
        for sid, lbl in skipped:
            print(f"    - [{lbl}] {sid}")

    print()
    print("  ミスマッチ概要:")
    for r in alignment["rows"]:
        rv = "?" if r["rv_mismatches"] < 0 else str(r["rv_mismatches"])
        print(f"    [{r['label']:9}] {r['name']:18} F mm={r['fw_mismatches']}, R mm={rv}")

    # テンプレート読み込み
    print("\n  テンプレート展開中...")
    with zipfile.ZipFile(template_path, "r") as zin:
        template_xml = zin.read("xml").decode("utf-8")
        template_history = zin.read("history").decode("utf-8")
        template_spw_bytes = zin.read("SPW")

    # SPW
    template_spw = javaobj.loads(template_spw_bytes)
    new_spw = build_spw(template_spw, alignment, primer_name)
    new_spw_bytes = javaobj.dumps(new_spw)

    # xml
    target_organism = detect_organism(target, fallback="Target organism")
    accessions_used = [r["name"] for r in alignment["rows"]]
    new_xml = build_xml(template_xml, alignment, primer_name, forward, reverse,
                        accessions_used, target_organism)

    # history
    annotations_xml_block = re.search(r"<annotationBox>(.*?)</annotationBox>",
                                      new_xml, re.DOTALL).group(1)
    accession_id = f"{primer_name}_alignment"
    new_history = build_history(template_history, accession_id,
                                alignment["consensus"], annotations_xml_block)

    # パッキング
    pack_cow(output_path, template_path, new_xml, new_history, new_spw_bytes)

    # 同名 FASTA を出力（手動ドラッグ＆ドロップ用）
    fasta_path = output_path.with_suffix(".fasta")
    write_alignment_fasta(fasta_path, alignment, primer_name)

    # 同名 GenBank を出力（"Open Sequences In Alignment Editor" 用）
    gb_path = output_path.with_suffix(".gb")
    write_alignment_genbank(gb_path, alignment, primer_name, forward, reverse, target_organism)

    print(f"\n  ✓ {output_path.name}  ({output_path.stat().st_size} bytes)")
    print(f"  ✓ {fasta_path.name}  ({fasta_path.stat().st_size} bytes, {len(alignment['rows'])} 配列)")
    print(f"  ✓ {gb_path.name}     ({gb_path.stat().st_size} bytes, consensus + {len(alignment['rows'])} 配列)")
    print()
    print("推奨ワークフロー（マルチレコード GenBank 経由）:")
    print(f"  open -a DNADynamo {gb_path}")
    print(f"  → ダイアログで全配列を選択 → 'Open Sequences In Alignment Editor' クリック")
    print(f"  → コンセンサスリファレンス + Forward/Reverse プライマー feature 付きで")
    print(f"     24配列がアラインメント表示される")


if __name__ == "__main__":
    main()
