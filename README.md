# GCA Primer Auto Design

**A macOS toolkit that fully automates the workflow from species name input → NCBI sequence fetch → Primer3 design → DNA Dynamo visualization in a single command.**

Built for researchers designing primers for **GCA (Gut Content Analysis)** — amplifying the target species DNA while avoiding non-target and related species.

> 日本語版: [README.ja.md](README.ja.md)

---

## Features

- **One-command interactive workflow** — run `./auto_full.sh`, answer 5 prompts, done
- **Two design modes**
  - **GCA mode (default)**: 100–160 bp products for degraded DNA / gut-content analysis
  - **General mode**: 200–1000 bp products for general specific PCR, cloning, sequencing
- **AT-rich auto-fallback** — automatic `--at-rich` retry when Primer3 fails (for Lepidoptera / Diptera mtDNA)
- **4-category visualization** — target / non-target / all sequences / excluded species, each as a separate DNA Dynamo `.cow` file
- **Full macOS DNA Dynamo automation** — Quartz API-based automation of the Drag Drop Assembly Window (optional)
- **Per-species isolation** — each species lands in its own `reports/<Species>_<Gene>/` folder

---

## Quick Start

```bash
git clone https://github.com/hirschmannielladiversa-bot/gca-primer-design.git
cd gca-primer-design

# Dependencies
pip3 install -r requirements.txt

# Optional (for DNA Dynamo automation)
pip3 install pyobjc-framework-Quartz

# Make scripts executable
chmod +x auto_full.sh auto_visualize.sh verify_slw_count.sh

# Run
./auto_full.sh
```

Interactive prompts:

```
Enter target species (e.g. Helicoverpa armigera): Helicoverpa armigera
Enter target gene (e.g. COI, matK): COI
Related species (comma-separated) [optional]: Helicoverpa assulta,Spodoptera exigua
Non-target organisms (comma-separated) [optional]: Spodoptera litura
Use AT-rich mode? (Lepidoptera etc. mtDNA) [y/N]: y
```

Output:

```
reports/Helicoverpa_armigera_COI/プライマー結合/
├── ターゲット配列のみ.cow      (Target sequences only)
├── 非標的配列のみ.cow          (Non-target sequences only)
├── すべての配列.cow            (All sequences)
└── 除外した種_Spodoptera_litura.cow  (Excluded species)
```

> **Note**: Output filenames are intentionally in Japanese for direct display in DNA Dynamo window titles. Category labels above are English translations.

---

## Documentation

For detailed usage, AT-rich explanation, troubleshooting, and a Q&A section, open **[manual.html](manual.html)** in your browser.

日本語マニュアルは [manual.ja.html](manual.ja.html) を参照してください。

---

## Command-line Usage

### Fully interactive

```bash
./auto_full.sh
```

### Non-interactive (CLI arguments)

```bash
# GCA mode (default, 100–160 bp)
./auto_full.sh "Helicoverpa armigera" COI
./auto_full.sh "Helicoverpa armigera" COI --at-rich

# General mode (200–1000 bp, for standard specific PCR)
./auto_full.sh "Helicoverpa armigera" COI --mode general
./auto_full.sh "Helicoverpa armigera" COI --mode general --product-size-min 300 --product-size-max 700

# Related / non-target species
./auto_full.sh "Bemisia tabaci" COI --related "Trialeurodes vaporariorum"

# Already designed or file-only
./auto_full.sh "Helicoverpa armigera" COI --skip-design   # Skip design, visualize only
./auto_full.sh "Helicoverpa armigera" COI --no-automation # Generate files, skip DNA Dynamo
```

### Verify generated files

```bash
./verify_slw_count.sh Helicoverpa_armigera_COI
```

---

## Requirements

- **macOS** 10.15+ (DNA Dynamo automation is macOS-only; file generation runs on Linux/Windows)
- **Python** 3.10+
- **DNA Dynamo** (only needed for the visualization/automation step)

---

## How "recommended primer" is chosen

The top row of `primer_candidates.csv` sorted by `Quality_Score` descending. Primer3 + crosscheck score combined. Use `--top-n 3` to process multiple top candidates.

---

## AT-rich Mode

Many insect mitochondrial DNA (especially Lepidoptera, Diptera, Hymenoptera) is highly AT-rich (>70%), making Primer3 default parameters unable to find candidates. The `--at-rich` flag relaxes:

| Parameter | Default | --at-rich |
|---|---|---|
| `PRIMER_MIN_TM` | 57°C | **50°C** |
| `PRIMER_MAX_TM` | 63°C | **60°C** |
| `PRIMER_MIN_GC` | 40% | **20%** |

**You don't have to decide upfront**: if Primer3 returns no candidates with default settings, `auto_full.py` automatically retries with `--at-rich`.

---

## Security

This tool was audited by two parallel security review agents (code-injection perspective + path-traversal / data-exfiltration perspective) before initial release. Fixes applied:

- **CRITICAL ×1**: Python heredoc injection in `verify_slw_count.sh` — fixed
- **HIGH ×5**: absolute path leakage, PATH hijack vector, zip-slip / zip-bomb defense, path traversal — all fixed
- **MEDIUM ×6**: input sanitization in species / gene / project names — all fixed

Subprocess calls use list-form `shell=False`, AppleScript arguments are passed via `argv` (not string interpolation), and all file paths are validated via `Path.resolve() + relative_to()` to stay within `reports/`.

See [manual.html](manual.html) for details.

---

## License

MIT License — see [LICENSE](LICENSE).

External libraries:

- [Biopython](https://biopython.org/) — Biopython License
- [primer3-py](https://github.com/libnano/primer3-py) — GPL v2
- [javaobj-py3](https://github.com/tcalmant/python-javaobj) — Apache 2.0
- [pyobjc](https://pyobjc.readthedocs.io/) — MIT

---

## Author

**TKG.M**

Bug reports and feature requests welcome via GitHub Issues.
