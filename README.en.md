# nanodelfy

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24498244.svg)](https://doi.org/10.5281/zenodo.24498244)

**nanodelfy** (from *NanoDel*, its original development) classifies long reads (Nanopore/PacBio) from
CRISPR amplicons with large deletions between two cut sites, with support for
**Cas9 and Cas12a**, amplicons of any size, and deletions from ~120 bp up to
several kb.

nanodelfy is the unification of the two NanoDel app versions originally developed
for the AESC7 and MmC7 amplicons (FRT-FRT deletions), now generalized so any
lab can analyze its own amplicons.

- Two references: you only need the **WT** amplicon FASTA; nanodelfy builds the
  EDIT reference (amplicon with the expected deletion) automatically.
- Per-read classification: WT, exact deletion, extended deletion, deletion with
  insertion, inversion, truncated (non-informative) reads...
- Microindels around guide cut sites, with Cas9-specific rules (NGG PAM, blunt
  cut) and Cas12a rules (TTTV PAM, staggered cut).
- Two interfaces: **local web GUI** (drag and drop your FASTQ) and a
  reproducible **CLI** (`nanodelfy run`) with JSON presets.

## Installation

Requirements: Python 3.10+.

```bash
pip install nanodelfy          # from PyPI (once published)
```

From source (this repository):

```bash
git clone https://github.com/Mencia11268/nanodelfy.git
cd nanodelfy
pip install .
```

Dependencies: `mappy` (minimap2 for Python), `pysam`, `pandas`, `openpyxl`,
`fastapi`, `uvicorn`.

## Quick start

### Web GUI

```bash
nanodelfy web            # opens http://127.0.0.1:8765
```

Upload your FASTQ and the WT amplicon FASTA, enter the expected deletion
coordinates and guides, and hit Run.

### Command line

```bash
nanodelfy run \
  --fastq sample.fastq.gz \
  --wt amplicon_WT.fa \
  --out output \
  --del-start 1438 --del-end 1557 \
  --min-sv-deletion 100 \
  --guides "guide1: CCGAAGTTCCTATTCTCTAG" \
  --nuclease cas9
```

See the [Spanish README](README.md) for the full parameter table, output
specification, per-read classes, presets, and limitations (the Spanish version
is the primary documentation; this file is a summary).

## Citation

If you use nanodelfy, please cite the version archived on Zenodo:

> Mencia, A., Murillas, R., & Garcia Diez, M. (2026). nanodelfy 1.1.2.
> Zenodo. https://doi.org/10.5281/zenodo.24498244

To refer to all versions, use the concept DOI:
https://doi.org/10.5281/zenodo.24498243. Citation metadata is available in
`CITATION.cff`.
