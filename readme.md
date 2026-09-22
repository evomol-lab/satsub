# Estimator of the substitution saturation rates (SatSub)

## Introduction

SatSub is a Python-based tool for estimating substitution saturation rates in
coding and non-coding nucleotide sequences, with a browser-based graphical
interface (built with Streamlit). It gives you a fast, visual estimate of
substitution saturation, helping you identify the most appropriate model of
nucleotide substitution for your data. SatSub also computes a range of
descriptive statistics for both nucleotide and amino-acid multiple sequence
alignments (MSAs) — GC content, base/amino-acid composition, gap percentage,
pairwise identity, conserved/variable/parsimony-informative site counts, and
more — to help you understand your dataset before downstream analysis.

## Features

**Alignment statistics (nucleotide and amino acid)**
- Sequence counts, alignment length, ungapped sequence lengths
- Gap and ambiguous-character percentages (overall and per sequence)
- Base composition, GC content and GC skew/AT skew (nucleotide)
- Codon-position GC content: GC1, GC2, GC3, GC12 mean (nucleotide, coding sequences)
- Amino-acid composition and a hydrophobic/polar/acidic/basic breakdown (amino acid)
- Conserved, variable, parsimony-informative and singleton site counts
- Pairwise identity matrix (both), plus a BLOSUM62-based similarity matrix (amino acid)
- Translation of an in-frame nucleotide alignment to amino acids, with a choice
  of NCBI genetic code tables (e.g. Standard, Vertebrate/Invertebrate Mitochondrial)

**Substitution saturation analysis (nucleotide)**
- Pairwise transition (s) and transversion (v) proportions
- Genetic distance estimation under four substitution models:
  - **JC69** — Jukes & Cantor (1969)
  - **K80** — Kimura two-parameter (1980)
  - **TN93** — Tamura & Nei (1993)
  - **GTR** — general time-reversible (rates estimated once per alignment/subset,
    then a per-pair maximum-likelihood distance under the fixed rate matrix)
- s and v vs. genetic distance saturation plots for each model
- The same analysis split by codon position (1st, 2nd, 3rd), with a
  configurable reading-frame start

**Interface**
- A Streamlit GUI: file upload, a bundled example dataset (8 vertebrate COI
  sequences), interactive Plotly charts, and CSV export for every table

## Installation

Requires Python 3.10+.

```bash
# with uv (recommended, no system python3-venv package needed)
uv venv .venv
uv pip install --python .venv/bin/python -e .

# or with a standard venv + pip
python3 -m venv .venv
.venv/bin/pip install -e .
```

## Running the GUI

```bash
.venv/bin/streamlit run satsub/app.py
# or, after installation:
.venv/bin/satsub-gui
```

This opens the app in your browser (default: http://localhost:8501). Use the
sidebar to upload a nucleotide alignment (FASTA, Clustal, Phylip, Nexus or
Stockholm), or click "Use bundled example" to try it immediately with the
included 8-taxon vertebrate COI dataset (`VertCOI.fas`). An amino-acid
alignment can either be translated from the nucleotide alignment (pick the
matching genetic code table, e.g. table 2 for vertebrate mitochondrial genes
such as COI) or uploaded separately.

## Running the tests

```bash
.venv/bin/pytest -q
```

## Project layout

```
satsub/
  io.py          # alignment reading, format/kind detection, translation
  codon.py       # codon-position column selection
  distances.py   # JC69 / K80 / TN93 / GTR pairwise distance estimators
  stats_nt.py    # nucleotide MSA statistics
  stats_aa.py    # amino-acid MSA statistics
  saturation.py  # pairwise s/v + distance tables, per codon position
  plotting.py    # Plotly figure builders
  app.py         # Streamlit GUI
  launch.py      # `satsub-gui` console-script entry point
tests/           # pytest suite
```

## Method notes

Closed-form distances (JC69, K80, TN93) use the standard textbook formulas
(Nei & Kumar, *Molecular Evolution and Phylogenetics*, 2000) and return an
undefined (NaN) result where the correction breaks down at high divergence —
itself a saturation signal, since it means the observed differences can no
longer be reliably corrected.

GTR has no closed-form pairwise distance. SatSub estimates its six
exchangeability rates once per alignment (or per codon-position subset) from
substitution counts pooled across every sequence pair, then computes each
pair's distance as the maximum-likelihood branch length under that fixed
rate matrix. This is the standard way to obtain pairwise GTR-like distances
without fitting a full tree, but it remains an aggregate/composite
approximation — treat the saturation plots as an exploratory diagnostic
rather than a publication-grade phylogenetic analysis.

Codon position is assigned from alignment-column index and the chosen
reading-frame start only; it is not re-derived around indels, so alignments
with frameshifting gaps will have codon positions shift for all downstream
columns. Use a codon-aware alignment for such data.
