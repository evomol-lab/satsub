<p align="center">
  <img src="satsub/assets/satsub-logo.png" alt="SatSub logo" width="500">
</p>

# Estimator of the substitution saturation rates (SatSub)

## Introduction

[SatSub](https://satsub.streamlit.app) is a Python-based tool for estimating substitution saturation rates in
coding and non-coding nucleotide sequences, with a browser-based graphical
interface (built with Streamlit). It gives you a fast, visual estimate of
substitution saturation, helping you identify the most appropriate model of
nucleotide substitution for your data. SatSub also computes a range of
descriptive statistics for both nucleotide and amino-acid multiple sequence
alignments (MSAs) — GC content, base/amino-acid composition, gap percentage,
pairwise identity, conserved/variable/parsimony-informative site counts, and
more — to help you understand your dataset before downstream analysis.

>*This tool was developed for teaching puporses only. For research datasets use it at your own risk!*

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
- s and v vs. genetic distance saturation plots for each model, each with a
  LOWESS-smoothed trend line per series so a plateauing (saturating) trend is
  visible rather than forced onto a straight line
- The same analysis split by codon position (1st, 2nd, 3rd), with a
  configurable reading-frame start

**Composition & codon usage (nucleotide)**
- Nucleotide (T/C/A/G) frequencies per sequence, overall and at each codon
  position, with a grouped bar chart of the cross-sequence mean
- Codon usage: mean codon count per sequence and relative synonymous codon
  usage (RSCU) for all 64 codons under a chosen NCBI genetic code table,
  shown as a usage bar chart and an RSCU chart colored on a diverging scale
  centered at 1.0 (equal usage within a synonymous family)
- Directional base-pair frequencies: identical/transition/transversion pair
  counts averaged per sequence comparison (ii, si, sv, R = si/sv), plus the
  full 4×4 directional substitution matrix as a heatmap, for all sites and
  for each codon position

**Interface**
- A Streamlit GUI: file upload, a bundled example dataset (8 vertebrate COI
  sequences), interactive Plotly charts, and CSV export for every table.

SatSub can be executed from the following link: https://satsub.streamlit.app/

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
alignment can either be translated from the nucleotide alignment or
uploaded separately. The sidebar's genetic code table (e.g. table 2 for
vertebrate mitochondrial genes such as COI) is a single shared choice: it
drives both that translation and the codon usage/RSCU calculation in the
Composition & codon usage tab, so the two are always consistent.

New to SatSub? [`TUTORIAL.md`](TUTORIAL.md) is a guided, hands-on walkthrough
of every tab using this same bundled example dataset, ending with a worked
example of spotting substitution saturation at the 3rd codon position.

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
  composition.py # nucleotide frequencies, codon usage/RSCU, directional pair frequencies
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

Codon usage frequencies are the mean raw codon count per sequence, averaged
across taxa (not a per-mille rate), matching the classic per-taxon-averaged
codon usage report. RSCU (Sharp, Tuohy & Mosurski, 1986) is then derived
from that averaged table: each codon's mean count divided by the average
count of its synonymous family under the selected genetic code, so 1.0 is
the value expected if every synonym were used equally.

Directional base-pair frequencies pool substitution counts across every
sequence pair and divide by the number of pairs. Pooling over unordered
pairs has no natural direction, so SatSub fixes one convention throughout:
the row is the base in whichever sequence is listed earlier in the
alignment, the column is the base in the later one. This is an arbitrary
but consistent labeling, not a claim about ancestral state or direction of
change.

Saturation-plot trend lines use LOWESS (Cleveland, 1979), a local weighted
regression, rather than a single straight-line fit, so a plateauing trend
shows up as a visible bend instead of being averaged into a straight slope.

## License

SatSub is released under the [MIT License](LICENSE), by the
[EvoMol Lab](https://evomol-lab.imd.ufrn.br). It is free to use, modify,
and redistribute, including for commercial purposes, provided the copyright
notice and license text are preserved.

MIT matches the convention for bioinformatics tooling (Biopython,
scikit-bio, DendroPy and most of the scientific Python stack use MIT or
BSD) and is OSI-approved, which matters for inclusion in package
ecosystems such as Bioconda/conda-forge and Debian Med that filter or
hesitate on non-OSI-approved licenses. It is also compatible with every
runtime dependency SatSub uses: Biopython (Biopython License Agreement, a
BSD-3-Clause-equivalent permissive license), NumPy, SciPy and pandas
(BSD-3-Clause), Plotly (MIT), and Streamlit (Apache-2.0) — none of these
impose copyleft/share-alike terms, so there is no license conflict.

## Disclosure on Generative AI Use

The developer team used generative AI tools for the following tasks. Throughout this process, the authors maintained full control over the research design and interpretation of results; the AI acted solely as a technical and linguistic aid.

- Code writing, revision, and optimization.

- Elaborate documentation topic structure.

- Review english language.