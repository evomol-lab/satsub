"""Nucleotide composition, codon usage, and directional base-pair frequency
statistics for nucleotide alignments, overall and by codon position.

These mirror the classic "Nucleotide Frequencies", "Codon Usage" and
"Nucleotide Pair Frequency" reports found in tools such as DAMBE: per-taxon
values averaged across the alignment, plus pooled pairwise substitution
counts averaged per sequence pair.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np
import pandas as pd
from Bio.Align import MultipleSeqAlignment
from Bio.Data import CodonTable

from satsub.codon import codon_position_columns, columns_for_subset
from satsub.distances import BASES, encode_sequence, pairwise_site_counts
from satsub.io import GENETIC_CODE_CHOICES

DNA_BASES = "ACGT"
PAIR_ORDER = ["T", "C", "A", "G"]  # display convention: pyrimidines then purines


# --------------------------------------------------------------------------
# Nucleotide frequencies (per sequence, overall + per codon position)
# --------------------------------------------------------------------------


def _pct_block(chars: str) -> tuple[dict, int]:
    counts = {b: chars.count(b) for b in DNA_BASES}
    total = sum(counts.values())
    pct = {b: (100 * counts[b] / total if total else float("nan")) for b in DNA_BASES}
    return pct, total


def nucleotide_composition_table(alignment: MultipleSeqAlignment, frame_start: int = 1) -> pd.DataFrame:
    """Per-sequence T/C/A/G percentages, overall and at each codon position
    (when the alignment is long enough for the chosen frame), with a
    trailing 'Avg.' row: the simple mean of each column across sequences."""
    names = [rec.id for rec in alignment]
    seqs = [str(rec.seq) for rec in alignment]
    aln_len = alignment.get_alignment_length()

    cols = None
    try:
        candidate = codon_position_columns(aln_len, frame_start)
        if candidate[1]:
            cols = candidate
    except ValueError:
        cols = None

    rows = []
    for name, seq in zip(names, seqs):
        row: dict = {"sequence": name}
        pct, total = _pct_block(seq)
        for b in DNA_BASES:
            row[f"{b}_percent"] = pct[b]
        row["Total_sites"] = total
        if cols is not None:
            for pos in (1, 2, 3):
                chars = "".join(seq[c] for c in cols[pos])
                pct_p, total_p = _pct_block(chars)
                for b in DNA_BASES:
                    row[f"{b}{pos}_percent"] = pct_p[b]
                row[f"Pos{pos}_sites"] = total_p
        rows.append(row)

    df = pd.DataFrame(rows)
    numeric_cols = [c for c in df.columns if c != "sequence"]
    avg_row = {"sequence": "Avg."}
    for c in numeric_cols:
        avg_row[c] = df[c].mean()
    df = pd.concat([df, pd.DataFrame([avg_row])], ignore_index=True)
    return df.round(2)


# --------------------------------------------------------------------------
# Codon usage
# --------------------------------------------------------------------------


@dataclass
class CodonUsageResult:
    table: pd.DataFrame  # codon (RNA notation), amino_acid, count_avg, count_total, rscu
    avg_codons_per_sequence: float
    table_id: int
    table_name: str


def codon_usage_table(
    alignment: MultipleSeqAlignment, frame_start: int = 1, table_id: int = 1
) -> CodonUsageResult:
    """Average codon usage (mean raw count per sequence, matching the
    classic per-taxon-averaged codon usage report) and relative synonymous
    codon usage (RSCU), computed from that averaged table and grouped by
    amino acid under the chosen genetic code."""
    aln_len = alignment.get_alignment_length()
    seqs = [str(rec.seq) for rec in alignment]
    start = frame_start - 1
    usable = ((aln_len - start) // 3) * 3
    if usable <= 0:
        raise ValueError("Alignment too short for the selected reading frame.")

    codon_table = CodonTable.unambiguous_dna_by_id[table_id]
    forward = dict(codon_table.forward_table)
    for stop in codon_table.stop_codons:
        forward[stop] = "*"

    all_codons = [a + b + c for a in DNA_BASES for b in DNA_BASES for c in DNA_BASES]

    per_seq_counts = []
    per_seq_totals = []
    for seq in seqs:
        counts = dict.fromkeys(all_codons, 0)
        total = 0
        for k in range(start, start + usable, 3):
            codon = seq[k : k + 3]
            if all(ch in DNA_BASES for ch in codon):
                counts[codon] += 1
                total += 1
        per_seq_counts.append(counts)
        per_seq_totals.append(total)

    avg_count = {c: float(np.mean([pc[c] for pc in per_seq_counts])) for c in all_codons}
    total_count = {c: int(sum(pc[c] for pc in per_seq_counts)) for c in all_codons}

    groups: dict[str, list[str]] = {}
    for c in all_codons:
        groups.setdefault(forward[c], []).append(c)

    rscu = {}
    for codons in groups.values():
        fam_mean = float(np.mean([avg_count[c] for c in codons]))
        for c in codons:
            rscu[c] = (avg_count[c] / fam_mean) if fam_mean > 0 else 0.0

    rows = [
        {
            "codon": c.replace("T", "U"),
            "amino_acid": forward[c],
            "count_avg": avg_count[c],
            "count_total": total_count[c],
            "rscu": rscu[c],
        }
        for c in sorted(all_codons, key=lambda c: (forward[c], c))
    ]
    table = pd.DataFrame(rows)
    table["count_avg"] = table["count_avg"].round(2)
    table["rscu"] = table["rscu"].round(3)

    return CodonUsageResult(
        table=table,
        avg_codons_per_sequence=float(np.mean(per_seq_totals)),
        table_id=table_id,
        table_name=GENETIC_CODE_CHOICES.get(table_id, str(table_id)),
    )


# --------------------------------------------------------------------------
# Directional base-pair frequencies
# --------------------------------------------------------------------------


@dataclass
class DirectionalPairResult:
    summary: pd.DataFrame  # subset, n_sites, ii, si, sv, R
    matrices: dict  # subset -> 4x4 DataFrame (rows/cols in T,C,A,G order), per-pair-averaged counts


def directional_pair_frequencies(
    alignment: MultipleSeqAlignment, frame_start: int = 1
) -> DirectionalPairResult:
    """Pooled directional substitution counts (row base = base in the
    earlier-listed sequence of each pair, column base = base in the
    later-listed sequence), averaged over every sequence pair and rounded,
    for the full alignment and for each codon position. Also reports mean
    identical (ii), transitional (si) and transversional (sv) pair counts
    per comparison and their ratio R = si/sv."""
    seqs = [str(rec.seq) for rec in alignment]
    encoded = [encode_sequence(s) for s in seqs]
    n = len(seqs)
    aln_len = alignment.get_alignment_length()
    base_idx = {b: BASES.index(b) for b in PAIR_ORDER}

    summary_rows = []
    matrices: dict = {}
    for subset in ("all", "1", "2", "3"):
        try:
            positions = columns_for_subset(aln_len, subset, frame_start)
        except ValueError:
            continue
        n_sites = len(positions) if positions is not None else aln_len

        pooled = np.zeros((4, 4), dtype=np.int64)
        ts = tv = ident = 0
        npairs = 0
        for i, j in combinations(range(n), 2):
            sc = pairwise_site_counts(encoded[i], encoded[j], positions)
            pooled += sc.counts4x4
            ts += sc.n_ts
            tv += sc.n_tv
            ident += sc.n_valid - sc.n_ts - sc.n_tv
            npairs += 1
        if npairs == 0:
            continue

        avg_matrix = pooled / npairs
        mean_si = ts / npairs
        mean_sv = tv / npairs
        mean_ii = ident / npairs
        R = (mean_si / mean_sv) if mean_sv else float("nan")

        summary_rows.append(
            {
                "subset": subset,
                "n_sites": n_sites,
                "ii": round(mean_ii),
                "si": round(mean_si),
                "sv": round(mean_sv),
                "R": round(R, 1) if mean_sv else float("nan"),
            }
        )
        reordered = np.array([[avg_matrix[base_idx[r], base_idx[c]] for c in PAIR_ORDER] for r in PAIR_ORDER])
        matrices[subset] = pd.DataFrame(np.round(reordered, 1), index=PAIR_ORDER, columns=PAIR_ORDER)

    return DirectionalPairResult(summary=pd.DataFrame(summary_rows), matrices=matrices)
