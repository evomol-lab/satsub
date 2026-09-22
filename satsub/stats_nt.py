"""Descriptive statistics for nucleotide multiple sequence alignments."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np
import pandas as pd
from Bio.Align import MultipleSeqAlignment

from satsub.codon import codon_position_columns
from satsub.distances import encode_sequence, pairwise_site_counts

GAP_CHARS = set("-.?")
CORE_BASES = "ACGT"


@dataclass
class NucleotideSummary:
    n_sequences: int
    alignment_length: int
    overall: dict
    per_sequence: pd.DataFrame
    identity_matrix: pd.DataFrame

    def overall_table(self) -> pd.DataFrame:
        return pd.DataFrame(list(self.overall.items()), columns=["statistic", "value"])


def _ungapped_length(seq: str) -> int:
    return sum(1 for c in seq if c not in GAP_CHARS)


def summarize_nucleotide_alignment(
    alignment: MultipleSeqAlignment, frame_start: int | None = 1
) -> NucleotideSummary:
    names = [rec.id for rec in alignment]
    seqs = [str(rec.seq) for rec in alignment]
    n = len(seqs)
    aln_len = alignment.get_alignment_length()
    encoded = [encode_sequence(s) for s in seqs]

    # --- per-sequence composition, gaps, GC ---
    per_seq_rows = []
    for name, seq in zip(names, seqs):
        counts = {b: seq.count(b) for b in CORE_BASES}
        gaps = sum(seq.count(c) for c in GAP_CHARS)
        other = len(seq) - sum(counts.values()) - gaps
        ungapped = _ungapped_length(seq)
        gc = counts["G"] + counts["C"]
        at = counts["A"] + counts["T"]
        per_seq_rows.append(
            {
                "sequence": name,
                "ungapped_length": ungapped,
                "gap_percent": 100 * gaps / aln_len if aln_len else float("nan"),
                "A": counts["A"],
                "C": counts["C"],
                "G": counts["G"],
                "T": counts["T"],
                "ambiguous_or_other": other,
                "GC_percent": 100 * gc / (gc + at) if (gc + at) else float("nan"),
                "GC_skew": (counts["G"] - counts["C"]) / gc if gc else float("nan"),
                "AT_skew": (counts["A"] - counts["T"]) / at if at else float("nan"),
            }
        )
    per_sequence = pd.DataFrame(per_seq_rows)

    # --- overall composition / gaps / GC ---
    total_counts = {b: per_sequence[b].sum() for b in CORE_BASES}
    total_gaps = sum(seq.count(c) for seq in seqs for c in GAP_CHARS)
    total_chars = n * aln_len
    total_other = total_chars - sum(total_counts.values()) - total_gaps
    gc_total = total_counts["G"] + total_counts["C"]
    at_total = total_counts["A"] + total_counts["T"]

    # --- site classification (variable / conserved / parsimony-informative) ---
    n_conserved = n_variable = n_parsimony_informative = n_fully_gapped = 0
    for col in range(aln_len):
        column = [seqs[r][col] for r in range(n)]
        states = [c for c in column if c in CORE_BASES]
        if not states:
            n_fully_gapped += 1
            continue
        distinct = set(states)
        if len(distinct) == 1:
            n_conserved += 1
            continue
        n_variable += 1
        counts_in_col = {s: states.count(s) for s in distinct}
        if sum(1 for c in counts_in_col.values() if c >= 2) >= 2:
            n_parsimony_informative += 1
    n_singleton = n_variable - n_parsimony_informative

    # --- pairwise identity + transition/transversion summary ---
    identity = np.full((n, n), np.nan)
    np.fill_diagonal(identity, 100.0)
    ts_total = tv_total = valid_total = 0
    pairwise_p = []
    for i, j in combinations(range(n), 2):
        sc = pairwise_site_counts(encoded[i], encoded[j])
        pct = 100 * (1 - sc.p_distance) if sc.n_valid else float("nan")
        identity[i, j] = identity[j, i] = pct
        ts_total += sc.n_ts
        tv_total += sc.n_tv
        valid_total += sc.n_valid
        if sc.n_valid:
            pairwise_p.append(sc.p_distance)
    identity_matrix = pd.DataFrame(identity, index=names, columns=names)

    off_diag = identity[~np.eye(n, dtype=bool)]
    off_diag = off_diag[~np.isnan(off_diag)]

    overall = {
        "n_sequences": n,
        "alignment_length": aln_len,
        "mean_ungapped_length": float(per_sequence["ungapped_length"].mean()),
        "min_ungapped_length": int(per_sequence["ungapped_length"].min()),
        "max_ungapped_length": int(per_sequence["ungapped_length"].max()),
        "overall_gap_percent": 100 * total_gaps / total_chars if total_chars else float("nan"),
        "ambiguous_or_other_percent": 100 * total_other / total_chars if total_chars else float("nan"),
        "overall_GC_percent": 100 * gc_total / (gc_total + at_total) if (gc_total + at_total) else float("nan"),
        "conserved_sites": n_conserved,
        "variable_sites": n_variable,
        "parsimony_informative_sites": n_parsimony_informative,
        "singleton_sites": n_singleton,
        "fully_gapped_sites": n_fully_gapped,
        "mean_pairwise_identity_percent": float(np.mean(off_diag)) if off_diag.size else float("nan"),
        "min_pairwise_identity_percent": float(np.min(off_diag)) if off_diag.size else float("nan"),
        "max_pairwise_identity_percent": float(np.max(off_diag)) if off_diag.size else float("nan"),
        "mean_pairwise_p_distance": float(np.mean(pairwise_p)) if pairwise_p else float("nan"),
        "aggregate_transitions": ts_total,
        "aggregate_transversions": tv_total,
        "aggregate_ts_tv_ratio": (ts_total / tv_total) if tv_total else float("nan"),
    }

    if frame_start is not None:
        try:
            gc_codon = codon_position_gc(alignment, frame_start)
            overall.update(gc_codon)
        except ValueError:
            pass

    return NucleotideSummary(
        n_sequences=n,
        alignment_length=aln_len,
        overall=overall,
        per_sequence=per_sequence,
        identity_matrix=identity_matrix,
    )


def codon_position_gc(alignment: MultipleSeqAlignment, frame_start: int = 1) -> dict:
    """GC content at each codon position (GC1, GC2, GC3) plus the GC1/GC2
    average, pooled across all sequences. Requires the alignment to contain
    at least one complete codon in the chosen frame."""
    aln_len = alignment.get_alignment_length()
    cols = codon_position_columns(aln_len, frame_start)
    seqs = [str(rec.seq) for rec in alignment]
    result = {}
    gc12_parts = []
    for pos in (1, 2, 3):
        idx = cols[pos]
        if not idx:
            raise ValueError("Alignment too short for codon-position GC content.")
        gc = at = 0
        for seq in seqs:
            for c in idx:
                ch = seq[c]
                if ch in "GC":
                    gc += 1
                elif ch in "AT":
                    at += 1
        pct = 100 * gc / (gc + at) if (gc + at) else float("nan")
        result[f"GC{pos}_percent"] = pct
        if pos in (1, 2):
            gc12_parts.append(pct)
    result["GC12_mean_percent"] = float(np.mean(gc12_parts)) if gc12_parts else float("nan")
    return result
