"""Descriptive statistics for amino-acid multiple sequence alignments."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np
import pandas as pd
from Bio.Align import MultipleSeqAlignment, substitution_matrices

GAP_CHARS = set("-.?")
STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"

# A simple, exhaustive 4-way physicochemical grouping of the 20 standard
# amino acids (each letter appears in exactly one group).
HYDROPHOBIC = set("AVLIMFWPG")
POLAR_UNCHARGED = set("STCYNQ")
ACIDIC = set("DE")
BASIC = set("KRH")

_BLOSUM62 = substitution_matrices.load("BLOSUM62")


@dataclass
class AminoAcidSummary:
    n_sequences: int
    alignment_length: int
    overall: dict
    per_sequence: pd.DataFrame
    identity_matrix: pd.DataFrame
    similarity_matrix: pd.DataFrame

    def overall_table(self) -> pd.DataFrame:
        return pd.DataFrame(list(self.overall.items()), columns=["statistic", "value"])


def _ungapped_length(seq: str) -> int:
    return sum(1 for c in seq if c not in GAP_CHARS)


def summarize_amino_acid_alignment(alignment: MultipleSeqAlignment) -> AminoAcidSummary:
    names = [rec.id for rec in alignment]
    seqs = [str(rec.seq) for rec in alignment]
    n = len(seqs)
    aln_len = alignment.get_alignment_length()

    per_seq_rows = []
    for name, seq in zip(names, seqs):
        counts = {aa: seq.count(aa) for aa in STANDARD_AA}
        gaps = sum(seq.count(c) for c in GAP_CHARS)
        ungapped = _ungapped_length(seq)
        other = len(seq) - sum(counts.values()) - gaps
        hydro = sum(counts[a] for a in HYDROPHOBIC)
        polar = sum(counts[a] for a in POLAR_UNCHARGED)
        acidic = sum(counts[a] for a in ACIDIC)
        basic = sum(counts[a] for a in BASIC)
        classified = hydro + polar + acidic + basic
        per_seq_rows.append(
            {
                "sequence": name,
                "ungapped_length": ungapped,
                "gap_percent": 100 * gaps / aln_len if aln_len else float("nan"),
                "ambiguous_or_other_percent": 100 * other / aln_len if aln_len else float("nan"),
                "hydrophobic_percent": 100 * hydro / classified if classified else float("nan"),
                "polar_uncharged_percent": 100 * polar / classified if classified else float("nan"),
                "acidic_percent": 100 * acidic / classified if classified else float("nan"),
                "basic_percent": 100 * basic / classified if classified else float("nan"),
            }
        )
    per_sequence = pd.DataFrame(per_seq_rows)

    total_counts = {aa: sum(seq.count(aa) for seq in seqs) for aa in STANDARD_AA}
    total_gaps = sum(seq.count(c) for seq in seqs for c in GAP_CHARS)
    total_chars = n * aln_len
    total_classified = sum(total_counts.values())

    # --- site classification ---
    n_conserved = n_variable = n_parsimony_informative = n_fully_gapped = 0
    for col in range(aln_len):
        column = [seqs[r][col] for r in range(n)]
        states = [c for c in column if c in STANDARD_AA]
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

    # --- pairwise identity and BLOSUM62-based similarity ---
    identity = np.full((n, n), np.nan)
    similarity = np.full((n, n), np.nan)
    np.fill_diagonal(identity, 100.0)
    np.fill_diagonal(similarity, 100.0)
    for i, j in combinations(range(n), 2):
        s1, s2 = seqs[i], seqs[j]
        n_valid = n_identical = n_similar = 0
        for a, b in zip(s1, s2):
            if a in STANDARD_AA and b in STANDARD_AA:
                n_valid += 1
                if a == b:
                    n_identical += 1
                    n_similar += 1
                else:
                    try:
                        score = _BLOSUM62[a, b]
                    except KeyError:
                        score = 0
                    if score > 0:
                        n_similar += 1
        if n_valid:
            identity[i, j] = identity[j, i] = 100 * n_identical / n_valid
            similarity[i, j] = similarity[j, i] = 100 * n_similar / n_valid
    identity_matrix = pd.DataFrame(identity, index=names, columns=names)
    similarity_matrix = pd.DataFrame(similarity, index=names, columns=names)

    off_diag = identity[~np.eye(n, dtype=bool)]
    off_diag = off_diag[~np.isnan(off_diag)]
    sim_off_diag = similarity[~np.eye(n, dtype=bool)]
    sim_off_diag = sim_off_diag[~np.isnan(sim_off_diag)]

    overall = {
        "n_sequences": n,
        "alignment_length": aln_len,
        "mean_ungapped_length": float(per_sequence["ungapped_length"].mean()),
        "min_ungapped_length": int(per_sequence["ungapped_length"].min()),
        "max_ungapped_length": int(per_sequence["ungapped_length"].max()),
        "overall_gap_percent": 100 * total_gaps / total_chars if total_chars else float("nan"),
        "ambiguous_or_other_percent": (
            100 * (total_chars - total_gaps - total_classified) / total_chars if total_chars else float("nan")
        ),
        "conserved_sites": n_conserved,
        "variable_sites": n_variable,
        "parsimony_informative_sites": n_parsimony_informative,
        "singleton_sites": n_singleton,
        "fully_gapped_sites": n_fully_gapped,
        "mean_pairwise_identity_percent": float(np.mean(off_diag)) if off_diag.size else float("nan"),
        "min_pairwise_identity_percent": float(np.min(off_diag)) if off_diag.size else float("nan"),
        "max_pairwise_identity_percent": float(np.max(off_diag)) if off_diag.size else float("nan"),
        "mean_pairwise_similarity_percent": float(np.mean(sim_off_diag)) if sim_off_diag.size else float("nan"),
    }

    return AminoAcidSummary(
        n_sequences=n,
        alignment_length=aln_len,
        overall=overall,
        per_sequence=per_sequence,
        identity_matrix=identity_matrix,
        similarity_matrix=similarity_matrix,
    )
