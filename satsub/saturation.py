"""Substitution-saturation analysis: pairwise s/v counts and genetic
distances (JC69, K80, TN93, GTR) for a nucleotide alignment, optionally
restricted to a single codon position.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Literal

import numpy as np
import pandas as pd
from Bio.Align import MultipleSeqAlignment

from satsub.codon import PositionLabel, columns_for_subset
from satsub.distances import (
    GTRModel,
    SiteCounts,
    base_frequencies,
    encode_sequence,
    jc69_distance,
    k80_distance,
    pairwise_site_counts,
    tn93_distance,
)

POSITION_LABELS: dict[PositionLabel, str] = {
    "all": "All positions",
    "1": "1st codon position",
    "2": "2nd codon position",
    "3": "3rd codon position",
}

MODEL_LABELS = {
    "jc69": "JC69 (Jukes & Cantor 1969)",
    "k80": "K80 (Kimura 1980)",
    "tn93": "TN93 (Tamura & Nei 1993)",
    "gtr": "GTR (general time-reversible)",
}


@dataclass
class SaturationResult:
    subset: PositionLabel
    frame_start: int
    n_sites_used: int
    freqs: np.ndarray
    gtr_model: GTRModel
    table: pd.DataFrame  # one row per sequence pair


def compute_saturation_table(
    alignment: MultipleSeqAlignment,
    subset: PositionLabel = "all",
    frame_start: int = 1,
) -> SaturationResult:
    """Compute pairwise s/v proportions and JC69/K80/TN93/GTR distances for
    every sequence pair, restricted to the requested codon-position subset.
    """
    names = [rec.id for rec in alignment]
    seqs = [str(rec.seq) for rec in alignment]
    encoded = [encode_sequence(s) for s in seqs]

    aln_len = alignment.get_alignment_length()
    positions = columns_for_subset(aln_len, subset, frame_start)
    n_sites_used = len(positions) if positions is not None else aln_len

    freqs = base_frequencies(encoded, positions)

    # Pool substitution counts across all pairs once, to fit a single GTR
    # rate matrix shared by every pair in this subset (see distances.py).
    pooled = np.zeros((4, 4), dtype=np.int64)
    pair_counts: list[tuple[int, int, SiteCounts]] = []
    for i, j in combinations(range(len(seqs)), 2):
        sc = pairwise_site_counts(encoded[i], encoded[j], positions)
        pair_counts.append((i, j, sc))
        pooled += sc.counts4x4

    gtr_model = GTRModel.fit(pooled, freqs)

    rows = []
    for i, j, sc in pair_counts:
        row = {
            "seq1": names[i],
            "seq2": names[j],
            "n_sites": sc.n_valid,
            "n_excluded": sc.n_excluded,
            "p_distance": sc.p_distance,
            "Ps": sc.Ps,
            "Pv": sc.Pv,
            "ts_tv_ratio": (sc.Ps / sc.Pv) if sc.Pv else float("nan"),
            "jc69": jc69_distance(sc),
            "k80": k80_distance(sc),
            "tn93": tn93_distance(sc, freqs),
            "gtr": gtr_model.pairwise_ml_distance(sc.counts4x4),
        }
        rows.append(row)

    table = pd.DataFrame(rows)
    return SaturationResult(
        subset=subset,
        frame_start=frame_start,
        n_sites_used=n_sites_used,
        freqs=freqs,
        gtr_model=gtr_model,
        table=table,
    )


def compute_all_subsets(
    alignment: MultipleSeqAlignment, frame_start: int = 1
) -> dict[PositionLabel, SaturationResult]:
    """Convenience wrapper: compute the saturation table for all positions
    plus each of the three codon positions."""
    results = {}
    for subset in ("all", "1", "2", "3"):
        try:
            results[subset] = compute_saturation_table(alignment, subset, frame_start)
        except ValueError:
            continue
    return results
