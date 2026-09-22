"""Codon-position column selection for nucleotide alignments.

Codon position is assigned purely from alignment-column index and a chosen
reading-frame start; it does not attempt to re-derive frames around indels.
For coding alignments without frameshifting gaps (the common case for
codon-aware/protein-back-translated alignments) this matches the classic
1st/2nd/3rd position convention used by tools such as MEGA and DAMBE.
"""

from __future__ import annotations

from typing import Literal

PositionLabel = Literal["all", "1", "2", "3"]


def codon_position_columns(alignment_length: int, frame_start: int = 1) -> dict[int, list[int]]:
    """Return 0-based column indices grouped by codon position (1, 2, 3).

    `frame_start` is 1-based (1, 2, or 3): the alignment column at which the
    first complete codon begins. Any leading columns before the frame start,
    and any trailing columns that do not complete a full codon, are excluded
    from all position groups (consistent with how codon position tables are
    computed in standard phylogenetics software).
    """
    if frame_start not in (1, 2, 3):
        raise ValueError("frame_start must be 1, 2, or 3")
    start = frame_start - 1
    usable = ((alignment_length - start) // 3) * 3
    cols = {1: [], 2: [], 3: []}
    for i in range(usable):
        col = start + i
        pos = (i % 3) + 1
        cols[pos].append(col)
    return cols


def columns_for_subset(alignment_length: int, subset: PositionLabel, frame_start: int = 1) -> list[int] | None:
    """Return the column indices for a requested subset ('all', '1', '2', '3').

    Returns None for 'all' (meaning: use every column, no filtering needed).
    """
    if subset == "all":
        return None
    positions = codon_position_columns(alignment_length, frame_start)
    return positions[int(subset)]
