"""Pairwise nucleotide substitution counting and genetic-distance estimators.

Implements four classic distance corrections used to assess substitution
saturation:

* JC69  - Jukes & Cantor (1969): equal base frequencies, equal rates.
* K80   - Kimura (1980) two-parameter: separate transition/transversion rates.
* TN93  - Tamura & Nei (1993): unequal base frequencies, separate purine and
          pyrimidine transition rates, single transversion rate.
* GTR   - General time-reversible: unequal frequencies, all six substitution
          types rated independently. GTR has no closed-form pairwise formula,
          so its parameters (frequencies + six exchangeabilities) are
          estimated once from data pooled across all sequence pairs in the
          alignment or codon-position subset being analysed, and then each
          pair's distance is obtained as the maximum-likelihood branch length
          under that fixed rate matrix. This "estimate once, profile per
          pair" strategy is the standard way pairwise GTR distances are
          computed without fitting a full tree (analogous to the composite-
          likelihood approach behind MEGA's Maximum Composite Likelihood
          method for pairwise distances).

All functions work on already-normalized nucleotide strings (uppercase,
DNA alphabet with U already converted to T by satsub.io).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
from scipy.optimize import minimize_scalar

BASES = "ACGT"
_BASE_TO_CODE = {b: i for i, b in enumerate(BASES)}
_PURINES = {0, 2}  # A, G
_PYRIMIDINES = {1, 3}  # C, T
RATE_PAIRS = [("A", "C"), ("A", "G"), ("A", "T"), ("C", "G"), ("C", "T"), ("G", "T")]


def encode_sequence(seq: str) -> np.ndarray:
    """Map a nucleotide string to int8 codes: A=0 C=1 G=2 T=3, anything else -1."""
    return np.array([_BASE_TO_CODE.get(c, -1) for c in seq], dtype=np.int8)


@dataclass
class SiteCounts:
    """Tally of aligned-site categories between one sequence pair."""

    n_valid: int
    n_excluded: int
    n_ts_purine: int  # A<->G
    n_ts_pyrimidine: int  # C<->T
    n_tv: int
    counts4x4: np.ndarray  # counts[i, j]: seq1 base i paired with seq2 base j (order A,C,G,T)

    @property
    def n_ts(self) -> int:
        return self.n_ts_purine + self.n_ts_pyrimidine

    @property
    def n_diff(self) -> int:
        return self.n_ts + self.n_tv

    @property
    def p_distance(self) -> float:
        return self.n_diff / self.n_valid if self.n_valid else float("nan")

    @property
    def Ps(self) -> float:
        """Observed transition proportion."""
        return self.n_ts / self.n_valid if self.n_valid else float("nan")

    @property
    def Pv(self) -> float:
        """Observed transversion proportion."""
        return self.n_tv / self.n_valid if self.n_valid else float("nan")

    @property
    def P1_purine(self) -> float:
        return self.n_ts_purine / self.n_valid if self.n_valid else float("nan")

    @property
    def P2_pyrimidine(self) -> float:
        return self.n_ts_pyrimidine / self.n_valid if self.n_valid else float("nan")


def pairwise_site_counts(
    encoded1: np.ndarray, encoded2: np.ndarray, positions: Sequence[int] | None = None
) -> SiteCounts:
    """Classify aligned sites between two pre-encoded sequences.

    `positions` restricts the comparison to the given 0-based column indices
    (e.g. a codon-position subset); pass None to use every column.
    Sites where either sequence has a gap or an ambiguous/non-ACGT
    character are excluded from all counts (`n_excluded`).
    """
    if positions is not None:
        a = encoded1[np.asarray(positions, dtype=np.intp)]
        b = encoded2[np.asarray(positions, dtype=np.intp)]
    else:
        a = encoded1
        b = encoded2

    valid = (a >= 0) & (b >= 0)
    n_excluded = int((~valid).sum())
    av, bv = a[valid], b[valid]
    n_valid = int(av.size)

    diff = av != bv
    a_purine = (av == 0) | (av == 2)
    b_purine = (bv == 0) | (bv == 2)
    ts_purine = diff & a_purine & b_purine
    a_pyr = (av == 1) | (av == 3)
    b_pyr = (bv == 1) | (bv == 3)
    ts_pyrimidine = diff & a_pyr & b_pyr
    tv = diff & ~ts_purine & ~ts_pyrimidine

    counts4x4 = np.zeros((4, 4), dtype=np.int64)
    if n_valid:
        flat_idx = av.astype(np.int64) * 4 + bv.astype(np.int64)
        bincount = np.bincount(flat_idx, minlength=16)
        counts4x4 = bincount.reshape(4, 4)

    return SiteCounts(
        n_valid=n_valid,
        n_excluded=n_excluded,
        n_ts_purine=int(ts_purine.sum()),
        n_ts_pyrimidine=int(ts_pyrimidine.sum()),
        n_tv=int(tv.sum()),
        counts4x4=counts4x4,
    )


def base_frequencies(
    encoded_seqs: Sequence[np.ndarray], positions: Sequence[int] | None = None
) -> np.ndarray:
    """Empirical A/C/G/T frequencies pooled across the given sequences.

    Falls back to equal (0.25 each) frequencies when there is no valid data,
    so downstream formulas never divide by zero.
    """
    counts = np.zeros(4, dtype=np.int64)
    for enc in encoded_seqs:
        arr = enc[np.asarray(positions, dtype=np.intp)] if positions is not None else enc
        valid = arr[arr >= 0]
        if valid.size:
            counts += np.bincount(valid, minlength=4)
    total = counts.sum()
    if total == 0:
        return np.full(4, 0.25)
    return counts / total


# --------------------------------------------------------------------------
# Closed-form distance corrections
# --------------------------------------------------------------------------


def jc69_distance(counts: SiteCounts) -> float:
    """Jukes & Cantor (1969) distance. NaN if the correction is undefined
    (p >= 0.75, i.e. sites are fully randomized / saturated)."""
    p = counts.p_distance
    if not counts.n_valid or math.isnan(p):
        return float("nan")
    x = 1 - (4 / 3) * p
    if x <= 0:
        return float("nan")
    return -0.75 * math.log(x)


def k80_distance(counts: SiteCounts) -> float:
    """Kimura (1980) two-parameter distance."""
    if not counts.n_valid:
        return float("nan")
    P, Q = counts.Ps, counts.Pv
    a = 1 - 2 * P - Q
    b = 1 - 2 * Q
    if a <= 0 or b <= 0:
        return float("nan")
    return -0.5 * math.log(a) - 0.25 * math.log(b)


def tn93_distance(counts: SiteCounts, freqs: np.ndarray) -> float:
    """Tamura & Nei (1993) distance.

    `freqs` are the A,C,G,T frequencies (order matching BASES) used for the
    correction; typically estimated once per alignment or codon-position
    subset and shared across all pairs.
    """
    if not counts.n_valid:
        return float("nan")
    gA, gC, gG, gT = freqs
    gR = gA + gG
    gY = gC + gT
    if gR <= 0 or gY <= 0 or gA <= 0 or gG <= 0 or gC <= 0 or gT <= 0:
        return float("nan")

    P1, P2, Q = counts.P1_purine, counts.P2_pyrimidine, counts.Pv

    k1 = 2 * gA * gG / gR
    k2 = 2 * gC * gT / gY
    k3 = 2 * (gR * gY - (gA * gG * gY) / gR - (gC * gT * gR) / gY)

    a1 = 1 - P1 / k1 - Q / (2 * gR) if k1 > 0 else float("nan")
    a2 = 1 - P2 / k2 - Q / (2 * gY) if k2 > 0 else float("nan")
    a3 = 1 - Q / (2 * gR * gY)

    if any(math.isnan(v) for v in (a1, a2, a3)) or a1 <= 0 or a2 <= 0 or a3 <= 0:
        return float("nan")

    return -k1 * math.log(a1) - k2 * math.log(a2) - k3 * math.log(a3)


# --------------------------------------------------------------------------
# GTR: aggregate-counting rate estimation + per-pair ML distance
# --------------------------------------------------------------------------


@dataclass
class GTRModel:
    freqs: np.ndarray  # order A,C,G,T
    rates: dict  # {(base, base): exchangeability}, symmetric pairs as keys
    _eigvals: np.ndarray = field(init=False, repr=False)
    _V: np.ndarray = field(init=False, repr=False)
    _Vinv: np.ndarray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        Q = self._build_rate_matrix()
        eigvals, V = np.linalg.eig(Q)
        self._eigvals = eigvals.real
        self._V = V.real
        self._Vinv = np.linalg.inv(V).real

    def _build_rate_matrix(self) -> np.ndarray:
        idx = {b: i for i, b in enumerate(BASES)}
        R = np.zeros((4, 4))
        for (x, y), r in self.rates.items():
            i, j = idx[x], idx[y]
            R[i, j] = r
            R[j, i] = r
        Q = R * self.freqs[np.newaxis, :]
        np.fill_diagonal(Q, 0.0)
        np.fill_diagonal(Q, -Q.sum(axis=1))
        mean_rate = -(self.freqs * np.diag(Q)).sum()
        if mean_rate > 0:
            Q = Q / mean_rate
        return Q

    def transition_matrix(self, t: float) -> np.ndarray:
        t = max(float(t), 0.0)
        P = (self._V @ np.diag(np.exp(self._eigvals * t)) @ self._Vinv)
        return np.clip(P, 1e-12, None)

    @classmethod
    def fit(cls, pooled_counts4x4: np.ndarray, freqs: np.ndarray) -> "GTRModel":
        """Estimate exchangeability rates from substitution counts pooled
        across every sequence pair in the alignment/subset (see module
        docstring for the rationale)."""
        idx = {b: i for i, b in enumerate(BASES)}
        sym = pooled_counts4x4 + pooled_counts4x4.T
        rates: dict = {}
        for x, y in RATE_PAIRS:
            i, j = idx[x], idx[y]
            denom = freqs[i] * freqs[j]
            rates[(x, y)] = (sym[i, j] / denom) if denom > 0 else 0.0
        ref = rates[("G", "T")]
        if ref > 0:
            rates = {k: v / ref for k, v in rates.items()}
        elif any(v > 0 for v in rates.values()):
            m = max(rates.values())
            rates = {k: v / m for k, v in rates.items()}
        else:
            rates = {k: 1.0 for k in RATE_PAIRS}
        return cls(freqs=np.asarray(freqs, dtype=float), rates=rates)

    def pairwise_ml_distance(self, pair_counts4x4: np.ndarray, t_max: float = 15.0) -> float:
        n = pair_counts4x4.sum()
        if n == 0:
            return float("nan")

        def neg_log_lik(t: float) -> float:
            P = self.transition_matrix(t)
            L = self.freqs[:, None] * P
            return -(pair_counts4x4 * np.log(L)).sum()

        result = minimize_scalar(neg_log_lik, bounds=(1e-8, t_max), method="bounded")
        if not result.success:
            return float("nan")
        if result.x >= 0.95 * t_max:
            # Likelihood still increasing at the search boundary: distance is
            # effectively unbounded (saturated) and not reliably estimable.
            return float("nan")
        return float(result.x)
