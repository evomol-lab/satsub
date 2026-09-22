import math

import numpy as np
import pytest

from satsub.distances import (
    GTRModel,
    encode_sequence,
    jc69_distance,
    k80_distance,
    pairwise_site_counts,
    tn93_distance,
)


def counts_from_pair(seq1: str, seq2: str):
    return pairwise_site_counts(encode_sequence(seq1), encode_sequence(seq2))


def test_identical_sequences_have_zero_distance():
    sc = counts_from_pair("ACGTACGT", "ACGTACGT")
    assert sc.p_distance == 0
    assert jc69_distance(sc) == pytest.approx(0.0)
    assert k80_distance(sc) == pytest.approx(0.0)
    assert tn93_distance(sc, np.full(4, 0.25)) == pytest.approx(0.0, abs=1e-9)


def test_gaps_and_ambiguous_sites_excluded():
    sc = counts_from_pair("ACGTN-AC", "ACGTAAAG")
    # positions: A-A id, C-C id, G-G id, T-T id, N-A excluded, ---A excluded, A-A id, C-G tv
    assert sc.n_valid == 6
    assert sc.n_excluded == 2
    assert sc.n_tv == 1


def test_jc69_matches_hand_calculation():
    # 1 transition out of 4 valid sites -> p = 0.25
    sc = counts_from_pair("AAAA", "AAAG")
    expected = -0.75 * math.log(1 - (4 / 3) * 0.25)
    assert jc69_distance(sc) == pytest.approx(expected)


def test_jc69_undefined_when_fully_saturated():
    # p >= 0.75 makes the JC69 log argument non-positive.
    sc = counts_from_pair("ACGTACGT", "TGCATGCA")  # all sites differ, p = 1.0
    assert math.isnan(jc69_distance(sc))


def test_k80_matches_hand_calculation():
    sc = counts_from_pair("AAAA", "AAAG")  # 1 transition, 0 transversions, n=4
    P, Q = 0.25, 0.0
    expected = -0.5 * math.log(1 - 2 * P - Q) - 0.25 * math.log(1 - 2 * Q)
    assert k80_distance(sc) == pytest.approx(expected)


def test_k80_pure_transversions_only_uses_second_term():
    sc = counts_from_pair("AAAA", "CCCC")  # transversions only
    P, Q = 0.0, 1.0
    # Q = 1.0 makes the second log argument (1 - 2Q) negative -> undefined.
    assert math.isnan(k80_distance(sc))


def test_tn93_reduces_to_k80_with_equal_frequencies_and_balanced_transitions():
    # Construct a pair with equal purine and pyrimidine transition counts and
    # equal base frequencies so TN93 must collapse exactly onto K80.
    seq1 = "AAAACCCC"
    seq2 = "AGAACCTC"  # 1 A<->G, 1 C<->T -> P1 == P2, low enough divergence to stay defined
    sc = counts_from_pair(seq1, seq2)
    assert sc.P1_purine == pytest.approx(sc.P2_pyrimidine)
    d_tn93 = tn93_distance(sc, np.full(4, 0.25))
    d_k80 = k80_distance(sc)
    assert d_tn93 == pytest.approx(d_k80, rel=1e-9)


def test_tn93_undefined_returns_nan_not_exception():
    sc = counts_from_pair("ACGTACGT", "TGCATGCA")
    d = tn93_distance(sc, np.full(4, 0.25))
    assert math.isnan(d)


def test_gtr_transition_matrix_is_row_stochastic():
    model = GTRModel.fit(
        pooled_counts4x4=np.array(
            [[0, 3, 5, 2], [3, 0, 2, 6], [5, 2, 0, 3], [2, 6, 3, 0]], dtype=np.int64
        ),
        freqs=np.array([0.3, 0.2, 0.2, 0.3]),
    )
    for t in (0.0, 0.1, 1.0, 5.0):
        P = model.transition_matrix(t)
        row_sums = P.sum(axis=1)
        assert row_sums == pytest.approx(np.ones(4), abs=1e-6)
    # At t=0 the transition matrix must be the identity.
    P0 = model.transition_matrix(0.0)
    assert P0 == pytest.approx(np.eye(4), abs=1e-6)


def test_gtr_equal_rates_equal_freqs_reduces_to_jc69_style_matrix():
    freqs = np.full(4, 0.25)
    rates = {pair: 1.0 for pair in [("A", "C"), ("A", "G"), ("A", "T"), ("C", "G"), ("C", "T"), ("G", "T")]}
    model = GTRModel(freqs=freqs, rates=rates)
    P = model.transition_matrix(0.5)
    # Under equal rates & frequencies every off-diagonal entry in a row must
    # be equal, by symmetry (this is exactly the JC69 substitution model).
    for row in P:
        off_diag = [row[k] for k in range(4)]
        diag = row.max()
        others = sorted(off_diag)[:3]
        assert max(others) - min(others) < 1e-9


def test_gtr_ml_distance_recovers_known_branch_length():
    # Build an "infinite-data" expected count matrix for a known true distance
    # under the model itself, then check the MLE recovers that distance.
    freqs = np.array([0.3, 0.2, 0.2, 0.3])
    rates = {("A", "C"): 1.0, ("A", "G"): 4.0, ("A", "T"): 1.0, ("C", "G"): 1.0, ("C", "T"): 4.0, ("G", "T"): 1.0}
    model = GTRModel(freqs=freqs, rates=rates)
    true_t = 0.35
    P = model.transition_matrix(true_t)
    n = 100_000
    expected_counts = np.outer(freqs, np.ones(4)) * P * n
    est_t = model.pairwise_ml_distance(expected_counts)
    assert est_t == pytest.approx(true_t, rel=1e-2)


def test_gtr_fit_rate_estimation_recovers_relative_rates_from_pooled_counts():
    # Pooled counts proportional to freq_i * freq_j * rate_ij should recover
    # the rate ratios used to generate them.
    idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    freqs = np.array([0.25, 0.25, 0.25, 0.25])
    true_rates = {("A", "C"): 1.0, ("A", "G"): 5.0, ("A", "T"): 1.0, ("C", "G"): 1.0, ("C", "T"): 5.0, ("G", "T"): 1.0}
    pooled = np.zeros((4, 4))
    for (x, y), r in true_rates.items():
        i, j = idx[x], idx[y]
        val = r * freqs[i] * freqs[j] * 1000
        pooled[i, j] = val
        pooled[j, i] = val
    model = GTRModel.fit(pooled, freqs)
    # G<->T was used as the normalization reference (=1); A<->G and C<->T should be ~5x that.
    assert model.rates[("A", "G")] == pytest.approx(5.0, rel=1e-2)
    assert model.rates[("C", "T")] == pytest.approx(5.0, rel=1e-2)
    assert model.rates[("A", "C")] == pytest.approx(1.0, rel=1e-2)
