import pytest

from satsub.stats_nt import codon_position_gc, summarize_nucleotide_alignment
from tests.helpers import make_alignment


def test_basic_summary_counts():
    aln = make_alignment(
        {
            "s1": "ACGTACGT",
            "s2": "ACGTACGA",  # differs from s1 at last site (T->A, transversion)
            "s3": "ACGTACGT",
        }
    )
    summary = summarize_nucleotide_alignment(aln, frame_start=None)
    assert summary.n_sequences == 3
    assert summary.alignment_length == 8
    # Only the last column varies (T,A,T) -> singleton (T appears twice, A once): not parsimony informative
    assert summary.overall["variable_sites"] == 1
    assert summary.overall["conserved_sites"] == 7
    assert summary.overall["parsimony_informative_sites"] == 0
    assert summary.overall["singleton_sites"] == 1


def test_parsimony_informative_site_detected():
    aln = make_alignment(
        {
            "s1": "AA",
            "s2": "AC",
            "s3": "CC",
            "s4": "CA",
        }
    )
    # Column 0: A,A,C,C -> 2 states each with count 2 -> parsimony informative
    # Column 1: A,C,C,A -> 2 states each with count 2 -> parsimony informative
    summary = summarize_nucleotide_alignment(aln, frame_start=None)
    assert summary.overall["parsimony_informative_sites"] == 2
    assert summary.overall["variable_sites"] == 2
    assert summary.overall["singleton_sites"] == 0


def test_gap_percent_and_identity():
    aln = make_alignment(
        {
            "s1": "ACGT",
            "s2": "AC-T",  # 1 gap, otherwise identical at valid sites
        }
    )
    summary = summarize_nucleotide_alignment(aln, frame_start=None)
    assert summary.overall["overall_gap_percent"] == pytest.approx(100 * 1 / 8)
    # identity computed over shared ungapped sites only -> 3/3 = 100%
    assert summary.identity_matrix.loc["s1", "s2"] == pytest.approx(100.0)


def test_gc_content():
    aln = make_alignment({"s1": "GGCC", "s2": "AATT"})
    summary = summarize_nucleotide_alignment(aln, frame_start=None)
    per_seq = summary.per_sequence.set_index("sequence")
    assert per_seq.loc["s1", "GC_percent"] == pytest.approx(100.0)
    assert per_seq.loc["s2", "GC_percent"] == pytest.approx(0.0)
    assert summary.overall["overall_GC_percent"] == pytest.approx(50.0)


def test_codon_position_gc():
    # 2 codons: "GGG AAA" per sequence -> pos1 all G, pos2 all G, pos3 all G for first codon
    aln = make_alignment({"s1": "GGGAAA", "s2": "GGGAAA"})
    result = codon_position_gc(aln, frame_start=1)
    # position 1 columns: 0,3 -> chars G,A -> 1 GC of 2 -> 50%
    assert result["GC1_percent"] == pytest.approx(50.0)
    assert result["GC2_percent"] == pytest.approx(50.0)
    assert result["GC3_percent"] == pytest.approx(50.0)


def test_aggregate_ts_tv_ratio():
    aln = make_alignment(
        {
            "s1": "AAAA",
            "s2": "GAAC",  # 1 transition (A->G), 1 transversion (A->C)
        }
    )
    summary = summarize_nucleotide_alignment(aln, frame_start=None)
    assert summary.overall["aggregate_transitions"] == 1
    assert summary.overall["aggregate_transversions"] == 1
    assert summary.overall["aggregate_ts_tv_ratio"] == pytest.approx(1.0)
