import pytest

from satsub.stats_aa import summarize_amino_acid_alignment
from tests.helpers import make_alignment


def test_basic_aa_summary():
    aln = make_alignment(
        {
            "s1": "MKV-LAG",
            "s2": "MKVWLAG",
        }
    )
    summary = summarize_amino_acid_alignment(aln)
    assert summary.n_sequences == 2
    assert summary.alignment_length == 7
    # col index 3: s1='-', s2='W' -> only one non-gap letter ('W') present -> counts as conserved (1 distinct state)
    assert summary.overall["conserved_sites"] == 7
    assert summary.overall["variable_sites"] == 0
    assert summary.overall["overall_gap_percent"] == pytest.approx(100 * 1 / 14)


def test_identity_and_similarity():
    aln = make_alignment(
        {
            "s1": "MKVL",
            "s2": "MKIL",  # V->I is a conservative substitution (positive BLOSUM62 score)
        }
    )
    summary = summarize_amino_acid_alignment(aln)
    assert summary.identity_matrix.loc["s1", "s2"] == pytest.approx(75.0)  # 3/4 identical
    # V<->I has a positive BLOSUM62 score, so similarity counts it too: 4/4
    assert summary.similarity_matrix.loc["s1", "s2"] == pytest.approx(100.0)


def test_physicochemical_grouping_sums_to_100():
    aln = make_alignment({"s1": "ADKSCW"})  # one of each group: hydrophobic A/W, acidic D, basic K, polar S/C
    summary = summarize_amino_acid_alignment(aln)
    row = summary.per_sequence.iloc[0]
    total = (
        row["hydrophobic_percent"]
        + row["polar_uncharged_percent"]
        + row["acidic_percent"]
        + row["basic_percent"]
    )
    assert total == pytest.approx(100.0)


def test_parsimony_informative_amino_acid_site():
    aln = make_alignment(
        {
            "s1": "A",
            "s2": "A",
            "s3": "D",
            "s4": "D",
        }
    )
    summary = summarize_amino_acid_alignment(aln)
    assert summary.overall["parsimony_informative_sites"] == 1
    assert summary.overall["variable_sites"] == 1
