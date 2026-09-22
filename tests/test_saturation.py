import math

from satsub.saturation import compute_saturation_table
from tests.helpers import make_alignment


def test_saturation_table_shape_and_columns():
    aln = make_alignment(
        {
            "s1": "ACGTACGTACGTACGT",
            "s2": "ACGAACGTACGCACGT",
            "s3": "AGGTACGTACTTACGT",
        }
    )
    result = compute_saturation_table(aln, subset="all", frame_start=1)
    df = result.table
    assert len(df) == 3  # 3 choose 2 pairs
    for col in ("seq1", "seq2", "n_sites", "p_distance", "Ps", "Pv", "jc69", "k80", "tn93", "gtr"):
        assert col in df.columns


def test_codon_position_subset_uses_fewer_sites():
    aln = make_alignment(
        {
            "s1": "ACGTACGTACGTACGT",
            "s2": "ACGAACGTACGCACGT",
        }
    )
    full = compute_saturation_table(aln, subset="all", frame_start=1)
    pos1 = compute_saturation_table(aln, subset="1", frame_start=1)
    assert pos1.n_sites_used < full.n_sites_used


def test_identical_sequences_give_zero_distance_under_every_model():
    aln = make_alignment({"s1": "ACGTACGTACGTACGT" * 3, "s2": "ACGTACGTACGTACGT" * 3})
    result = compute_saturation_table(aln, subset="all", frame_start=1)
    row = result.table.iloc[0]
    for model in ("jc69", "k80", "tn93", "gtr"):
        # Closed-form models are exact; the GTR branch-length optimizer has a
        # small numerical tolerance near its lower search bound.
        assert math.isclose(row[model], 0.0, abs_tol=1e-3)
