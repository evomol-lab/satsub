import pytest

from satsub.composition import codon_usage_table, directional_pair_frequencies, nucleotide_composition_table
from satsub.io import example_alignment_path, read_alignment


@pytest.fixture(scope="module")
def vertcoi_alignment():
    return read_alignment(example_alignment_path()).alignment


def test_nucleotide_composition_matches_reference_avg_row(vertcoi_alignment):
    df = nucleotide_composition_table(vertcoi_alignment, frame_start=1)
    avg = df.set_index("sequence").loc["Avg."]
    # Reference values from VertCOI_nucleotide-composition.txt (DAMBE-style report)
    assert avg["T_percent"] == pytest.approx(27.9, abs=0.1)
    assert avg["C_percent"] == pytest.approx(28.4, abs=0.1)
    assert avg["A_percent"] == pytest.approx(27.1, abs=0.1)
    assert avg["G_percent"] == pytest.approx(16.6, abs=0.1)
    assert avg["Total_sites"] == 1509
    assert avg["T3_percent"] == pytest.approx(20.4, abs=0.1)
    assert avg["C3_percent"] == pytest.approx(37.3, abs=0.1)
    assert avg["A3_percent"] == pytest.approx(36.6, abs=0.1)
    assert avg["G3_percent"] == pytest.approx(5.8, abs=0.1)
    assert avg["Pos3_sites"] == 503


def test_nucleotide_composition_per_sequence_row(vertcoi_alignment):
    df = nucleotide_composition_table(vertcoi_alignment, frame_start=1)
    homo = df.set_index("sequence").loc["HomoSapiensCOX1"]
    assert homo["T_percent"] == pytest.approx(26.7, abs=0.1)
    assert homo["G1_percent"] == pytest.approx(28.6, abs=0.1)


def test_codon_usage_matches_reference_counts_and_rscu(vertcoi_alignment):
    result = codon_usage_table(vertcoi_alignment, frame_start=1, table_id=2)
    t = result.table.set_index("codon")
    # Reference values from VertCOI_codon-usage.txt (RNA notation, table 2)
    assert t.loc["UUU", "count_avg"] == pytest.approx(13.4, abs=0.05)
    assert t.loc["UUU", "amino_acid"] == "F"
    assert t.loc["UUU", "rscu"] == pytest.approx(0.65, abs=0.01)
    assert t.loc["UUC", "count_avg"] == pytest.approx(27.6, abs=0.05)
    assert t.loc["UUC", "rscu"] == pytest.approx(1.35, abs=0.01)
    assert t.loc["CUA", "count_avg"] == pytest.approx(26.0, abs=0.05)
    assert t.loc["CUA", "rscu"] == pytest.approx(2.55, abs=0.01)
    assert t.loc["UGA", "amino_acid"] == "W"  # vertebrate mitochondrial: UGA -> Trp, not stop
    assert t.loc["UGA", "rscu"] == pytest.approx(1.77, abs=0.01)
    assert t.loc["AGA", "amino_acid"] == "*"  # vertebrate mitochondrial: AGA/AGG are stops
    assert result.avg_codons_per_sequence == pytest.approx(503.0, abs=0.01)
    assert result.table["count_avg"].sum() == pytest.approx(503.4, abs=1.0)


def test_directional_pair_frequencies_match_reference(vertcoi_alignment):
    result = directional_pair_frequencies(vertcoi_alignment, frame_start=1)
    s = result.summary.set_index("subset")
    # Reference values from VertCOI-Directional-Pair-Freqs.txt
    assert s.loc["all", "ii"] == 1180
    assert s.loc["all", "si"] == 173
    assert s.loc["all", "sv"] == 156
    assert s.loc["all", "R"] == pytest.approx(1.1)
    assert s.loc["1", "ii"] == 446
    assert s.loc["1", "si"] == 36
    assert s.loc["1", "sv"] == 21
    assert s.loc["2", "R"] == pytest.approx(1.4)
    assert s.loc["3", "ii"] == 246
    assert s.loc["3", "si"] == 129
    assert s.loc["3", "sv"] == 128

    m = result.matrices["all"]
    # The file's directional convention is transposed relative to ours (an
    # arbitrary, order-dependent choice with no canonical direction when
    # pooling over unordered pairs); diagonal (identical-pair) entries and
    # the transpose of the off-diagonal entries must match the reference.
    assert m.loc["T", "T"] == pytest.approx(331, abs=1)
    assert m.loc["C", "T"] == pytest.approx(54, abs=1)  # reference "TC"
    assert m.loc["T", "C"] == pytest.approx(68, abs=1)  # reference "CT"
    assert m.values.sum() == pytest.approx(1509, abs=1)


def test_directional_pair_frequencies_subset_sites():
    from satsub.io import example_alignment_path, read_alignment

    aln = read_alignment(example_alignment_path()).alignment
    result = directional_pair_frequencies(aln, frame_start=1)
    s = result.summary.set_index("subset")
    assert s.loc["all", "n_sites"] == 1509
    assert s.loc["1", "n_sites"] == 503
    assert s.loc["2", "n_sites"] == 503
    assert s.loc["3", "n_sites"] == 503
