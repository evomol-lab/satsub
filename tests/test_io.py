import pytest

from satsub.io import example_alignment_path, read_alignment, translate_alignment
from tests.helpers import make_alignment


def test_reads_bundled_example_as_nucleotide_and_normalizes_rna():
    loaded = read_alignment(example_alignment_path())
    assert loaded.kind == "nucleotide"
    assert loaded.had_rna is True
    assert loaded.n_sequences == 8
    # Normalized: no 'U' should remain in the working alignment.
    for rec in loaded.alignment:
        assert "U" not in str(rec.seq)
    assert loaded.length == 1509


def test_detects_amino_acid_alignment():
    aln = make_alignment({"s1": "MKVLAGDEFH", "s2": "MKVLAGDEFY"})
    from satsub.io import detect_kind

    assert detect_kind(aln) == "amino_acid"


def test_translate_vertebrate_mitochondrial_no_internal_stops():
    loaded = read_alignment(example_alignment_path())
    aa_aln, warnings = translate_alignment(loaded.alignment, frame_start=1, table=2)
    assert len(aa_aln) == 8
    assert not any("internal stop" in w for w in warnings)
    # Translated length should be exactly one third of the nucleotide length.
    assert aa_aln.get_alignment_length() == loaded.length // 3


def test_translate_wrong_frame_can_trigger_stop_warning():
    aln = make_alignment({"s1": "TAAGGGCCC"})  # frame 1: TAA(stop) GGG CCC
    _, warnings = translate_alignment(aln, frame_start=1, table=1)
    assert any("internal stop" in w for w in warnings)
