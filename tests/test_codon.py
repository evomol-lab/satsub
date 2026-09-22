import pytest

from satsub.codon import codon_position_columns, columns_for_subset


def test_codon_columns_frame_start_1_basic():
    cols = codon_position_columns(9, frame_start=1)
    assert cols[1] == [0, 3, 6]
    assert cols[2] == [1, 4, 7]
    assert cols[3] == [2, 5, 8]


def test_codon_columns_drops_incomplete_trailing_codon():
    cols = codon_position_columns(10, frame_start=1)  # 10 % 3 == 1 leftover column, dropped entirely
    assert cols[1] == [0, 3, 6]
    assert cols[2] == [1, 4, 7]
    assert cols[3] == [2, 5, 8]


def test_codon_columns_frame_start_2():
    cols = codon_position_columns(10, frame_start=2)
    # First usable column is index 1 (0-based); usable length = ((10-1)//3)*3 = 9
    assert cols[1] == [1, 4, 7]
    assert cols[2] == [2, 5, 8]
    assert cols[3] == [3, 6, 9]


def test_invalid_frame_start_raises():
    with pytest.raises(ValueError):
        codon_position_columns(9, frame_start=4)


def test_columns_for_subset_all_returns_none():
    assert columns_for_subset(9, "all", frame_start=1) is None


def test_columns_for_subset_matches_codon_position_columns():
    assert columns_for_subset(9, "2", frame_start=1) == [1, 4, 7]
