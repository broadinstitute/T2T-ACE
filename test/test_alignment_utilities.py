import pytest
from T2T_ACE.alignment_utilities import load_reference, ReferenceLoadError, sum_cigar_events


# Test load_reference function
def test_load_non_existent_reference():
    with pytest.raises(ReferenceLoadError):
        load_reference("non_existent_file.fasta")


class TestSumCigarEvents:
    @pytest.mark.parametrize("cigar_str,expected_output", [
        ("10M3I5D", "10M3I5D"),
        ("20M", "20M"),
        ("10M3I2M3I4D1M", "13M6I4D"),
        ("10I5D", "10I5D"),
        ("10M10M", "20M"),  # Repeated operations
    ])
    def test_valid_cigar_strings(self, cigar_str, expected_output):
        assert sum_cigar_events(cigar_str) == expected_output

    @pytest.mark.parametrize("cigar_str", [
        "",  # Empty string
        "M10",  # No digit before operation
        "10",  # No operation after digit
        "MID",  # Only operations, no digits
    ])
    def test_invalid_cigar_strings(self, cigar_str):
        with pytest.raises(ValueError):
            sum_cigar_events(cigar_str)

    @pytest.mark.parametrize("cigar_str", [
        None,
        123,
        []
    ])
    def test_invalid_input_types(self, cigar_str):
        with pytest.raises(TypeError):
            sum_cigar_events(cigar_str)




