import os
import pytest
from T2T_ACE.genomic_queries import get_sequence_from_interval


class TestGetSequenceFromInterval:
    current_dir = os.path.dirname(os.path.abspath(__file__))
    fasta_path = os.path.join(current_dir, 'mock_reference.fasta')
    non_existent_file = os.path.join(current_dir, 'non_existent_file.fasta')

    @pytest.mark.parametrize(
        "fasta_path,interval,expected", [
            (fasta_path, 'chr1:1-4', 'AGCT'),
            (fasta_path, 'chr2:1-3', 'TTG'),
            (fasta_path, 'chr3:1-4', None),
            # Add more parameter sets as needed
        ]
    )
    def test_valid_interval(self, fasta_path, interval, expected):
        """Test valid intervals."""
        assert get_sequence_from_interval(fasta_path, interval) == expected

    @pytest.mark.parametrize(
        "fasta_path,interval", [
            (fasta_path, 'chr1:1_4'),
            (fasta_path, 'chr1-1-4')
        ]
    )
    def test_invalid_interval_format(self, fasta_path, interval):
        """Test invalid interval formats."""
        with pytest.raises(TypeError):
            get_sequence_from_interval(fasta_path, interval)

    def test_invalid_file(self):
        """Test invalid FASTA file."""
        with pytest.raises(IOError):
            get_sequence_from_interval(self.non_existent_file, 'chr1:1-4')

    def test_invalid_chr(self):
        assert get_sequence_from_interval(self.fasta_path, 'chr3:1-4') is None

