import os
import pytest

from T2T_ACE.genomic_queries import get_sequence_from_interval, get_flanking_regions


class TestGetSequenceFromInterval:
    current_dir = os.path.dirname(os.path.abspath(__file__))
    fasta_path = os.path.join(current_dir, 'mock_reference.fasta')
    non_existent_file = os.path.join(current_dir, 'non_existent_file.fasta')

    @pytest.mark.parametrize(
        "fasta_path,interval,expected", [
            (fasta_path, 'chr1:1-4', 'AGCT'),
            (fasta_path, 'chr2:1-3', 'TTG'),
            # Add more parameter sets as needed
        ]
    )
    def test_valid_interval(self, fasta_path, interval, expected):
        """Test valid intervals."""
        assert get_sequence_from_interval(self.fasta_path, interval) == expected

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
        with pytest.raises(KeyError):
            get_sequence_from_interval(self.fasta_path, 'chr3:1-4')

    def test_invalid_interval(self):
        with pytest.raises(IndexError):
            assert get_sequence_from_interval(self.fasta_path, 'chr1:100-400')


class TestGetFlankingRegions:

    current_dir = os.path.dirname(os.path.abspath(__file__))
    fasta_path = os.path.join(current_dir, 'mock_reference.fasta')
    non_existent_file = os.path.join(current_dir, 'non_existent_file.fasta')

    def test_valid_intervals(self):
        left, right = get_flanking_regions(self.fasta_path, "chr1:4-8", padding=2)
        assert left == "GC"
        assert right == "AG"

    def test_valid_intervals_no_padding(self):
        left, right = get_flanking_regions(self.fasta_path, "chr1:4-8", padding=0)
        assert left == ""
        assert right == ""

    def test_invalid_intervals_same_chr(self):
        with pytest.raises(IndexError):
            get_flanking_regions(self.fasta_path, "chr1:20-25", padding=2)

    def test_invalid_intervals_diff_chr(self):
        with pytest.raises(TypeError): # Should be KeyError
            get_flanking_regions(self.fasta_path, "chr3:4-8")

    def test_invalid_padding(self):
        with pytest.raises(ValueError):
            get_flanking_regions(self.fasta_path, "chr1:4-8", padding=-1)

    def test_file_not_found(self): # Should be IOError
        with pytest.raises(OSError):
            get_flanking_regions("nonexistent.fasta", "chr1:4-8", padding=2)
