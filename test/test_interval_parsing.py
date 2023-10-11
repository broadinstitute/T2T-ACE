import pytest
from T2T_ACE.interval_parsing import parse_interval, create_interval, closest_interval, find_overlapping_intervals


class TestParseInterval:
    @pytest.mark.parametrize("interval,expected", [
        ("chr1:200-300", ("chr1", 200, 300)),
        ("chrX:1000-2000", ("chrX", 1000, 2000)),
        ("chrY:50-100", ("chrY", 50, 100))
    ])
    def test_valid_intervals(self, interval: str, expected: tuple):
        assert parse_interval(interval) == expected

    def test_invalid_interval_format(self):
        assert parse_interval("invalid") is None
        assert parse_interval("chr1:200-") is None
        assert parse_interval("chr1:200") is None

    def test_empty_chromosome(self):
        with pytest.raises(ValueError):
            assert parse_interval(":200-300")

    def test_interval_with_negative_numbers(self):
        assert parse_interval("chr1:-200-300") is None
        assert parse_interval("chr1:200--300") is None


class TestCreateInterval:
    def test_valid_interval(self):
        assert create_interval("chr1", 200, 300) == "chr1:200-300"

    def test_negative_start_position(self):
        with pytest.raises(ValueError):
            create_interval("chr1", -1, 300)

    def test_negative_end_position(self):
        with pytest.raises(ValueError):
            create_interval("chr1", 200, -1)

    def test_both_positions_negative(self):
        with pytest.raises(ValueError):
            create_interval("chr1", -1, -1)


class TestFindOverlappingIntervals:
    def test_simple_overlap(self):
        assert (find_overlapping_intervals(["chr1:200-300", "chr1:250-350"]) ==
                {"chr1:200-300": ["chr1:250-350"], "chr1:250-350": ["chr1:200-300"]})

    def test_no_overlap(self):
        assert find_overlapping_intervals(["chr1:200-300", "chr1:350-450"]) == {}

    def test_multiple_overlaps(self):
        intervals = ["chr1:200-300", "chr1:250-350", "chr1:325-350"]
        expected = {
            "chr1:200-300": ["chr1:250-350"],
            "chr1:250-350": ["chr1:200-300", "chr1:325-350"],
            "chr1:325-350": ["chr1:250-350"]
        }
        assert find_overlapping_intervals(intervals) == expected

    def test_different_chromosomes(self):
        assert find_overlapping_intervals(["chr1:200-300", "chr2:250-350"]) == {}

    def test_empty_list(self):
        assert find_overlapping_intervals([]) == {}

    def test_invalid_format(self):
        with pytest.raises(IndexError):
            find_overlapping_intervals(["invalid"])

    def test_single_interval(self):
        assert find_overlapping_intervals(["chr1:200-300"]) == {}


class TestClosestInterval:
    def test_same_chromosome(self):
        assert closest_interval("chr1:200-300", ["chr1:100-200", "chr1:250-350"]) == ["chr1:250-350"]

    def test_multiple_closest_intervals(self):
        assert closest_interval("chr1:200-300",
                                ["chr1:100-200", "chr1:210-310", "chr1:210-320"]) == ["chr1:210-310", "chr1:210-320"]

    def test_no_matching_chromosome(self):
        assert closest_interval("chr1:200-300", ["chr2:100-200", "chr3:250-350"]) == []

    def test_empty_interval_list(self):
        assert closest_interval("chr1:200-300", []) == []

    def test_invalid_target_format(self):
        with pytest.raises(ValueError):
            closest_interval("invalid", ["chr1:100-200"])

    def test_invalid_interval_list_format(self):
        with pytest.raises(ValueError):
            closest_interval("chr1:200-300", ["invalid"])

    def test_single_interval_same_chromosome(self):
        assert closest_interval("chr1:200-300", ["chr1:100-200"]) == ["chr1:100-200"]

    def test_single_interval_different_chromosome(self):
        assert closest_interval("chr1:200-300", ["chr2:100-200"]) == []