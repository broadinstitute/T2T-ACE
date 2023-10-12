import pytest
from T2T_ACE.interval_parsing import (parse_interval, create_interval, closest_interval, find_overlapping_intervals,
                                      sort_intervals, find_next_interval, distance_between_intervals)


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


class TestSortIntervals:
    def test_basic_sort(self):
        intervals = ["chr1:5000-6000", "chr1:1000-2000"]
        expected = [("chr1", 1000, 2000), ("chr1", 5000, 6000)]
        assert sort_intervals(intervals) == expected

    def test_same_chromosome_diff_position(self):
        intervals = ["chr1:5000-6000", "chr1:1000-2000", "chr1:7000-8000"]
        expected = [("chr1", 1000, 2000), ("chr1", 5000, 6000), ("chr1", 7000, 8000)]
        assert sort_intervals(intervals) == expected

    def test_multiple_chromosomes(self):
        intervals = ["chr2:3000-4000", "chr1:1000-2000", "chr1:5000-6000"]
        expected = [("chr1", 1000, 2000), ("chr1", 5000, 6000), ("chr2", 3000, 4000)]
        assert sort_intervals(intervals) == expected

    def test_empty_list(self):
        assert sort_intervals([]) == []

    def test_single_interval(self):
        intervals = ["chr1:1000-2000"]
        expected = [("chr1", 1000, 2000)]
        assert sort_intervals(intervals) == expected


class TestFindNextInterval:
    def setup_class(self):
        self.intervals_list = ["chr1:5000-6000", "chr1:1000-2000", "chr1:7000-8000", "chr2:3000-4000"]

    def test_immediate_next(self):
        assert find_next_interval("chr1:2342-4332", self.intervals_list) == "chr1:5000-6000"

    def test_no_next_interval(self):
        assert find_next_interval("chr1:9000-10000", self.intervals_list) is None

    def test_same_chromosome_diff_position(self):
        assert find_next_interval("chr1:4000-4500", self.intervals_list) == "chr1:5000-6000"

    def test_multiple_chromosomes(self):
        assert find_next_interval("chr2:2000-2500", self.intervals_list) == "chr2:3000-4000"

    def test_exact_match_next(self):
        assert find_next_interval("chr1:5000-6000", self.intervals_list) == "chr1:7000-8000"

    def test_empty_list(self):
        assert find_next_interval("chr1:4000-4500", []) is None

    def test_diff_chromosomes(self):
        assert find_next_interval("chr3:1000-2000", self.intervals_list) is None


class TestDistanceBetweenIntervals:

    def test_same_chr_non_overlapping(self):
        assert distance_between_intervals("chr1:100-200", "chr1:300-400") == 99

    def test_same_chr_overlapping(self):
        assert distance_between_intervals("chr1:100-250", "chr1:250-300") == -1
        assert distance_between_intervals("chr1:100-250", "chr1:200-300") == -51

    def test_same_chr_touching(self):
        assert distance_between_intervals("chr1:100-199", "chr1:200-300") == 0
        assert distance_between_intervals("chr1:200-300", "chr1:100-199") == 0

    def test_different_chr(self):
        with pytest.raises(ValueError):
            distance_between_intervals("chr1:100-200", "chr2:300-400")
