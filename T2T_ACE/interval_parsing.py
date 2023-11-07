from typing import Optional, Tuple, List, Dict
from collections import defaultdict
from bisect import bisect_right

import logging


def parse_interval(interval: str) -> Optional[Tuple[str, int, int]]:
    """
    Parse a genomic interval string into its components: chromosome, start, and end positions.

    Parameters:
    interval (str): Genomic interval in the format "contig:start-end". E.g., "chr1:200-300".

    Returns:
    Optional[Tuple[str, int, int]]: A tuple containing contig (str), start (int), and end (int) positions.
                                    Returns None if the interval is not valid.

    Raises:
    ValueError: If the input interval string is malformed.

    Example:
    parse_interval("chr1:200-300")  # Returns ("chr1", 200, 300)
    parse_interval("invalid")       # Raises ValueError
    """
    if interval == "":
        logging.error("Interval is empty", extra={'interval': interval})
        raise ValueError("Interval cannot be empty")

    # Split the interval into chromosome, start, and end
    try:
        chrom, positions = interval.split(":")
        start, end = map(int, positions.split("-"))
    except (ValueError, AttributeError):  # Catching errors related to splitting and conversion to int
        logging.exception("Error in parsing interval", extra={'interval': interval})
        return None

    # Validate the parsed values
    if chrom == "":
        logging.error("Chromosome name is empty", extra={'interval': interval})
        raise ValueError("Chromosome name cannot be empty")

    if start < 0 or end < 0:
        logging.error("Start or end position is negative", extra={'interval': interval})
        raise ValueError("Start and end positions must be non-negative")

    return chrom, start, end


def create_interval(chrom: str, start: int, end: int) -> str:
    """
    Create a genomic interval string from its components: chromosome, start, and end positions.

    Parameters:
    chrom (str): Chromosome name. E.g., "chr1".
    start (int): Start position of the interval. E.g., 200.
    end (int): End position of the interval. E.g., 300.

    Returns:
    str: Genomic interval in the format "chr:start-end". E.g., "chr1:200-300".

    Raises:
    ValueError: If the start or end positions are negative.

    Example:
    create_interval("chr1", 200, 300)  # Returns "chr1:200-300"
    create_interval("chr1", -1, 300)   # Raises ValueError
    """
    if start < 0 or end < 0:
        raise ValueError("Start and end positions must be non-negative")

    return f"{chrom}:{start}-{end}"


def find_overlapping_intervals(interval_list: List[str]) -> Dict[str, List[str]]:
    """
     Finds overlapping genomic intervals in a list and returns them as a dictionary.

     This function sorts the given list of genomic intervals and identifies overlapping
     intervals. Overlaps are determined based on the start and end positions within the
     same chromosome. The function returns a dictionary where keys are intervals that
     have at least one overlapping interval, and the corresponding value is a list of
     overlapping intervals.

     Parameters:
     -----------
     interval_list : List[str]
         A list of genomic intervals in the format "chr:start-end".
         E.g., ["chr1:200-300", "chr1:250-350", "chr2:400-500"].

     Returns:
     --------
     Dict[str, List[str]]
         A dictionary where each key is a genomic interval that has at least one overlap,
         and the corresponding value is a list of overlapping intervals.
         E.g., {"chr1:200-300": ["chr1:250-350"], "chr1:250-350": ["chr1:200-300"]}

     Example:
     --------
     >>> find_overlapping_intervals(["chr1:200-300", "chr1:250-350", "chr2:400-500"])
     {"chr1:200-300": ["chr1:250-350"], "chr1:250-350": ["chr1:200-300"]}
     """
    sorted_intervals = sorted(
        interval_list,
        key=lambda x: (x.split(":")[0], int(x.split(":")[1].split("-")[0]))
    )
    overlapping_dict = defaultdict(list)

    for i, interval_1 in enumerate(sorted_intervals[:-1]):
        chr1, range1 = interval_1.split(":")
        start1, end1 = map(int, range1.split("-"))

        for interval_2 in sorted_intervals[i + 1:]:
            chr2, range2 = interval_2.split(":")
            start2, end2 = map(int, range2.split("-"))

            if chr1 != chr2:
                break

            if start2 <= end1:
                overlapping_dict[interval_1].append(interval_2)
                overlapping_dict[interval_2].append(interval_1)

    return overlapping_dict


def closest_interval(target: str, interval_list: List[str]) -> List[str]:
    """
    Finds the genomic interval(s) in the list closest to the target interval based on start positions.

    Parameters:
    -----------
    target : str
        The target genomic interval in the format "chr:start-end".
    interval_list : List[str]
        List of genomic intervals to compare against the target.

    Returns:
    --------
    List[str]
        List of closest genomic interval(s). Returns an empty list if no intervals are on the same chromosome.

    Example:
    --------
    >>> closest_interval("chr1:200-300", ["chr1:100-200", "chr1:210-310"])
    ["chr1:210-310"]
    """
    closest_intervals = []
    min_distance = float('inf')

    target_chr, target_range = target.split(":")
    target_start, target_end = map(int, target_range.split("-"))

    for interval in interval_list:
        chr, range = interval.split(":")
        start, end = map(int, range.split("-"))

        if chr != target_chr:
            continue

        distance = abs(target_start - start)

        if distance < min_distance:
            closest_intervals = [interval]
            min_distance = distance
        elif distance == min_distance:
            closest_intervals.append(interval)

    return closest_intervals


def sort_intervals(intervals):
    """
    Sort a list of genomic intervals.

    Args:
        intervals (list): List of genomic intervals in string format.

    Returns:
        list: A list of sorted genomic intervals in tuple format.
    """
    parsed_intervals = [parse_interval(interval) for interval in intervals]
    return sorted(parsed_intervals, key=lambda x: (x[0], x[1]))


def find_next_interval(target, intervals):
    """
    Find the genomic interval immediately to the right of the target interval, if it exists.

    The function first filters intervals based on the chromosome of the target, sorts them,
    and then finds the interval that is immediately to the right of the target interval.

    Args:
        target (str): The target genomic interval in the format "chrX:start-end".
        intervals (list): A list of genomic intervals in string format.

    Returns:
        str or None: The genomic interval immediately to the right of the target, or None if no such interval exists.

    Example:
        >>> find_next_interval("chr1:2342-4332", ["chr1:5000-6000", "chr1:1000-2000", "chr1:7000-8000"])
        'chr1:5000-6000'
    """
    parsed_target = parse_interval(target)
    target_chrom = parsed_target[0]

    # Filter intervals by chromosome
    filtered_intervals = [x for x in sort_intervals(intervals) if x[0] == target_chrom]

    index = bisect_right(filtered_intervals, parsed_target)

    if index >= len(filtered_intervals):
        return None

    next_interval = filtered_intervals[index]
    return f"{next_interval[0]}:{next_interval[1]}-{next_interval[2]}"


def distance_between_intervals(interval1: str, interval2: str) -> int:
    """
    Find the distance between two genomic intervals.

    The distance between two intervals is defined as follows:
    - For non-overlapping intervals, it's the minimum distance between any two positions in the two intervals.
    - For overlapping intervals, the distance is negative and indicates the size of the overlap.

    Args:
        interval1 (str): The first genomic interval in the format "chrX:start-end".
        interval2 (str): The second genomic interval in the format "chrX:start-end".

    Returns:
        int: The distance between the two intervals. Negative if overlapping.

    Raises:
        ValueError: If the intervals are not on the same chromosome.

    Example:
        >>> distance_between_intervals("chr1:100-200", "chr1:300-400")
        99
        >>> distance_between_intervals("chr1:100-200", "chr1:200-300")
        -1
        yg_update: I think this is wrong. Should be 0
    """
    if interval1 is None or interval2 is None:
        return 0

    #parse_interval(sorted([interval1,interval2], key=lambda x: parse_interval(x)[1]))
    parsed_interval1 = parse_interval(interval1)
    parsed_interval2 = parse_interval(interval2)

    if parsed_interval1[0] != parsed_interval2[0]:
        raise ValueError("Intervals must be on the same chromosome")

    start1, end1 = parsed_interval1[1], parsed_interval1[2]
    start2, end2 = parsed_interval2[1], parsed_interval2[2]

    if start1 > end2:
        return start1 - end2 - 1
    elif start2 > end1:
        return start2 - end1 - 1
    else:
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)
        return overlap_start - overlap_end - 1  # Negative value for overlap


def interval_size(interval: str) -> int:
    """
    Returns the size of the interval.

    Args:
        interval (str): The genomic interval in the format "chrX:start-end".

    Returns:
        int: The size of the interval.
    """
    chrom, start, end = parse_interval(interval)
    return end - start + 1


def interval_between_intervals(interval1: str, interval2: str) -> str:
    """
    Constructs the genomic interval that lies between two non-overlapping intervals.

    Args:
        interval1 (str): The first genomic interval in the format "chrX:start-end".
        interval2 (str): The second genomic interval in the format "chrX:start-end".

    Returns:
        str: The genomic interval that lies between the two given intervals.

    Raises:
        ValueError: If the intervals are not on the same chromosome or are overlapping.

    Example:
        >>> interval_between_intervals("chr1:100-200", "chr1:300-400")
        "chr1:201-299"
    """
    chrom1, start1, end1 = parse_interval(interval1)
    chrom2, start2, end2 = parse_interval(interval2)

    if chrom1 != chrom2:
        raise ValueError("Intervals must be on the same chromosome")

    chr_name = chrom1

    if start1 - end2 == 1 or start2 - end1 == 1:
        raise ValueError("Intervals should not be touching")

    if start1 > end2:
        return f"{chr_name}:{end2 + 1}-{start1 - 1}"
    elif start2 > end1:
        return f"{chr_name}:{end1 + 1}-{start2 - 1}"
    else:
        raise ValueError("Intervals should not be overlapping")

