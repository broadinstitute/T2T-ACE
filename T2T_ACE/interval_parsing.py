from typing import Optional, Tuple, List, Dict
from collections import defaultdict

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
