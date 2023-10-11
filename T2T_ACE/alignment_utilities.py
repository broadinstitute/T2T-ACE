import logging
import mappy as mp
from collections import Counter
import re

from typing import List

# Initialize logging
logging.basicConfig(level=logging.INFO)


class Config:
    # [match, mismatch, gap_open, gap_extension_penalty,
    # long_gap_open_penalty, long_gap_extension_penalty, ambiguous_match_penalty]
    #SCORING: tuple = (20, 20, 40, 20, 60, 1, 10)
    SCORING: tuple = (2, 2, 4, 2, 6, 1, 1)


class ReferenceLoadError(Exception):
    pass


def load_reference(reference: str) -> mp.Aligner:
    logging.info(f"Loading reference from: {reference}")
    aligner = mp.Aligner(reference, scoring=Config.SCORING)
    if not aligner:
        raise ReferenceLoadError(f"ERROR: failed to load/build {reference} reference file.")
    return aligner


def is_good_hit(hit: mp.Alignment, qlen: int, threshold: float = 0.95) -> bool:
    """
    Is the hit a good hit?
    :param hit: Alignment to check
    :param qlen: Length of the query
    :param threshold: Threshold for a good hit
    :return: True if the hit is good, False otherwise
    """
    query_match_length = hit.q_en - hit.q_st
    return query_match_length / qlen > threshold


def print_hits(name: str, query_length: int, alignment_hits: List[mp.Alignment]) -> None:
    """
    Print the hits
    :param name: Name of the sequence
    :param query_length: Length of the query sequence
    :param alignment_hits: List of alignment hits
    :return: None
    """
    for hit in alignment_hits:
        if is_good_hit(hit, query_length):
            print("{}{} {}: {} {}-{}\t({}-{})\t{} {}".format("+", name, query_length, hit.ctg, hit.r_st,
                                                       hit.r_en, hit.q_st,
                                                       hit.q_en, hit.mlen, sum_cigar_events(hit.cigar_str)))
        else:
            print("{} {}: {} {}-{}\t({}-{})\t{} {}".format(name, query_length, hit.ctg, hit.r_st,
                                                       hit.r_en, hit.q_st,
                                                       hit.q_en, hit.mlen, sum_cigar_events(hit.cigar_str)))


def print_hits1(alignment_hits: List[mp.Alignment]) -> None:
    for hit in alignment_hits:
        print("{} {}-{}\t({}-{})\t{} {}".format(hit.ctg, hit.r_st, hit.r_en, hit.q_st, hit.q_en, hit.mlen,
                                                    sum_cigar_events(hit.cigar_str)))


# Precompiled regex pattern
CIGAR_PATTERN = re.compile(r'(\d+)([MID])')


def sum_cigar_events(cigar_str: str) -> str:
    """
    Summarize the events in a CIGAR string used for genomic alignments.

    Parameters:
    cigar_string (str): The CIGAR string representing alignment events. E.g., "10M3I5D".

    Returns:
    str: A string summarizing the counts of each event type (e.g., "M:10, I:3, D:5").

    Raises:
    ValueError: If the CIGAR string contains invalid characters or is empty.
    TypeError: If the input is not a string.

    Example:
    sum_cigar_events("10M3I5D")  # Returns "10M3I5D
    sum_cigar_events("Z10")  # Raises ValueError
    sum_cigar_events(100)  # Raises TypeError
    """

    counts = Counter()
    found_pairs = CIGAR_PATTERN.findall(cigar_str)

    # Validate empty or invalid string
    if not cigar_str or not found_pairs:
        raise ValueError(f"Empty or malformed CIGAR string: {cigar_str}")

    for num, op in found_pairs:
        counts[op] += int(num)

    # Validate if entire string was parsed
    parsed_str = "".join(f"{num}{op}" for num, op in found_pairs)
    if parsed_str != cigar_str:
        raise ValueError(f"Unparsed segments in CIGAR string: {cigar_str}")

    # Validate if counts are filled correctly
    if not all(op in 'MID' for op in counts.keys()):
        raise ValueError(f"Invalid operations in CIGAR string: {cigar_str}")

    # Create the aggregated CIGAR string
    aggregated_cigar_str = "".join(f"{counts[op]}{op}" for op in 'MID' if counts.get(op, 0) > 0)

    return aggregated_cigar_str

# %%