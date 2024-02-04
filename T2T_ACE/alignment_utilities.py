import logging
import os

import mappy as mp
from collections import Counter
import re

from typing import List
from Bio.Align.Applications import MafftCommandline

# Initialize logging
logging.basicConfig(level=logging.INFO)
mafft_exe = "/opt/homebrew/bin/mafft"


class Config:
    # [match, mismatch, gap_open, gap_extension_penalty,
    # long_gap_open_penalty, long_gap_extension_penalty, ambiguous_match_penalty]
    # SCORING: tuple = (20, 20, 40, 20, 60, 1, 10)
    SCORING: tuple = (2, 2, 4, 2, 6, 1, 1)


class ReferenceLoadError(Exception):
    pass


def load_reference(reference: str) -> mp.Aligner:
    """
    Description:
        Load a reference genome or sequence file and create an mp.Aligner object for sequence alignment.

    Parameters:
        reference (str): The file path or URI to the reference genome or sequence data.

    Returns:
        mp.Aligner: Initialized mp.Aligner object for sequence alignment tasks.

    Exceptions:
        ReferenceLoadError: Raised when the reference genome or sequence fails to load or build.

    Examples:
        try:
            aligner = load_reference("/path/to/reference/genome")
        except ReferenceLoadError as e:
            print(e)
    """
    logging.info(f"Loading reference from: {reference}")
    aligner = mp.Aligner(reference, scoring=Config.SCORING, best_n=10, n_threads=4)
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
        if is_good_hit(hit, query_length) and 'alt' not in hit.ctg:
            print("{}{} {}: {} {}-{}\t({}-{})\t{} {} {}".format("+", name, query_length, hit.ctg, hit.r_st,
                                                             hit.r_en, hit.q_st,
                                                             hit.q_en, hit.mlen, sum_cigar_events(hit.cigar_str),
                                                             hit.strand))
        else:
            print("{} {}: {} {}-{}\t({}-{})\t{} {} {}".format(name, query_length, hit.ctg, hit.r_st,
                                                           hit.r_en, hit.q_st,
                                                           hit.q_en, hit.mlen, sum_cigar_events(hit.cigar_str),
                                                           hit.strand))


def print_hits1(alignment_hits: List[mp.Alignment]) -> None:
    """
    Description:
        Prints alignment hits in a custom format.

    Parameters:
        alignment_hits (List[mp.Alignment]): List of alignment hits as mp.Alignment objects.

    Returns:
        None: This function does not return any value.

    Example Usage:
        >>> print_hits1([hit1, hit2])
        ctg1 10-20 (5-10) 10 50
        ctg2 30-40 (10-20) 12 60
    """
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

def get_multiseq_alignment(seq_list, seq_names):
    multiseq_fasta = f"{seq_names[0]}_intervals.fasta"
    multiseq_aligned_file = f"{seq_names[0]}_intervals_mafft_aligned.fasta"

    # Remove old files if they exist
    if os.path.exists(multiseq_fasta):
        os.remove(multiseq_fasta)
    elif os.path.exists(multiseq_aligned_file):
        os.remove(multiseq_aligned_file)

    # Start new aLignment
    # Get the sequences in FASTA format
    for seq, name in zip(seq_list, seq_names):
        with open(multiseq_fasta, "a") as f:
            f.write(f">{name}\n{seq}\n")
    # Perform multiseq alignment using MAFFT
    print(f"Running MAFFT alignment for {multiseq_fasta}")
    mafft_cline = MafftCommandline(mafft_exe, input=multiseq_fasta, clustalout="on")
    mafft_cline(stdout=multiseq_aligned_file)



