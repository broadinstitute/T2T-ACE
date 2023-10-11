from typing import Optional, Tuple
from interval_parsing import parse_interval
import pysam
import logging


def get_sequence_from_interval(reference_fasta: str, interval: str) -> Optional[str]:
    """
    Retrieve a DNA sequence for a given genomic interval.

    Parameters:
    chromosome (str): The chromosome where the interval is located. E.g., "chr1".
    start (int): The start position of the interval. E.g., 200.
    end (int): The end position of the interval. E.g., 300.

    Returns:
    str: The DNA sequence corresponding to the genomic interval.

    Raises:
    ValueError: If the start or end positions are negative.
    IndexError: If the interval is out of bounds for the given chromosome.
    IOError: If there's a problem reading from the genome file.

    Example:
    get_sequence_from_interval("hg38.fasta", "chr1:200-300")  # Returns "ATCGGCTA..."
    get_sequence_from_interval("hg38.fasta", "chr1:-1, 300")  # Raises ValueError
    get_sequence_from_interval("hg38.fasta", "chr1:200-1000000")  # Raises IndexError
    """

    try:
        # Split the interval into chromosome, start, and end
        chrom, start, end = parse_interval(interval)

        with pysam.FastaFile(reference_fasta) as ref_genome:
            # Fetch the sequence
            sequence = ref_genome.fetch(chrom, start-1, end)

        return sequence

    except KeyError as e:  # Catching errors related to chromosome name
        logging.exception("Error in fetching sequence", extra={'interval': interval, 'reference_fasta': reference_fasta})
        return None

    except IndexError as e:
        print("Index Error")
        logging.exception("Error in fetching sequence", extra={'interval': interval, 'reference_fasta': reference_fasta})
        raise

    except ValueError as e:  # Catching errors related to splitting and conversion to int
        logging.error(f"Error: Malformed input or could not convert data: {e}")
        return None

    except IOError as e:  # Catching errors related to file operations
        logging.error(f"Error: Failed to open or read the file: {e}")
        raise


def get_region_around_deletion(reference_fasta: str, interval: str, padding: int = 50) -> str:
    """
    Retrieve the genomic region around a given deletion site, extending by a specified flank length on either side.

    Parameters:
    chromosome (str): The chromosome where the deletion is located. E.g., "chr1".
    deletion_start (int): The start position of the deletion. E.g., 200.
    deletion_end (int): The end position of the deletion. E.g., 300.
    flank_length (int): Length of sequence to retrieve on either side of the deletion. E.g., 50.

    Returns:
    Tuple[str, int, int]: A tuple containing chromosome (str), start (int), and end (int) positions of the region.

    Raises:
    ValueError: If the start or end positions are negative, or if flank_length is negative.
    IndexError: If the region is out of bounds for the given chromosome.
    IOError: If there's a problem reading from the genome file.

    Example:
    get_region_around_deletion("chr1", 200, 300, 50)  # Returns DNA sequence in the interval chr1:150-350
    get_region_around_deletion("chr1", -1, 300, 50)  # Raises ValueError
    get_region_around_deletion("chr1", 200, 1000000, 50)  # Raises IndexError
    """

    chrom, start, end = parse_interval(interval)
    left = get_sequence_from_interval(reference_fasta, f"{chrom}:{start-padding}-{start-1}")
    right = get_sequence_from_interval(reference_fasta, f"{chrom}:{end+1}-{end+padding}")

    return left + right


def get_flanking_regions(reference_fasta: str, interval: str, padding: int = 50) -> Tuple[str, str]:
    """
    Retrieve the genomic regions around a given interval, extending by a specified flank length on either side.

    Parameters:
    reference_fasta (str): Path to the reference genome in FASTA format.
    interval (str): The interval to pad
    padding (int): Length of sequence to retrieve on either side of the interval. E.g., 50.

    Returns:
    Tuple[str, int, int]: A tuple containing chromosome (str), start (int), and end (int) positions of the region.

    Raises:
    ValueError: If the start or end positions are negative, or if flank_length is negative.
    IndexError: If the region is out of bounds for the given chromosome.
    IOError: If there's a problem reading from the genome file.

    Example:
    get_flanking_regions("chr1", 200, 300, 50)  # Returns DNA sequence in the interval chr1:150-350
    get_flanking_regions("chr1", -1, 300, 50)  # Raises ValueError
    get_flanking_regions("chr1", 200, 1000000, 50)  # Raises IndexError
    """

    chrom, start, end = parse_interval(interval)
    left = get_sequence_from_interval(reference_fasta, f"{chrom}:{start-padding}-{start-1}")
    right = get_sequence_from_interval(reference_fasta, f"{chrom}:{end+1}-{end+padding}")

    return left, right
