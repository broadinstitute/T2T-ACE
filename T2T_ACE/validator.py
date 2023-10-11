import logging

import mappy

from interval_parsing import create_interval
from genomic_queries import get_sequence_from_interval, get_flanking_regions, get_region_around_deletion
from interval_parsing import parse_interval
import alignment_utilities as au



def log_error(level, msg, *args):
    logging.log(level, f"Error: {msg}", *args)


def evaluate_region(calling_reference: str, interval: str, calling_aligner, t2t_aligner):
    raw_sequence = get_sequence_from_interval(calling_reference, interval)
    left_flank, right_flank = get_flanking_regions(calling_reference, interval, padding=500)

    hg38_hits = [_ for _ in calling_aligner.map(raw_sequence)]
    hg002t2t_hits = [_ for _ in t2t_aligner.map(raw_sequence)]

    left_flank_hg38_hits = [_ for _ in calling_aligner.map(left_flank)]
    left_flank_hg002t2t_hits = [_ for _ in t2t_aligner.map(left_flank)]

    right_flank_hg38_hits = [_ for _ in calling_aligner.map(right_flank)]
    right_flank_hg002t2t_hits = [_ for _ in t2t_aligner.map(right_flank)]

    return [hg38_hits, hg002t2t_hits,
            left_flank_hg38_hits, left_flank_hg002t2t_hits,
            right_flank_hg38_hits, right_flank_hg002t2t_hits]


def evaluate_deletion(calling_reference: str, truth_reference: str, interval: str, calling_aligner, t2t_aligner,
                      make_plots: bool = False) -> list:

    # If HET deletion, we expect to see the deleted_sequence once in the hg002 reference
    # If HOM deletion, we expect to not see the deleted sequence in the hg002 reference
    raw_sequence = get_sequence_from_interval(calling_reference, interval)
    left_flank, right_flank = get_flanking_regions(calling_reference, interval, padding=500)

    hg38_hits = [_ for _ in calling_aligner.map(raw_sequence)]
    hg002t2t_hits = [_ for _ in t2t_aligner.map(raw_sequence)]

    left_flank_hg38_hits = [_ for _ in calling_aligner.map(left_flank)]
    left_flank_hg002t2t_hits = [_ for _ in t2t_aligner.map(left_flank)]

    right_flank_hg38_hits = [_ for _ in calling_aligner.map(right_flank)]
    right_flank_hg002t2t_hits = [_ for _ in t2t_aligner.map(right_flank)]

    au.print_hits("Raw Sequence hg38", len(raw_sequence), hg38_hits)
    au.print_hits("Raw Sequence hg002", len(raw_sequence), hg002t2t_hits)

    au.print_hits("Left Flank Sequence hg38", len(left_flank), left_flank_hg38_hits)
    au.print_hits("Left Flank Sequence hg002", len(left_flank), left_flank_hg002t2t_hits)

    au.print_hits("Right Flank Sequence hg38", len(right_flank), right_flank_hg38_hits)
    au.print_hits("Right Flank Sequence hg002", len(right_flank), right_flank_hg002t2t_hits)

    if make_plots is True:
        #for hit in hg38_hits:
        #    avu.dot_plot(seq1=raw_sequence,
        #                 seq2=augmented_sequence,
        #                 seq1_name=interval,
        #                 seq2_name=hit.ctg)

        for hit in left_flank_hg002t2t_hits:
            print("Left flank: ", create_interval(hit.ctg, hit.r_st+1, hit.r_en))
            #avu.dot_plot(seq1=augmented_sequence,
            #             seq2=get_sequence_from_interval(truth_reference, create_interval(hit.ctg, hit.r_st+1, hit.r_en)),
            #             seq1_name=interval,
            #             seq2_name=create_interval(hit.ctg, hit.r_st+1, hit.r_en))

        for hit in right_flank_hg002t2t_hits:
            print("Right flank: ", create_interval(hit.ctg, hit.r_st+1, hit.r_en))
            #avu.dot_plot(seq1=augmented_sequence,
            #             seq2=get_sequence_from_interval(truth_reference, create_interval(hit.ctg, hit.r_st+1, hit.r_en)),
            #             seq1_name=interval,
            #             seq2_name=create_interval(hit.ctg, hit.r_st+1, hit.r_en))

    return [hg38_hits, hg002t2t_hits]


def check_interval(interval, calling_reference_fasta: str, called_ref_aligner, truth_ref_aligner, upper_bound=1000000):
    _, start, stop = parse_interval(interval)

    raw_sequence = get_sequence_from_interval(calling_reference_fasta, interval)
    if stop - start > upper_bound:
        return "Interval too large", [0, 0, 0, 0]

    joined_sequence = get_region_around_deletion(calling_reference_fasta, interval)

    hg38_raw_alignments = [_ for _ in called_ref_aligner.map(raw_sequence)]
    hg002_alignments = [_ for _ in truth_ref_aligner.map(raw_sequence)]

    hg38_count = [au.is_good_hit(_, len(raw_sequence)) for _ in called_ref_aligner.map(raw_sequence)].count(True)
    hg002_count = [au.is_good_hit(_, len(raw_sequence)) for _ in truth_ref_aligner.map(raw_sequence)].count(True)

    hg38_joined_alignments = [_ for _ in called_ref_aligner.map(joined_sequence)]
    hg002_joined_alignments = [_ for _ in truth_ref_aligner.map(joined_sequence)]

    hg38_joined_count = [au.is_good_hit(_, len(joined_sequence))
                         for _ in called_ref_aligner.map(joined_sequence)].count(True)
    hg002_joined_count = [au.is_good_hit(_, len(joined_sequence))
                          for _ in truth_ref_aligner.map(joined_sequence)].count(True)

    left_flank, right_flank = get_flanking_regions(calling_reference_fasta, interval, padding=500)
    left_flank_hg38_hits = [_ for _ in called_ref_aligner.map(left_flank)]
    left_flank_hg002t2t_hits = [_ for _ in truth_ref_aligner.map(left_flank)]

    right_flank_hg38_hits = [_ for _ in called_ref_aligner.map(right_flank)]
    right_flank_hg002t2t_hits = [_ for _ in truth_ref_aligner.map(right_flank)]

    result = [hg38_count, hg002_count, hg38_joined_count, hg002_joined_count]
    alignment_results = [hg38_raw_alignments, hg002_alignments, hg38_joined_alignments, hg002_joined_alignments]

    if hg38_count == 0:
        return "Interval likely has Ns", result, alignment_results
    elif (hg38_count == 1 and
          hg002_count == 1 and
          hg38_joined_count == 0 and
          hg002_joined_count == 1):
        return "Het Deletion", result, alignment_results
    elif result == [1, 0, 0, 2]:
        return "Hom Deletion", result, alignment_results
    elif (result[1] > result[0]) and (result[2] == 0) and (result[3] == 0):
        return "Duplication Copy number: " + str(result[1] / result[0]), result, alignment_results


def extract_interval_from_hit(hit: mappy.Alignment):
    return hit.ctg + ":" + str(hit.r_st) + "-" + str(hit.r_en)


def align_flanking_sequences(interval, calling_reference_fasta: str, called_ref_aligner, truth_ref_aligner, upper_bound=1000000):
    _, start, stop = parse_interval(interval)

    if stop - start > upper_bound:
        raise ValueError("Interval too large")

    left_flank, right_flank = get_flanking_regions(calling_reference_fasta, interval, padding=500)
    left_flank_hg38_hits = [extract_interval_from_hit(_) for _ in called_ref_aligner.map(left_flank)]
    left_flank_hg002t2t_hits = [extract_interval_from_hit(_) for _ in truth_ref_aligner.map(left_flank)]
    right_flank_hg38_hits = [extract_interval_from_hit(_) for _ in called_ref_aligner.map(right_flank)]
    right_flank_hg002t2t_hits = [extract_interval_from_hit(_) for _ in truth_ref_aligner.map(right_flank)]

    return {"left_ref_flank": left_flank_hg38_hits,
            "left_truth_flank": left_flank_hg002t2t_hits,
            "right_ref_flank": right_flank_hg38_hits,
            "right_truth_flank": right_flank_hg002t2t_hits}
