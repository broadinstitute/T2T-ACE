import logging

import mappy

from T2T_ACE.interval_parsing import create_interval
from T2T_ACE.genomic_queries import get_sequence_from_interval, get_flanking_regions, get_region_around_deletion
from T2T_ACE.interval_parsing import (parse_interval, find_next_interval, find_previous_interval,
                                      distance_between_intervals, interval_between_intervals, interval_within_interval,
                                      interval_size)
import T2T_ACE.alignment_utilities as au
import T2T_ACE.alignment_visualization_utilities as avu
import sys


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


def check_interval(interval: str, calling_reference_fasta: str,
                   called_ref_aligner: mappy.Aligner, truth_ref_aligner: mappy.Aligner, upper_bound=1e6):
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

    left_flank, right_flank = get_flanking_regions(calling_reference_fasta, interval, padding=5000)
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


def align_flanking_sequences(interval, calling_reference_fasta: str, called_ref_aligner, truth_ref_aligner, upper_bound=1000000, padding_length=500):
    _, start, stop = parse_interval(interval)

    if stop - start > upper_bound:
        raise ValueError("Interval too large")

    left_flank, right_flank = get_flanking_regions(calling_reference_fasta, interval, padding=padding_length)
    left_flank_hg38_hits = [extract_interval_from_hit(_) for _ in called_ref_aligner.map(left_flank) if 'alt' not in _.ctg]
    left_flank_hg002t2t_hits = [extract_interval_from_hit(_) for _ in truth_ref_aligner.map(left_flank)]
    right_flank_hg38_hits = [extract_interval_from_hit(_) for _ in called_ref_aligner.map(right_flank) if 'alt' not in _.ctg]
    right_flank_hg002t2t_hits = [extract_interval_from_hit(_) for _ in truth_ref_aligner.map(right_flank)]

    return {"left_ref_flank": left_flank_hg38_hits,
            "left_truth_flank": left_flank_hg002t2t_hits,
            "right_ref_flank": right_flank_hg38_hits,
            "right_truth_flank": right_flank_hg002t2t_hits}


def get_flanking_pairs(interval, calling_reference_fasta: str, called_ref_aligner, truth_ref_aligner, upper_bound=1000000, padding_length=500):
    flanking_alignments = align_flanking_sequences(interval, calling_reference_fasta,
                                                   called_ref_aligner, truth_ref_aligner, upper_bound, padding_length)

    stuff = {}
    try:
        stuff["ref_flank"] = []
        for flanking_alignment in flanking_alignments["left_ref_flank"]:
            next_interval = find_next_interval(flanking_alignment, flanking_alignments["right_ref_flank"])
            if next_interval is not None:
                # in_between = interval_between_intervals(flanking_alignment, next_interval)
                print(flanking_alignment, next_interval, distance_between_intervals(flanking_alignment, next_interval))
                stuff["ref_flank"].append((flanking_alignment, next_interval, distance_between_intervals(flanking_alignment, next_interval)))

        stuff["truth_flank"] = []
        for flanking_alignment in flanking_alignments["left_truth_flank"]:
            next_interval = find_next_interval(flanking_alignment, flanking_alignments["right_truth_flank"])
            if next_interval is not None:
                # in_between = interval_between_intervals(flanking_alignment, next_interval)
                print(flanking_alignment, next_interval, distance_between_intervals(flanking_alignment, next_interval))
                stuff["truth_flank"].append((flanking_alignment, next_interval, distance_between_intervals(flanking_alignment, next_interval)))


    except:
        raise ValueError("No flanking alignments found")

    return stuff

def stuff2list(interval, stuff:dict):

    ref_intervals = [interval, stuff['ref_flank'][0][0], stuff['ref_flank'][0][1]]
    truth_intervals = [stuff['truth_flank'][0][0], stuff['truth_flank'][0][1], stuff['truth_flank'][1][0], stuff['truth_flank'][1][1]]
    return ref_intervals, truth_intervals

def align_interval(interval, calling_reference_fasta: str, called_ref_aligner, truth_ref_aligner) -> list:
    interval_seq = get_sequence_from_interval(calling_reference_fasta, interval)
    interval_hg2_hits = [[extract_interval_from_hit(_), _.strand] for _ in truth_ref_aligner.map(interval_seq)]
    # Collect all the hg38 alignments that are not alt contigs
    interval_hg38_hits = [[extract_interval_from_hit(_), _.strand] for _ in called_ref_aligner.map(interval_seq) if 'alt' not in _.ctg]
    return interval_hg38_hits, interval_hg2_hits

def eval_del_in_dup(del_interval, dup_interval, calling_reference_fasta: str, called_ref_aligner, truth_ref_aligner, plot=False, plot_ratio=70, save_plot=False):
    # Check if the input DEL interval is within the DUP interval
    if interval_within_interval(del_interval, dup_interval):
        # First Align the DUP sequence to HG2 and hg38
        dup_alignments = align_interval(dup_interval, calling_reference_fasta, called_ref_aligner, truth_ref_aligner)
        del_alignments = align_interval(del_interval, calling_reference_fasta, called_ref_aligner, truth_ref_aligner)

        # Check the number of alignments for DUP and DEL in HG2 and hg38
        hg38_dup_count = len(dup_alignments[0])
        hg2_dup_count = len(dup_alignments[1])
        hg38_del_count = len(del_alignments[0])
        hg2_del_count = len(del_alignments[1])

        # Check if the DUP is a DUP (*2 bc HG2 is diploid)
        if hg2_dup_count > hg38_dup_count*2:
            print(f"{dup_interval} is a real DUP. It has {hg2_dup_count} copy(ies) in HG2 and {hg38_dup_count} copy(ies) in hg38\n")
        else:
            sys.exit(f"{dup_interval} is not a DUP. It has {hg2_dup_count} copy(ies) in HG2 and {hg38_dup_count} copy(ies) in hg38\n")

        # print the alignments
        print(f"DUP copies in hg38:{hg38_dup_count}")
        print(f"DUP copies in HG2:{hg2_dup_count}\n")
        print(f"DEL copies in hg38:{hg38_del_count}")
        print(f"DEL copies in HG2:{hg2_del_count}\n")

        # If DEL has copies in HG2, check if they occur in the DUP copy intervals in HG2
        if hg2_del_count != 0:
            for dup_alignment in dup_alignments:
                # print(alignment)
                for dup_aln_interval in dup_alignment:
                    if '_' in dup_aln_interval[0]:
                        # HG2-T2T alignments intervals of DUP
                        DUP_copy = dup_aln_interval[0]
                        for del_alignment in del_alignments:
                            for del_aln_interval in del_alignment:
                                if '_' in del_aln_interval[0]:
                                    # HG2-T2T alignments intervals of DEL
                                    DEL_copy = del_aln_interval[0]
                                    # If they are on the same chromosome, check if the DEL is within the DUP copy
                                    if DEL_copy.split(':')[0] == DUP_copy.split(':')[0]:
                                        if interval_within_interval(DEL_copy, DUP_copy):
                                            print('DEL in DUP copy:', DEL_copy, DUP_copy)
        # Collect Flanking regions
        dup_chrom, dup_pos, dup_end = parse_interval(dup_interval)
        del_chrom, del_pos, del_end = parse_interval(del_interval)

        left_flanking_interval = create_interval(dup_chrom, dup_pos, del_pos - 1)
        right_flanking_interval = create_interval(dup_chrom, del_end, dup_end)

        # Align the flanking regions to HG2 and hg38
        left_flanking_hg2_alignments = align_interval(left_flanking_interval, calling_reference_fasta, called_ref_aligner, truth_ref_aligner)[1]
        right_flanking_hg2_alignments = align_interval(right_flanking_interval, calling_reference_fasta, called_ref_aligner, truth_ref_aligner)[1]
        # Grab intervals of right flanking region alignments
        right_intervals = [i[0] for i in right_flanking_hg2_alignments]
        left_intervals = [i[0] for i in left_flanking_hg2_alignments]

        print(f'\nLeft flanking region: {left_flanking_interval}; HG2 Copies: {len(left_intervals)}')
        print(f'DEL length: {interval_size(del_interval)} bp')
        print(f'Right flanking region: {right_flanking_interval}; HG2 Copies: {len(right_intervals)}')
        print('\n')

        # Check the distance between the flanking regions' alignments
        # If the distance is less than half of the DEL interval, print the flanking regions' alignments
        for left_interval, strand in left_flanking_hg2_alignments:
            # If the alignment is on the same strand as the DUP, check if the next interval of left flanking region's alignment
            if strand == 1:
                distance_between_flankings = distance_between_intervals(left_interval, find_next_interval(left_interval, right_intervals))
                # print('+', left_interval, ip.find_next_interval(left_interval,right_intervals), distance_between_flankings)
                if distance_between_flankings < int(interval_size(del_interval)) * 0.5:
                    print('*', left_interval, find_next_interval(left_interval, right_intervals),
                          distance_between_flankings)
            # If the alignment is on the opposite strand as the DUP, check if the previous interval of left flanking region's alignment
            elif strand == -1:
                distance_between_flankings = distance_between_intervals(left_interval, find_previous_interval(left_interval, right_intervals))
                # print('-', left_interval, ip.find_previous_interval(left_interval,right_intervals), distance_between_flankings)
                if distance_between_flankings < int(interval_size(del_interval)) * 0.5:
                    print('*', left_interval, find_previous_interval(left_interval, right_intervals),
                          distance_between_flankings)

        # If plot is True, plot the alignments of the flanking regions
        if plot:
            hg38_flanking_intervals = [left_flanking_interval, right_flanking_interval]
            hg2_flanking_alignment_intervals = [i[0] for i in right_flanking_hg2_alignments] + [i[0] for i in left_flanking_hg2_alignments]
            avu.PlotIntervals(hg38_flanking_intervals, hg2_flanking_alignment_intervals).plot_intervals_comparison(flanking=False, ratio=plot_ratio, save=save_plot)
    else:
        sys.exit(f"DEL interval {del_interval} not within DUP interval {dup_interval}")

