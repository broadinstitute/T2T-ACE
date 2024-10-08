import logging
import numpy as np
import pandas as pd
import mappy
from . import alignment_utilities as au
from . import alignment_visualization_utilities as avu
import sys
import os
from.interval_parsing import create_interval
from.genomic_queries import get_sequence_from_interval, get_flanking_regions, get_region_around_deletion
from.interval_parsing import (parse_interval, find_next_interval, find_previous_interval,
                                      distance_between_intervals, interval_between_intervals, interval_within_interval,
                                      interval_size, get_reversed_sequence, interval_overlapping_percentage)

# Load the chromosome sizes
current_dir = os.path.dirname(os.path.abspath(__file__))
hg38_chrom_size_filepath = os.path.join(current_dir, 'resources/hg38_chrom_size.txt')
chrom_size_df = pd.read_csv(hg38_chrom_size_filepath, sep='\t', header=None, names=['chr', 'size'])
chrom_size_dict = {row['chr']: row['size'] for index, row in chrom_size_df.iterrows()}

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
    interval_length = interval_size(interval)
    interval_seq = get_sequence_from_interval(calling_reference_fasta, interval)
    # Collect all alignments that are at least 95% of the interval length
    interval_hg2_hits = [[extract_interval_from_hit(_), _.strand, _.q_st, _.q_en] for _ in truth_ref_aligner.map(interval_seq) if (_.q_en - _.q_st + 1) / interval_length > 0.5]
    # Collect all the hg38 alignments that are not alt contigs
    interval_hg38_hits = [[extract_interval_from_hit(_), _.strand, _.q_st, _.q_en] for _ in called_ref_aligner.map(interval_seq) if 'alt' not in _.ctg and 'random' not in _.ctg and 'chrUn' not in _.ctg and (_.q_en - _.q_st + 1) / interval_length > 0.5]
    return interval_hg38_hits, interval_hg2_hits

def align_interval_no_restriction(interval, calling_reference_fasta: str, called_ref_aligner, truth_ref_aligner) -> list:
    interval_length = interval_size(interval)
    interval_seq = get_sequence_from_interval(calling_reference_fasta, interval)
    # Collect all alignments that are at least 95% of the interval length
    interval_hg2_hits = [[extract_interval_from_hit(_), _.strand, _.q_st, _.q_en] for _ in truth_ref_aligner.map(interval_seq)]
    # Collect all the hg38 alignments that are not alt contigs
    interval_hg38_hits = [[extract_interval_from_hit(_), _.strand, _.q_st, _.q_en] for _ in called_ref_aligner.map(interval_seq) if 'alt' not in _.ctg and 'random' not in _.ctg and 'chrUn' not in _.ctg]
    return interval_hg38_hits, interval_hg2_hits

# Report the alignment of the DUP interval to the two references (hg38 and HG2)
def reportAlignment(dup_interval, calling_reference_fasta, called_ref_aligner, truth_ref_aligner, print_alignments=False, alignment_filter=True):
        # First Align the DUP sequence to HG2 and hg38
        if alignment_filter:
            dup_alignments = align_interval(dup_interval, calling_reference_fasta, called_ref_aligner, truth_ref_aligner)
        else:
            dup_alignments = align_interval_no_restriction(dup_interval, calling_reference_fasta, called_ref_aligner, truth_ref_aligner)

        # Check the number of alignments for DUP and DEL in HG2 and hg38
        hg38_dup_count = len(dup_alignments[0])
        hg2_dup_count = len(dup_alignments[1])
        hg2_mat_copies = [[interval, strand, q_start, q_end] for interval, strand, q_start, q_end in dup_alignments[1] if 'MAT' in interval]
        hg2_pat_copies = [[interval, strand, q_start, q_end] for interval, strand, q_start, q_end in dup_alignments[1] if 'PAT' in interval]

        if print_alignments:
            print(f"input dup interval: {dup_interval}")
            print(f"hg38 dup count: {hg38_dup_count}")
            for interval, strand, q_start, q_end in dup_alignments[0]:
                print(f"interval: {interval}\tstrand: {strand}, start: {q_start}, end: {q_end}")
            print(f"hg2 dup count: {hg2_dup_count}")
            for interval, strand, q_start, q_end in hg2_mat_copies:
                print(f"interval: {interval}\tstrand: {strand}, start: {q_start}, end: {q_end}")
            for interval, strand, q_start, q_end in hg2_pat_copies:
                print(f"interval: {interval}\tstrand: {strand}, start: {q_start}, end: {q_end}")

        return(hg38_dup_count, len(hg2_mat_copies), len(hg2_pat_copies), hg2_dup_count)

def classifyDupInterval(dup_interval, calling_reference_fasta, called_ref_aligner, truth_ref_aligner, alignment_filter=True):
    if alignment_filter:
        dup_alignments = align_interval(dup_interval, calling_reference_fasta, called_ref_aligner, truth_ref_aligner)
    else:
        dup_alignments = align_interval_no_restriction(dup_interval, calling_reference_fasta, called_ref_aligner,
                                                       truth_ref_aligner)

    # Check the number of alignments for DUP and DEL in HG2 and hg38
    hg38_dup_count = len(dup_alignments[0])
    hg2_dup_count = len(dup_alignments[1])
    hg2_mat_copies = [[interval, strand, q_start, q_end] for interval, strand, q_start, q_end in dup_alignments[1] if
                      'MAT' in interval]
    hg2_pat_copies = [[interval, strand, q_start, q_end] for interval, strand, q_start, q_end in dup_alignments[1] if
                      'PAT' in interval]

    hg2_pat_count = len(hg2_pat_copies)
    hg2_mat_count = len(hg2_mat_copies)
    # Gather the chromosomes of the DUP copies in HG2
    hg2_mat_interval_chrs = [interval.split(':')[0].split("_")[0] for interval, strand, q_start, q_end in hg2_mat_copies]
    hg2_pat_interval_chrs = [interval.split(':')[0].split("_")[0] for interval, strand, q_start, q_end in
                             hg2_pat_copies]

    chr, pos, end = parse_interval(dup_interval)

    # Classify the DUP interval
    # HG2 doesn't have the copy of the DUP interval in the chromosome that was called in hg38
    if chr not in hg2_mat_interval_chrs and chr not in hg2_pat_interval_chrs:
        major_classification = "Reference Error" # Might update this name to "Translocation"
        sub_classification = "hg38 Reference Error"
    elif hg38_dup_count == 0:
        major_classification = "Minimap2 Error"
        sub_classification = "Minimap2 Error"
    elif hg2_mat_count > hg38_dup_count and hg2_pat_count > hg38_dup_count:
        major_classification = "Duplication"
        sub_classification = "Homozygous Duplication"
    elif hg2_mat_count > hg38_dup_count and hg2_pat_count <= hg38_dup_count:
        major_classification = "Duplication"
        sub_classification = "Maternal Heterozygous Duplication"
    # Taking account of Sex chromosomes (This is a temporary fix for HG2/Male sample)
    elif chr == 'chrX' and hg2_mat_count > hg38_dup_count and hg2_pat_count == 0:
        major_classification = "Duplication"
        sub_classification = "Maternal Heterozygous Duplication"
    elif hg2_pat_count > hg38_dup_count and hg2_mat_count <= hg38_dup_count:
        major_classification = "Duplication"
        sub_classification = "Paternal Heterozygous Duplication"
    elif chr == 'chrY' and hg2_pat_count > hg38_dup_count and hg2_mat_count == 0:
        major_classification = "Duplication"
        sub_classification = "Paternal Heterozygous Duplication"
    elif hg2_mat_count == hg38_dup_count and hg2_pat_count == hg38_dup_count:
        major_classification = "Copy Neutral"
        sub_classification = "Biallelic Copy Neutral"
    elif hg2_mat_count == hg38_dup_count and hg2_pat_count ==0:
        major_classification = "Copy Neutral"
        sub_classification = "Maternal Copy Neutral"
    elif hg2_pat_count == hg38_dup_count and hg2_mat_count == 0:
        major_classification = "Copy Neutral"
        sub_classification = "Paternal Copy Neutral"
    elif hg38_dup_count == 1 and hg2_mat_count == 0 and hg2_pat_count == 0:
        major_classification = "Reference Error"
        sub_classification = "hg38 Reference Error"
    else:
        major_classification = "Unknown"
        sub_classification = "Unknown"
    return major_classification, sub_classification

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

        # print(right_flanking_hg2_alignments)
        # print(left_flanking_hg2_alignments)

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
                else:
                    print(left_interval, find_next_interval(left_interval, right_intervals),
                          distance_between_flankings)
            # If the alignment is on the opposite strand as the DUP, check if the previous interval of left flanking region's alignment
            elif strand == -1:
                distance_between_flankings = distance_between_intervals(left_interval, find_previous_interval(left_interval, right_intervals))
                # print('-', left_interval, ip.find_previous_interval(left_interval,right_intervals), distance_between_flankings)
                if distance_between_flankings < int(interval_size(del_interval)) * 0.5:
                    print('*', left_interval, find_previous_interval(left_interval, right_intervals),
                          distance_between_flankings)
                else:
                    print(left_interval, find_previous_interval(left_interval, right_intervals),
                          distance_between_flankings)

        # If plot is True, plot the alignments of the flanking regions
        if plot:
            hg38_flanking_intervals = [left_flanking_interval, right_flanking_interval]
            hg2_flanking_alignment_intervals = [i[0] for i in right_flanking_hg2_alignments] + [i[0] for i in left_flanking_hg2_alignments]
            avu.PlotIntervals(hg38_flanking_intervals, hg2_flanking_alignment_intervals).plot_intervals_comparison(flanking=False, ratio=plot_ratio, save=save_plot)
    else:
        sys.exit(f"DEL interval {del_interval} not within DUP interval {dup_interval}")

def collect_del_flankings(del_interval, calling_reference_fasta: str, called_ref_aligner, truth_ref_aligner, flanking_size=None):
    del_flankings_sum_dict = {}
    del_chrom, del_pos, del_end = parse_interval(del_interval)
    del_interval_size = interval_size(del_interval)
    # Add interval size check
    if del_interval_size > 1000000:
        print("DEL interval is too large")
        return None

    # Align the DEL interval to HG2 and hg38
    # DEL interval sequence should align to HG2 once if it's a HET DEL and zero times if it's a HOM DEL
    # If it aligns to HG2 more than once, it's likely to be a False DEL
    del_seq_alignments = align_interval(del_interval, calling_reference_fasta, called_ref_aligner, truth_ref_aligner)
    del_seq_hg38_alignment_intervals = [interval for interval, strand, q_start, q_end in
                                              del_seq_alignments[0] if f"{del_chrom}:" in interval]
    del_seq_hg2_alignment_intervals = [interval for interval, strand, q_start, q_end in del_seq_alignments[1] if f"{del_chrom}_" in interval]
    del_seq_hg38_hits = len(del_seq_hg38_alignment_intervals)
    del_seq_hg2_hits = len(del_seq_hg2_alignment_intervals)
    print('DEL interval:', del_interval, del_interval_size)
    print('DEL interval Alignments in hg38:', del_seq_hg38_hits)
    print('DEL interval Alignments in HG2:', del_seq_hg2_hits)
    """
    Flanking logic: if DEL is less than 1000 bp, flanking regions are 1000 bp on both sides. If 1000bp is
    not enough to have alignments in HG2, then increase the flanking region size by 1000 bp until there are
    alignments in HG2. If DEL is more than 1000 bp, flanking regions are 10% of the DEL interval on both sides.
    flanking size can also be specified by the user.
    """
    if flanking_size:
        flanking_size = flanking_size
    else:
        if del_interval_size < 6500:
            flanking_size = 2000
        else:
            flanking_size = int(del_interval_size * 0.3)
    # initiate of flanking intervals
    hg38_left_flanking_interval = None
    hg38_right_flanking_interval = None
    left_flanking_hg38_alignment_intervals = []
    right_flanking_hg38_alignment_intervals = []
    left_flanking_hg2_alignment_intervals = []
    right_flanking_hg2_alignment_intervals = []
    left_flanking_interval_aligned = []

    # Use while loop to increase the flanking size until there are at least two alignments in HG2
    # Assuming the alignment is one MAT nad one PAT
    # If one flanking sequence is heterozygous which is unlikely, at least another flanking seq should
    #  have two alignments in HG2 on one chromosome
    # Taking account of sex chrom
    if del_chrom == 'chrX' or del_chrom == 'chrY':
        copy_threshold = 1
    else:
        copy_threshold = 2
    while len(left_flanking_hg2_alignment_intervals) < copy_threshold or len(right_flanking_hg2_alignment_intervals) < copy_threshold:
        print(f"flanking_size: {flanking_size}, copy_threshold: {copy_threshold}")
        left_flanking_pos = del_pos - flanking_size
        left_flanking_end = del_pos
        right_flanking_pos = del_end
        right_flanking_end = del_end + flanking_size
        if right_flanking_end >= chrom_size_dict[del_chrom]:
            print("Extension exceeds the length of the chromosome")
            break
        if left_flanking_pos < 0:
            print("Extension exceeds the length of the chromosome")
            break
        hg38_left_flanking_interval = create_interval(del_chrom, left_flanking_pos, left_flanking_end)
        hg38_right_flanking_interval = create_interval(del_chrom, right_flanking_pos, right_flanking_end)

        # Align the flanking intervals to the reference genome and HG2
        left_flanking_interval_aligned = align_interval(hg38_left_flanking_interval, calling_reference_fasta, called_ref_aligner,
                                                          truth_ref_aligner)
        right_flanking_interval_aligned = align_interval(hg38_right_flanking_interval, calling_reference_fasta, called_ref_aligner,
                                                          truth_ref_aligner)
        # Gather the intervals of flanking region alignments in reference genome (hg38)
        # The alignments in hg38 and HG2 should be on the same chromosome as the DEL interval
        # If more than 1, it means this DEL interval might exist in a DUP region
        left_flanking_hg38_alignment_intervals = [interval for interval, strand, q_start, q_end in
                                                  left_flanking_interval_aligned[0] if f"{del_chrom}:" in interval]
        right_flanking_hg38_alignment_intervals = [interval for interval, strand, q_start, q_end in
                                                   right_flanking_interval_aligned[0] if f"{del_chrom}:" in interval]
        # Gather the intervals of flanking region alignments in HG2
        left_flanking_hg2_alignment_intervals = [interval for interval, strand, q_start, q_end in left_flanking_interval_aligned[1] if f"{del_chrom}_" in interval]
        right_flanking_hg2_alignment_intervals = [interval for interval, strand, q_start, q_end in right_flanking_interval_aligned[1] if f"{del_chrom}_" in interval]

        # Print out the flanking intervals, their sizes, and the number of alignments in HG002
        print(
            f"left flanking interval: {hg38_left_flanking_interval}, {interval_size(hg38_left_flanking_interval)}, {len(left_flanking_hg38_alignment_intervals)}, {len(left_flanking_hg2_alignment_intervals)}")
        print(
            f"right flanking interval: {hg38_right_flanking_interval}, {interval_size(hg38_right_flanking_interval)}, {len(right_flanking_hg38_alignment_intervals)}, {len(right_flanking_hg2_alignment_intervals)}")

        # Increase the flanking size by 1000bp if there is less than two alignments in HG002
        # 1000 is arbitrary. But it seems to work well for most of the cases.
        if len(left_flanking_hg2_alignment_intervals) < copy_threshold or len(right_flanking_hg2_alignment_intervals) < copy_threshold:
            # If the flanking size is larger than the DEL interval, break the loop
            if flanking_size > del_interval_size and flanking_size > 10000:
                print(f"Flanking size ({flanking_size}bp) exceeds the size of the DEL interval ({del_interval_size}bp)")
                break
            elif del_interval_size > 100000:
                 flanking_size = flanking_size + int(np.round(0.1 * del_interval_size,0))
            else:
                flanking_size = flanking_size + 1000


    # After flanking size is determined, add general information to the dictionary
    del_flankings_sum_dict['del_interval'] = del_interval
    del_flankings_sum_dict['del_interval_size'] = del_interval_size
    del_flankings_sum_dict['flanking_size'] = flanking_size
    del_flankings_sum_dict['left_flanking_interval'] = hg38_left_flanking_interval
    del_flankings_sum_dict['right_flanking_interval'] = hg38_right_flanking_interval
    del_flankings_sum_dict['left_flanking_hg38_hits'] = len(left_flanking_hg38_alignment_intervals)
    del_flankings_sum_dict['right_flanking_hg38_hits'] = len(right_flanking_hg38_alignment_intervals)
    del_flankings_sum_dict['left_flanking_hg2_hits'] = len(left_flanking_hg2_alignment_intervals)
    del_flankings_sum_dict['right_flanking_hg2_hits'] = len(right_flanking_hg2_alignment_intervals)

    # Calculate the distance between the flanking regions' alignments
    # If the distance is less than half of the DEL interval means this is evidence for DEL event
    # Otherwise, no evidence for DEL event
    distance_between_flankings_list = []
    flanking_connection_strand_list = []
    left_aligned_interval_list = []
    right_aligned_interval_list = []
    interval_between_matching_flankings_list = []
    for left_flanking_interval, strand, q_start, q_end in left_flanking_interval_aligned[1]:
        if left_flanking_interval in left_flanking_hg2_alignment_intervals:
            if strand == 1:
                if find_next_interval(left_flanking_interval, right_flanking_hg2_alignment_intervals):
                    matching_right_flanking_interval = find_next_interval(left_flanking_interval,
                                                                          right_flanking_hg2_alignment_intervals)
                    distance_between_flankings = distance_between_intervals(left_flanking_interval,
                                                                              matching_right_flanking_interval)
                    if distance_between_flankings > 0:
                        interval_between_matching_flankings = interval_between_intervals(left_flanking_interval, matching_right_flanking_interval)
                    else:
                        interval_between_matching_flankings = None
                    distance_between_flankings_list.append(distance_between_flankings)
                    flanking_connection_strand_list.append("POS")
                    left_aligned_interval_list.append(left_flanking_interval)
                    right_aligned_interval_list.append(matching_right_flanking_interval)
                    interval_between_matching_flankings_list.append(interval_between_matching_flankings)
                    if distance_between_flankings < int(del_interval_size) * 0.5:
                        print('********** Potential DEL **********')
                    else:
                        print("----------- No DEL Evidence -----------")
                    print(
                        f"{left_flanking_interval} ({interval_size(left_flanking_interval)}bp), {matching_right_flanking_interval} ({interval_size(matching_right_flanking_interval)}bp), {distance_between_flankings}")
                elif find_previous_interval(left_flanking_interval, right_flanking_hg2_alignment_intervals):
                    matching_right_flanking_interval = find_previous_interval(left_flanking_interval,
                                                                              right_flanking_hg2_alignment_intervals)
                    distance_between_flankings = distance_between_intervals(left_flanking_interval,
                                                                              matching_right_flanking_interval)
                    if distance_between_flankings > 0:
                        interval_between_matching_flankings = interval_between_intervals(left_flanking_interval, matching_right_flanking_interval)
                    else:
                        interval_between_matching_flankings = None
                    distance_between_flankings_list.append(distance_between_flankings)
                    flanking_connection_strand_list.append("NEG")
                    left_aligned_interval_list.append(left_flanking_interval)
                    right_aligned_interval_list.append(matching_right_flanking_interval)
                    interval_between_matching_flankings_list.append(interval_between_matching_flankings)
                    if distance_between_flankings < int(del_interval_size) * 0.5:
                        print('********** Potential DEL **********')
                    else:
                        print("----------- No DEL Evidence -----------")
                    print(
                        f"{left_flanking_interval} ({interval_size(left_flanking_interval)}bp), {matching_right_flanking_interval} ({interval_size(matching_right_flanking_interval)}bp), {distance_between_flankings}")
                else:
                    distance_between_flankings_list.append(None)
                    flanking_connection_strand_list.append(None)
                    left_aligned_interval_list.append(None)
                    right_aligned_interval_list.append(None)
                    interval_between_matching_flankings_list.append(None)
                    print(f"No matching right flanking interval for {left_flanking_interval}")
            elif strand == -1:
                if find_previous_interval(left_flanking_interval, right_flanking_hg2_alignment_intervals):
                    matching_right_flanking_interval = find_previous_interval(left_flanking_interval,
                                                                                right_flanking_hg2_alignment_intervals)
                    distance_between_flankings = distance_between_intervals(left_flanking_interval,
                                                                              matching_right_flanking_interval)
                    if distance_between_flankings > 0:
                        interval_between_matching_flankings = interval_between_intervals(left_flanking_interval, matching_right_flanking_interval)
                    else:
                        interval_between_matching_flankings = None
                    distance_between_flankings_list.append(distance_between_flankings)
                    flanking_connection_strand_list.append("POS")
                    left_aligned_interval_list.append(left_flanking_interval)
                    right_aligned_interval_list.append(matching_right_flanking_interval)
                    interval_between_matching_flankings_list.append(interval_between_matching_flankings)
                    if distance_between_flankings < int(del_interval_size) * 0.5:
                        print('********** Potential DEL**********')
                    else:
                        print("----------- No DEL Evidence -----------")
                    print(
                        f"{left_flanking_interval} ({interval_size(left_flanking_interval)}bp), {matching_right_flanking_interval} ({interval_size(matching_right_flanking_interval)}bp), {distance_between_flankings}")
                elif find_next_interval(left_flanking_interval, right_flanking_hg2_alignment_intervals):
                    matching_right_flanking_interval = find_next_interval(left_flanking_interval,
                                                                              right_flanking_hg2_alignment_intervals)
                    distance_between_flankings = distance_between_intervals(left_flanking_interval,
                                                                              matching_right_flanking_interval)
                    if distance_between_flankings > 0:
                        interval_between_matching_flankings = interval_between_intervals(left_flanking_interval, matching_right_flanking_interval)
                    else:
                        interval_between_matching_flankings = None
                    distance_between_flankings_list.append(distance_between_flankings)
                    flanking_connection_strand_list.append("NEG")
                    left_aligned_interval_list.append(left_flanking_interval)
                    right_aligned_interval_list.append(matching_right_flanking_interval)
                    interval_between_matching_flankings_list.append(interval_between_matching_flankings)
                    if distance_between_flankings < int(del_interval_size) * 0.5:
                        print('********** Potential DEL **********')
                    else:
                        print("----------- No DEL Evidence -----------")
                    print(
                        f"{left_flanking_interval} ({interval_size(left_flanking_interval)}bp), {matching_right_flanking_interval} ({interval_size(matching_right_flanking_interval)}bp), {distance_between_flankings}")
                else:
                    distance_between_flankings_list.append(None)
                    flanking_connection_strand_list.append(None)
                    left_aligned_interval_list.append(None)
                    right_aligned_interval_list.append(None)
                    interval_between_matching_flankings_list.append(None)
                    print(f"No matching right flanking interval for {left_flanking_interval}")

    distance_between_flankings_list = [i for i in distance_between_flankings_list if i != None]
    flanking_connection_strand_list = [i for i in flanking_connection_strand_list if i != None]
    # Add the distance between the flanking regions' alignments and strand to the dictionary
    del_flankings_sum_dict['distance_between_flankings'] = distance_between_flankings_list
    del_flankings_sum_dict['flanking_connection_strand'] = flanking_connection_strand_list
    del_flankings_sum_dict['hg38_plotting_flanking_intervals'] = [del_interval, hg38_left_flanking_interval, hg38_right_flanking_interval]
    del_flankings_sum_dict['hg2_plotting_flanking_intervals'] = left_aligned_interval_list + right_aligned_interval_list
    print("PRINTING DISTANCE BETWEEN FLANKINGS LIST")
    print(distance_between_flankings_list)

    # Use the obtained details about the flanking regions to classify the DEL interval
    if len(left_flanking_hg38_alignment_intervals) > 1 and len(right_flanking_hg38_alignment_intervals) > 1 and len(distance_between_flankings_list) > 2:
        if len(distance_between_flankings_list)>0 and min(distance_between_flankings_list) <= del_interval_size * 0.5:
            major_classification = 'DEL'
            minor_classification = 'DEL in DUP'
        else:
            major_classification = 'False DEL'
            minor_classification = 'False DEL'
            interval_between_matching_flankings_list = [interval for interval in interval_between_matching_flankings_list if interval != None]
            # print(interval_between_matching_flankings_list)
            print('----------- Checking DEL Sequence Alignment in HG2-----------')
            for distance_interval in interval_between_matching_flankings_list:
                    distance_interval_chr, distance_interval_start, distance_interval_end = parse_interval(
                        distance_interval)
                    for del_hg2_seq in del_seq_hg2_alignment_intervals:
                        del_hg2_chr, del_hg2_start, del_hg2_end = parse_interval(del_hg2_seq)
                        if distance_interval_chr == del_hg2_chr:
                            if min(interval_overlapping_percentage(distance_interval, del_hg2_seq)) > 0:
                                overlap_to_distance_between_flankings, overlap_to_del_seq = interval_overlapping_percentage(distance_interval, del_hg2_seq)
                                print(f"DEL HG2 alignment:{del_hg2_seq}; Interval Between HG2 Aligned Flankings:{distance_interval}\nThe overlap between the DEL sequence and the interval between the HG2 aligned flankings is {overlap_to_del_seq} of the DEL sequence and {overlap_to_distance_between_flankings} of the interval between the HG2 aligned flankings")
    else:
        del_flanking_dist_list = [distance for distance in distance_between_flankings_list if (distance <= del_interval_size * 0.5)]
        print("PRINTING DEL FLANKING DISTANCE LIST")
        print(del_flanking_dist_list)
        if len(distance_between_flankings_list) == len(del_flanking_dist_list) and len(del_flanking_dist_list) > 0:
            if del_chrom != 'chrX' or del_chrom != 'chrY':
                if len(del_flanking_dist_list) == 2:
                    major_classification = 'DEL'
                    minor_classification = 'Homozygous DEL'
                elif len(del_flanking_dist_list) == 1:
                    major_classification = 'DEL'
                    minor_classification = 'Heterozygous DEL'
                else:
                    major_classification = 'DEL'
                    minor_classification = 'Unclassified DEL'
            else:
                major_classification = 'DEL'
                minor_classification = 'Heterozygous DEL'
        elif len(distance_between_flankings_list) > len(del_flanking_dist_list) and len(del_flanking_dist_list) > 0:
            major_classification = 'DEL'
            minor_classification = 'Heterozygous DEL'
        elif len(distance_between_flankings_list) == 1 and len(del_flanking_dist_list) == 0:
            major_classification = 'False DEL'
            minor_classification = 'Potential SV'
        else:
            major_classification = 'False DEL'
            minor_classification = 'False DEL'
            interval_between_matching_flankings_list = [interval for interval in interval_between_matching_flankings_list if interval != None]
            # print(interval_between_matching_flankings_list)
            print('----------- Checking DEL Sequence Alignment in HG2 -----------')
            for distance_interval in interval_between_matching_flankings_list:
                distance_interval_chr, distance_interval_start, distance_interval_end = parse_interval(distance_interval)
                for del_hg2_seq in  del_seq_hg2_alignment_intervals:
                    del_hg2_chr, del_hg2_start, del_hg2_end = parse_interval(del_hg2_seq)
                    if distance_interval_chr == del_hg2_chr:
                        if min(interval_overlapping_percentage(distance_interval, del_hg2_seq)) >0:
                            overlap_to_distance_between_flankings, overlap_to_del_seq = interval_overlapping_percentage(
                                distance_interval, del_hg2_seq)
                            print(
                                f"DEL HG2 alignment:{del_hg2_seq}; Interval Between HG2 Aligned Flankings:{distance_interval}\nThe overlap between the DEL sequence and the interval between the HG2 aligned flankings is {overlap_to_del_seq} of the DEL sequence and {overlap_to_distance_between_flankings} of the interval between the HG2 aligned flankings")
                            # print(f"DEL HG2 alignment:{del_hg2_seq}; Interval Between HG2 Aligned Flankings:{distance_interval}; Overlapping:{interval_overlapping_percentage(distance_interval, del_hg2_seq)}")
            # else:
            #     major_classification = 'Unknown'
            #     minor_classification = 'Unknown'

    # Add the classification to the dictionary
    del_flankings_sum_dict['classification'] = major_classification
    del_flankings_sum_dict['minor_classification'] = minor_classification
    return del_flankings_sum_dict

def plot_del_flankings(del_interval, calling_reference_fasta: str, called_ref_aligner, truth_ref_aligner, plot_ratio=70, save_plot=False, save_plot_path=None):
    del_flanking_dict = collect_del_flankings(del_interval, calling_reference_fasta, called_ref_aligner, truth_ref_aligner)
    hg38_flanking_intervals = del_flanking_dict['hg38_plotting_flanking_intervals']
    hg2_flanking_alignment_intervals = del_flanking_dict['hg2_plotting_flanking_intervals']
    avu.PlotIntervals(hg38_flanking_intervals, hg2_flanking_alignment_intervals).plot_intervals_comparison(flanking=True, ratio=plot_ratio, save=save_plot, savepath=save_plot_path)




