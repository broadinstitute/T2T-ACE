import numpy as np
import pandas as pd
from .validator import align_interval, align_interval_no_restriction, reportAlignment, classifyDupInterval
from .interval_parsing import parse_interval, interval_size, create_interval
from . import alignment_visualization_utilities as avu
from .dup_basepair_correction import extend_2_left, extend_2_right
import matplotlib.pyplot as plt
import os

# Load unresolved intervals within hg38
current_dir = os.path.dirname(os.path.abspath(__file__))
ucsc_gap_filepath = os.path.join(current_dir, 'resources/ucsc_hg38_gap_region.txt')
gap_df = pd.read_csv(ucsc_gap_filepath, sep="\t", index_col=[0],header=None)
# Gather gap intervals from gap_df
gap_intervals = [f"{row[1]}:{row[2]}-{row[3]}" for index,row in gap_df.iterrows()]

def gap_percentage(dup_interval):
    chr, pos, end = parse_interval(dup_interval)
    dup_interval_size = interval_size(dup_interval)
    # Collect all the gaps that have overlap with the DUP interval
    overlap_gap_intervals = []
    for gap_interval in gap_intervals:
        gap_chr, gap_pos, gap_end = parse_interval(gap_interval)
        if chr == gap_chr and pos <= gap_pos and end >= gap_end:
            overlap_gap_intervals.append(gap_interval)
    overlap_gap_intervals_size = sum([interval_size(i) for i in overlap_gap_intervals])
    gap_overlap_percentage = np.round((overlap_gap_intervals_size / dup_interval_size) * 100, 2)
    return gap_overlap_percentage
# This class is to evaluate the DUP intervals
class eval_dup_interval:
    def __init__(self, dup_interval, calling_reference_fasta, truth_ref_fasta, called_ref_aligner, truth_ref_aligner):
        self.dup_interval = dup_interval
        self.calling_reference_fasta = calling_reference_fasta
        self.truth_ref_fasta = truth_ref_fasta
        self.called_ref_aligner = called_ref_aligner
        self.truth_ref_aligner = truth_ref_aligner

    # Plot the alignment of the DUP interval to the two references (hg38 and HG2)
    def plot_alignment(self, ratio=15, save=False, save_path=None):
        hg38_alignments, hg2_alignments = align_interval(self.dup_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)
        hg38_alignment_intervals = [i[0] for i in hg38_alignments]
        hg2_alignment_intervals = [i[0] for i in hg2_alignments]
        avu.PlotIntervals(hg38_alignment_intervals, hg2_alignment_intervals).plot_intervals_comparison(flanking=False, ratio=ratio, save=save, savepath=save_path)

    # This function is to create sliding window alignment for DUPs
    # I am keeping the binning related alignment functions in this class for now
    def bin_alignment(self, window_size=1000):
        # Use SVLEN to determine the window size
        svlen = interval_size(self.dup_interval)
        if svlen < 10000:
            window_size = int(np.round(svlen / 15))
        else:
            window_size = window_size

        sliding_window_interval_list = []
        sliding_window_hg38_hit_list = []
        sliding_window_hg2_mat_hit_list = []
        sliding_window_hg2_pat_hit_list = []
        sliding_window_hg2_hit_list = []
        chr, pos, end = parse_interval(self.dup_interval)

        # Make sure the last sliding window is not out of the original interval
        while pos <= end - window_size:
            sliding_interval = create_interval(chr, pos, pos + window_size)
            pos += window_size
            sliding_window_interval_list.append(sliding_interval)
            # Number of copies in hg38
            sliding_window_hg38_hit_list.append(
                len(align_interval(sliding_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)[0]))
            # Number of copies in HG2-T2T MAT
            sliding_window_hg2_mat_hit_list.append(
                len([i for i, j, s, e in align_interval(sliding_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)[1] if
                     'MAT' in i]))
            # Number of copies in HG2-T2T PAT
            sliding_window_hg2_pat_hit_list.append(
                len([i for i, j, s, e in align_interval(sliding_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)[1] if
                     'PAT' in i]))
            # Number of copies in HG2-T2T
            sliding_window_hg2_hit_list.append(
                len(align_interval(sliding_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)[1]))

        return window_size, sliding_window_interval_list, sliding_window_hg38_hit_list, sliding_window_hg2_mat_hit_list, sliding_window_hg2_pat_hit_list, sliding_window_hg2_hit_list

    def classify_bin_alignment(self):
        window_size, sliding_window_interval_list, sliding_window_hg38_hit_list, sliding_window_hg2_mat_hit_list, sliding_window_hg2_pat_hit_list, sliding_window_hg2_hit_list = self.bin_alignment()
        sliding_window_classification_list = []
        for i in sliding_window_interval_list:
            if sliding_window_hg38_hit_list[sliding_window_interval_list.index(i)] < sliding_window_hg2_mat_hit_list[sliding_window_interval_list.index(i)] and sliding_window_hg38_hit_list[sliding_window_interval_list.index(i)] < sliding_window_hg2_pat_hit_list[sliding_window_interval_list.index(i)]:
                 sliding_window_classification_list.append("Homozygous Duplication")
            elif sliding_window_hg38_hit_list[sliding_window_interval_list.index(i)] < sliding_window_hg2_mat_hit_list[sliding_window_interval_list.index(i)] and sliding_window_hg38_hit_list[sliding_window_interval_list.index(i)] >= sliding_window_hg2_pat_hit_list[sliding_window_interval_list.index(i)]:
                 sliding_window_classification_list.append("Maternal Heterozygous Duplication")
            elif sliding_window_hg38_hit_list[sliding_window_interval_list.index(i)] >= sliding_window_hg2_mat_hit_list[sliding_window_interval_list.index(i)] and sliding_window_hg38_hit_list[sliding_window_interval_list.index(i)] < sliding_window_hg2_pat_hit_list[sliding_window_interval_list.index(i)]:
                 sliding_window_classification_list.append("Paternal Heterozygous Duplication")
            else:
                 sliding_window_classification_list.append("False Duplication")

        print(f"input dup interval: {self.dup_interval}")
        print(f"interval size: {interval_size(self.dup_interval)}")
        print(f"window size: {window_size}, number of windows: {len(sliding_window_interval_list)}")

        bin_classification_dict = {}
        for i in set(sliding_window_classification_list):
            bin_classification_dict[i] = sliding_window_classification_list.count(i)
            classification_length = sliding_window_classification_list.count(i)*window_size
            classification_percentage = np.round(100*(classification_length/interval_size(self.dup_interval)),2)
            print(f"{i}: {classification_percentage}%")

        return bin_classification_dict

    def plot_bin_alignment(self, save=False):
        # Get the number of copies of the whole interval
        whole_interval_hg38_hit = len(align_interval(self.dup_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)[0])
        whole_interval_hg2_mat_hit = len(
            [i for i, j in align_interval(self.dup_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)[1] if 'MAT' in i])
        whole_interval_hg2_pat_hit = len(
            [i for i, j in align_interval(self.dup_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)[1] if 'PAT' in i])
        whole_interval_hg2_hit = len(align_interval(self.dup_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)[1])
        window_size, sliding_window_interval_list, sliding_window_hg38_hit_list, sliding_window_hg2_mat_hit_list, sliding_window_hg2_pat_hit_list, sliding_window_hg2_hit_list = self.bin_alignment()

        fig, axs = plt.subplots(1, 4, figsize=(12, 3), sharey=True)
        for i in sliding_window_interval_list:
            if sliding_window_hg38_hit_list[sliding_window_interval_list.index(i)] < sliding_window_hg2_mat_hit_list[
                sliding_window_interval_list.index(i)] and sliding_window_hg38_hit_list[
                sliding_window_interval_list.index(i)] < sliding_window_hg2_pat_hit_list[
                sliding_window_interval_list.index(i)]:
                axs[0].plot(i, sliding_window_hg38_hit_list[sliding_window_interval_list.index(i)], '*', color='blue')
                axs[1].plot(i, sliding_window_hg2_mat_hit_list[sliding_window_interval_list.index(i)], '*', color='red')
                axs[2].plot(i, sliding_window_hg2_pat_hit_list[sliding_window_interval_list.index(i)], '*', color='red')
                axs[3].plot(i, sliding_window_hg2_hit_list[sliding_window_interval_list.index(i)], '*', color='red')
            elif sliding_window_hg38_hit_list[sliding_window_interval_list.index(i)] < sliding_window_hg2_mat_hit_list[
                sliding_window_interval_list.index(i)] and sliding_window_hg38_hit_list[
                sliding_window_interval_list.index(i)] >= sliding_window_hg2_pat_hit_list[
                sliding_window_interval_list.index(i)]:
                axs[0].plot(i, sliding_window_hg38_hit_list[sliding_window_interval_list.index(i)], '*', color='blue')
                axs[1].plot(i, sliding_window_hg2_mat_hit_list[sliding_window_interval_list.index(i)], '*', color='red')
                axs[2].plot(i, sliding_window_hg2_pat_hit_list[sliding_window_interval_list.index(i)], '*',
                            color='grey')
                axs[3].plot(i, sliding_window_hg2_hit_list[sliding_window_interval_list.index(i)], '*', color='red')
            elif sliding_window_hg38_hit_list[sliding_window_interval_list.index(i)] >= sliding_window_hg2_mat_hit_list[
                sliding_window_interval_list.index(i)] and sliding_window_hg38_hit_list[
                sliding_window_interval_list.index(i)] < sliding_window_hg2_pat_hit_list[
                sliding_window_interval_list.index(i)]:
                axs[0].plot(i, sliding_window_hg38_hit_list[sliding_window_interval_list.index(i)], '*', color='blue')
                axs[1].plot(i, sliding_window_hg2_mat_hit_list[sliding_window_interval_list.index(i)], '*',
                            color='grey')
                axs[2].plot(i, sliding_window_hg2_pat_hit_list[sliding_window_interval_list.index(i)], '*', color='red')
                axs[3].plot(i, sliding_window_hg2_hit_list[sliding_window_interval_list.index(i)], '*', color='red')
            else:
                axs[0].plot(i, sliding_window_hg38_hit_list[sliding_window_interval_list.index(i)], '*', color='grey')
                axs[1].plot(i, sliding_window_hg2_mat_hit_list[sliding_window_interval_list.index(i)], '*',
                            color='grey')
                axs[2].plot(i, sliding_window_hg2_pat_hit_list[sliding_window_interval_list.index(i)], '*',
                            color='grey')
                axs[3].plot(i, sliding_window_hg2_hit_list[sliding_window_interval_list.index(i)], '*', color='grey')

        axs[0].set_title(f'hg38 ({whole_interval_hg38_hit})')
        axs[0].set_ylim(0, 10)
        axs[0].set_xticks('')

        axs[1].set_title(f'HG2-MAT ({whole_interval_hg2_mat_hit})')
        axs[1].set_ylim(0, 10)
        axs[1].set_xticks('')

        axs[2].set_title(f'HG2-PAT ({whole_interval_hg2_pat_hit})')
        axs[2].set_ylim(0, 10)
        axs[2].set_xticks('')

        axs[3].set_title(f'HG2-T2T ({whole_interval_hg2_hit})')
        axs[3].set_ylim(0, 10)
        axs[3].set_xticks('')

        fig.suptitle(f'{self.dup_interval} ({interval_size(self.dup_interval)}bp)')
        plt.tight_layout()
        if save:
            plt.savefig(f'{self.dup_interval}.png', dpi=300)

    def analyzeDuplicationIntervals(self):
        dup_summary_dict = {}
        print("----------------------------------------------------------------")
        print("Alignment Report of Original Interval")
        print("----------------------------------------------------------------")

        chr, pos, end = parse_interval(self.dup_interval)
        dup_interval_size = interval_size(self.dup_interval)
        if dup_interval_size > 1000000:
            print("DUP interval is too large")
            return None

        # Check if the DUP interval contains more than 10% unresolved region in hg38
        big_gap_dup = False
        # If the overlap is more than 10%, then the DUP interval will be aligned without any restriction
        # If the overlap is less than 10%, then the DUP interval will be aligned with the alignment length restriction
        if gap_percentage(self.dup_interval) >= 10:
            big_gap_dup = True
            print(f"The DUP interval {self.dup_interval} contains {gap_percentage(self.dup_interval)}% unresolved region in hg38")
            hg38_dup_count, hg2_mat_count, hg2_pat_count, hg2_dup_count = reportAlignment(self.dup_interval,
                                                                                          self.calling_reference_fasta,
                                                                                          self.called_ref_aligner,
                                                                                          self.truth_ref_aligner,
                                                                                          print_alignments=True,
                                                                                          alignment_filter=False)
            # Get alignment output of the gap contained interval (without alignment length restriction)
            hg38_alignments, hg2_alignments = align_interval_no_restriction(self.dup_interval, self.calling_reference_fasta,
                                                             self.called_ref_aligner, self.truth_ref_aligner)
            # Classify the DUP interval
            dup_interval_classification = classifyDupInterval(self.dup_interval, self.calling_reference_fasta,
                                                                  self.called_ref_aligner, self.truth_ref_aligner,
                                                                  alignment_filter=False)
            major_classification = dup_interval_classification[0]
            sub_classification = dup_interval_classification[1]
        else:
            hg38_dup_count, hg2_mat_count, hg2_pat_count, hg2_dup_count = reportAlignment(self.dup_interval,
                                                                                          self.calling_reference_fasta,
                                                                                          self.called_ref_aligner,
                                                                                          self.truth_ref_aligner,
                                                                                          print_alignments=True,
                                                                                          alignment_filter=True)
            # Get alignment output of the original interval (with alignment length restriction)
            hg38_alignments, hg2_alignments = align_interval(self.dup_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)
            # Classify the DUP interval (with alignment length restriction)
            dup_interval_classification = classifyDupInterval(self.dup_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner, alignment_filter=True)
            major_classification = dup_interval_classification[0]
            sub_classification = dup_interval_classification[1]

        # Save the alignment report of the original interval to the summary dictionary
        dup_summary_dict['dup_interval'] = self.dup_interval
        dup_summary_dict['dup_interval_size'] = dup_interval_size
        dup_summary_dict['original_hg38_hit_count'] = hg38_dup_count
        dup_summary_dict['original_hg2_hit_count'] = hg2_dup_count
        dup_summary_dict['original_hg2_mat_hit_count'] = hg2_mat_count
        dup_summary_dict['original_hg2_pat_hit_count'] = hg2_pat_count
        dup_summary_dict['original_dup_interval_major_classification'] = major_classification
        dup_summary_dict['original_dup_interval_sub_classification'] = sub_classification
        dup_summary_dict['original_dup_interval_contain_big_gap'] = big_gap_dup

        if major_classification == "Copy Neutral" or major_classification == "Unknown" or major_classification == "Reference Error" or big_gap_dup == True or dup_interval_size > 1000000:
            print(f"The DUP interval's classification is {major_classification} and it will not be corrected")
            # Assign the attributes associated with corrected interval to NA
            dup_summary_dict['corrected_interval'] = np.nan
            dup_summary_dict['corrected_interval_size'] = np.nan
            dup_summary_dict['expanded_length'] = np.nan
            dup_summary_dict['corrected_hg38_hit_count'] = np.nan
            dup_summary_dict['corrected_hg2_hit_count'] = np.nan
            dup_summary_dict['corrected_hg2_mat_hit_count'] = np.nan
            dup_summary_dict['corrected_hg2_pat_hit_count'] = np.nan
            dup_summary_dict['corrected_dup_interval_major_classification'] = np.nan
            dup_summary_dict['corrected_dup_interval_sub_classification'] = np.nan
            dup_summary_dict['corrected_pos_accuracy'] = np.nan
            dup_summary_dict['corrected_end_accuracy'] = np.nan
        # If the DUP interval is not a false duplication, then it might need to be corrected
        else:
            # Use hg2 alignments starts and ends to determine if the DUP interval needs to be adjusted
            hg2_alignments_starts = [i[2] for i in hg2_alignments]
            hg2_alignments_ends = [i[3] for i in hg2_alignments]
            left_side_correction = ""
            right_side_correction = ""
            left_basepair_accuracy = ""
            right_basepair_accuracy = ""
            print("----------------------------------------------------------------")
            print(f'Analyzing the DUP interval {self.dup_interval}({dup_interval_size}bp)')
            print("----------------------------------------------------------------")
            print("ANALYZING THE POS OF THE DUP INTERVAL")
            print("----------------------------------------------------------------")
            if max(hg2_alignments_starts)>1:
                print(f"The called interval pos needs to be moved to the right by {max(hg2_alignments_starts)}bp")
                left_side_correction = pos + max(hg2_alignments_starts)
                left_basepair_accuracy = 1
            elif min(hg2_alignments_starts)<=1:
                print("The called interval pos needs to be checked by extend_2_left")
                print("**Using extend_2_left to check the pos**")
                interval_2_extend = create_interval(chr, pos + min(hg2_alignments_starts), end)
                extend_2_left_interval, left_basepair_accuracy = extend_2_left(interval_2_extend, self.calling_reference_fasta, self.called_ref_aligner,
                              self.truth_ref_aligner)
                extend_2_left_chr, extend_2_left_pos, extend_2_left_end = parse_interval(extend_2_left_interval)
                left_side_correction = int(extend_2_left_pos)
            print("----------------------------------------------------------------")
            print("ANALYZING THE END OF THE DUP INTERVAL")
            print("----------------------------------------------------------------")
            if min(hg2_alignments_ends)<dup_interval_size and hg2_alignments_ends.count(max(hg2_alignments_ends)) <= hg38_dup_count:
                print(f"The called interval end needs to be moved to the left by {max(hg2_alignments_ends)-min(hg2_alignments_ends)}bp")
                right_side_correction = pos + min(hg2_alignments_ends)
                right_basepair_accuracy = 1
            # If there are enough alignments to support the original interval end, then the end of the interval will not be moved
            elif min(hg2_alignments_ends)<dup_interval_size and hg2_alignments_ends.count(max(hg2_alignments_ends)) > hg38_dup_count:
                print(f"No need to move the called interval end to the left. Left movement might cause more alignments than expected")
                right_side_correction = pos + max(hg2_alignments_ends)
                right_basepair_accuracy = 1
            elif max(hg2_alignments_ends)>=dup_interval_size:
                print("The called interval end needs to be checked by extend_2_right")
                print("**Using extend_2_right to check the end**")
                interval_2_extend = create_interval(chr, pos, pos + max(hg2_alignments_ends))
                extended_2_right_interval, right_basepair_accuracy = extend_2_right(interval_2_extend, self.calling_reference_fasta, self.called_ref_aligner,
                               self.truth_ref_aligner)
                extend_2_right_chr, extend_2_right_pos, extend_2_right_end = parse_interval(extended_2_right_interval)
                right_side_correction = int(extend_2_right_end)

            corrected_interval = create_interval(chr, left_side_correction, right_side_correction)
            corrected_interval_size = interval_size(corrected_interval)
            print("----------------------------------------------------------------")
            print(f"Corrected Interval Details")
            print("----------------------------------------------------------------")
            print(f"The original interval is {self.dup_interval} ({dup_interval_size}bp)")
            print(f"The corrected interval is {corrected_interval} ({corrected_interval_size}bp)")
            extended_length = corrected_interval_size - dup_interval_size
            print(f"The corrected interval is {extended_length}bp longer than the original interval")
            print("----------------------------------------------------------------")
            print("Alignment Report of Corrected Interval")
            print("----------------------------------------------------------------")
            corrected_hg38_dup_count, corrected_hg2_mat_count, corrected_hg2_pat_count, corrected_hg2_dup_count = reportAlignment(corrected_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner, print_alignments=True)
            dup_summary_dict['corrected_interval'] = corrected_interval
            dup_summary_dict['corrected_interval_size'] = corrected_interval_size
            dup_summary_dict['expanded_length'] = extended_length
            refined_hg38_alignments, refined_hg2_alignments = align_interval(corrected_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)
            corrected_dup_interval_classification = classifyDupInterval(corrected_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)
            corrected_interval_major_classification = corrected_dup_interval_classification[0]
            corrected_interval_sub_classification = corrected_dup_interval_classification[1]
            dup_summary_dict['corrected_hg38_hit_count'] = corrected_hg38_dup_count
            dup_summary_dict['corrected_hg2_hit_count'] = corrected_hg2_dup_count
            dup_summary_dict['corrected_hg2_mat_hit_count'] = corrected_hg2_mat_count
            dup_summary_dict['corrected_hg2_pat_hit_count'] = corrected_hg2_pat_count
            dup_summary_dict['corrected_dup_interval_major_classification'] = corrected_interval_major_classification
            dup_summary_dict['corrected_dup_interval_sub_classification'] = corrected_interval_sub_classification
            dup_summary_dict['corrected_pos_accuracy'] = left_basepair_accuracy
            dup_summary_dict['corrected_end_accuracy'] = right_basepair_accuracy
        print("----------------------------------------------------------------")
        print(f"End of Analysis of DUP interval {self.dup_interval}({dup_interval_size}bp)")
        print("----------------------------------------------------------------")
        return dup_summary_dict












