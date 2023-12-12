import numpy as np
from T2T_ACE.validator import align_interval, reportAlignment
from T2T_ACE.interval_parsing import parse_interval, interval_size, create_interval
import T2T_ACE.alignment_visualization_utilities as avu
from T2T_ACE.dup_basepair_correction import extend_2_left, extend_2_right
import matplotlib.pyplot as plt

class eval_dup_interval:
    def __init__(self, dup_interval, calling_reference_fasta, truth_ref_fasta, called_ref_aligner, truth_ref_aligner):
        self.dup_interval = dup_interval
        self.calling_reference_fasta = calling_reference_fasta
        self.truth_ref_fasta = truth_ref_fasta
        self.called_ref_aligner = called_ref_aligner
        self.truth_ref_aligner = truth_ref_aligner

    # Classify the interval as either a true duplication or a false duplication
    def classifyInterval(self):
        hg38_dup_count, hg2_mat_count, hg2_pat_count, hg2_dup_count = reportAlignment(self.dup_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)

        classification = ""
        if hg2_mat_count > hg38_dup_count and hg2_pat_count > hg38_dup_count:
            classification = "Homozygous Duplication"
        elif hg2_mat_count > hg38_dup_count and hg2_pat_count == hg38_dup_count:
            classification = "Maternal Heterozygous Duplication"
        elif hg2_pat_count > hg38_dup_count and hg2_mat_count == hg38_dup_count:
            classification = "Paternal Heterozygous Duplication"
        else:
            classification = "False Duplication"
        return classification

    # Plot the alignment of the DUP interval to the two references (hg38 and HG2)
    def plot_alignment(self, ratio=15, save=False):
        hg38_alignments, hg2_alignments = align_interval(self.dup_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)
        hg38_alignment_intervals = [i[0] for i in hg38_alignments]
        hg2_alignment_intervals = [i[0] for i in hg2_alignments]
        avu.PlotIntervals(hg38_alignment_intervals, hg2_alignment_intervals).plot_intervals_comparison(flanking=False, ratio=ratio, save=save)

    # This function is to create sliding window alignment for DUPs
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
        reportAlignment(self.dup_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)
        # Analyze the duplication interval and check its relationship with the real event
        hg38_alignments, hg2_alignments = align_interval(self.dup_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)
        # Save the alignment report of the original interval to the summary dictionary
        dup_summary_dict['dup_interval'] = self.dup_interval
        dup_summary_dict['dup_interval_size'] = interval_size(self.dup_interval)
        dup_summary_dict['original_hg38_hit_count'] = len(hg38_alignments)
        dup_summary_dict['original_hg2_hit_count'] = len(hg2_alignments)
        if len(hg38_alignments) > 1:
            print(f"WARNING: the DUP interval {self.dup_interval} has more than one alignment in hg38")
        chr, pos, end = parse_interval(self.dup_interval)
        dup_interval_size = interval_size(self.dup_interval)
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
            extend_2_left_chr, extend_2_left_pos, extend_2_left_end = parse_interval(extend_2_left(interval_2_extend, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)[0])
            left_basepair_accuracy = extend_2_left(interval_2_extend, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)[1]
            left_side_correction = int(extend_2_left_pos)
        print("----------------------------------------------------------------")
        print("ANALYZING THE END OF THE DUP INTERVAL")
        print("----------------------------------------------------------------")
        if min(hg2_alignments_ends)<dup_interval_size:
            print(f"The called interval end needs to be moved to the left by {max(hg2_alignments_ends)-min(hg2_alignments_ends)}bp")
            right_side_correction = pos + min(hg2_alignments_ends)
            right_basepair_accuracy = 1
        elif max(hg2_alignments_ends)>=dup_interval_size:
            print("The called interval end needs to be checked by extend_2_right")
            print("**Using extend_2_right to check the end**")
            interval_2_extend = create_interval(chr, pos, pos + max(hg2_alignments_ends))
            extend_2_right_chr, extend_2_right_pos, extend_2_right_end = parse_interval(extend_2_right(interval_2_extend, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)[0])
            right_basepair_accuracy = extend_2_right(interval_2_extend, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)[1]
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
        reportAlignment(corrected_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)
        dup_summary_dict['corrected_interval'] = corrected_interval
        dup_summary_dict['corrected_interval_size'] = corrected_interval_size
        dup_summary_dict['expanded_length'] = extended_length
        refined_hg38_alignments, refined_hg2_alignments = align_interval(corrected_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)
        dup_summary_dict['corrected_hg38_hit_count'] = len(refined_hg38_alignments)
        dup_summary_dict['corrected_hg2_hit_count'] = len(refined_hg2_alignments)
        dup_summary_dict['corrected_pos_accuracy'] = left_basepair_accuracy
        dup_summary_dict['corrected_end_accuracy'] = right_basepair_accuracy
        return dup_summary_dict

    # This function will determine if the DUP interval needs further extension to match with reality
    # def extend_interval(self, extension_alignment=False):
    #     window_size, sliding_window_interval_list, sliding_window_hg38_hit_list, sliding_window_hg2_mat_hit_list, sliding_window_hg2_pat_hit_list, sliding_window_hg2_hit_list = self.bin_alignment()
    #     # Check the duplication status of the leftest window
    #     if sliding_window_hg38_hit_list[0]*2 >= sliding_window_hg2_hit_list[0]:
    #         print(f"The leftest window {sliding_window_interval_list[0]} is a false duplication, no need to extend the interval")
    #     elif sliding_window_hg38_hit_list[0]*2 < sliding_window_hg2_hit_list[0]:
    #         print(f"The leftest window {sliding_window_interval_list[0]} is a real duplication, extend the interval to the left one window at a time")
    #         # Extend the interval to the left by one window
    #         leftest_interval = sliding_window_interval_list[0]
    #         leftest_chr, leftest_pos, leftest_end = parse_interval(leftest_interval)
    #         leftest_interval_hg38_copies = align_interval(leftest_interval, self.calling_reference_fasta, self.called_ref_aligner,
    #                        self.truth_ref_aligner)[0]
    #         leftest_interval_hg38_hits = len(leftest_interval_hg38_copies)
    #         leftest_interval_hg2_copies = align_interval(leftest_interval, self.calling_reference_fasta, self.called_ref_aligner,
    #                        self.truth_ref_aligner)[1]
    #         leftest_interval_hg2_hits = len(leftest_interval_hg2_copies)
    #         print(f"leftest interval: {leftest_interval},{leftest_interval_hg38_hits}, {leftest_interval_hg2_hits}")
    #         if extension_alignment:
    #             print(f"{leftest_interval_hg38_copies}, {leftest_interval_hg2_copies}")
    #         extended_leftest_intervals = []
    #         while leftest_interval_hg38_hits*2 < leftest_interval_hg2_hits:
    #             leftest_pos = leftest_pos - window_size
    #             leftest_end = leftest_end - window_size
    #             new_leftest_interval = create_interval(leftest_chr, leftest_pos, leftest_end)
    #             leftest_interval_hg38_copies = align_interval(new_leftest_interval, self.calling_reference_fasta, self.called_ref_aligner,
    #                                self.truth_ref_aligner)[0]
    #             leftest_interval_hg38_hits = len(leftest_interval_hg38_copies)
    #             leftest_interval_hg2_copies = align_interval(new_leftest_interval, self.calling_reference_fasta, self.called_ref_aligner,
    #                                self.truth_ref_aligner)[1]
    #             leftest_interval_hg2_hits = len(leftest_interval_hg2_copies)
    #             print(f"new leftest interval: {new_leftest_interval}, {leftest_interval_hg38_hits}, {leftest_interval_hg2_hits}")
    #             if extension_alignment:
    #                 print(f"{leftest_interval_hg38_copies}, {leftest_interval_hg2_copies}")
    #             extended_leftest_intervals.append(new_leftest_interval)
    #         # Check the duplication status of the rightest window
    #         if sliding_window_hg38_hit_list[-1] * 2 >= sliding_window_hg2_hit_list[-1]:
    #             print(f"The rightest window {sliding_window_interval_list[-1]} is a false duplication, no need to extend the interval")
    #         elif sliding_window_hg38_hit_list[-1] * 2 < sliding_window_hg2_hit_list[-1]:
    #             print(f"The rightest window {sliding_window_interval_list[-1]} is a real duplication, extend the interval to the right one window at a time")
    #             # Extend the interval to the left by one window
    #             rightest_interval = sliding_window_interval_list[-1]
    #             rightest_chr, rightest_pos, rightest_end = parse_interval(rightest_interval)
    #             rightest_interval_hg38_copies = align_interval(rightest_interval, self.calling_reference_fasta, self.called_ref_aligner,
    #                            self.truth_ref_aligner)[0]
    #             rightest_interval_hg38_hits = len(rightest_interval_hg38_copies)
    #             rightest_interval_hg2_copies = align_interval(rightest_interval, self.calling_reference_fasta, self.called_ref_aligner,
    #                            self.truth_ref_aligner)[1]
    #             rightest_interval_hg2_hits = len(rightest_interval_hg2_copies)
    #             print(
    #                 f"rightest interval: {rightest_interval},{rightest_interval_hg38_hits}, {rightest_interval_hg2_hits}")
    #             if extension_alignment:
    #                 print(f"{rightest_interval_hg38_copies}, {rightest_interval_hg2_copies}")
    #             extended_rightest_intervals = []
    #             while rightest_interval_hg38_hits * 2 < rightest_interval_hg2_hits:
    #                 rightest_pos = rightest_pos + window_size
    #                 rightest_end = rightest_end + window_size
    #                 new_rightest_interval = create_interval(rightest_chr, rightest_pos, rightest_end)
    #                 rightest_interval_hg38_copies = align_interval(new_rightest_interval, self.calling_reference_fasta, self.called_ref_aligner,
    #                                self.truth_ref_aligner)[0]
    #                 rightest_interval_hg38_hits = len(rightest_interval_hg38_copies)
    #                 rightest_interval_hg2_copies = align_interval(new_rightest_interval, self.calling_reference_fasta, self.called_ref_aligner,
    #                                self.truth_ref_aligner)[1]
    #                 rightest_interval_hg2_hits = len(rightest_interval_hg2_copies)
    #                 print(
    #                     f"new rightest interval: {new_rightest_interval}, {rightest_interval_hg38_hits}, {rightest_interval_hg2_hits}")
    #                 if extension_alignment:
    #                     print(f"{rightest_interval_hg38_copies}, {rightest_interval_hg2_copies}")
    #                 extended_rightest_intervals.append(new_rightest_interval)












