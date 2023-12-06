import numpy as np
from T2T_ACE.validator import align_interval, check_pairwise_alignment
from T2T_ACE.interval_parsing import parse_interval, interval_size, create_interval
import T2T_ACE.alignment_visualization_utilities as avu
import matplotlib.pyplot as plt

class eval_dup_interval:
    def __init__(self, dup_interval, calling_reference_fasta, truth_ref_fasta, called_ref_aligner, truth_ref_aligner):
        self.dup_interval = dup_interval
        self.calling_reference_fasta = calling_reference_fasta
        self.truth_ref_fasta = truth_ref_fasta
        self.called_ref_aligner = called_ref_aligner
        self.truth_ref_aligner = truth_ref_aligner

    # Report the alignment of the DUP interval to the two references (hg38 and HG2)
    def reportAlignment(self):
        # First Align the DUP sequence to HG2 and hg38
        dup_alignments = align_interval(self.dup_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)

        # Check the number of alignments for DUP and DEL in HG2 and hg38
        hg38_dup_count = len(dup_alignments[0])
        hg2_dup_count = len(dup_alignments[1])
        hg2_mat_copies = [[i,j] for i, j in dup_alignments[1] if 'MAT' in i]
        hg2_pat_copies = [[i,j] for i, j in dup_alignments[1] if 'PAT' in i]

        print(f"input dup interval: {self.dup_interval}")
        print(f"hg38 dup count: {hg38_dup_count}")
        for i, j in dup_alignments[0]:
            print(f"interval: {i}\tstrand: {j}")
        print(f"hg2 dup count: {hg2_dup_count}")
        for i, j in hg2_mat_copies:
            print(f"interval: {i}\tstrand: {j}")
        for i, j in hg2_pat_copies:
            print(f"interval: {i}\tstrand: {j}")

        return(hg38_dup_count, len(hg2_mat_copies), len(hg2_pat_copies), hg2_dup_count)

    # Classify the interval as either a true duplication or a false duplication
    def classifyInterval(self):
        hg38_dup_count, hg2_mat_count, hg2_pat_count, hg2_dup_count = self.reportAlignment()

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
                len([i for i, j in align_interval(sliding_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)[1] if
                     'MAT' in i]))
            # Number of copies in HG2-T2T PAT
            sliding_window_hg2_pat_hit_list.append(
                len([i for i, j in align_interval(sliding_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)[1] if
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

    def extend_interval(self, extend_size=1000):
        window_size, sliding_window_interval_list, sliding_window_hg38_hit_list, sliding_window_hg2_mat_hit_list, sliding_window_hg2_pat_hit_list, sliding_window_hg2_hit_list = self.bin_alignment()
        # Check the duplication status of the leftest window
        if sliding_window_hg38_hit_list[0]*2 >= sliding_window_hg2_hit_list[0]:
            print("The leftest window is a false duplication, no need to extend the interval")
        elif sliding_window_hg38_hit_list[0]*2 < sliding_window_hg2_hit_list[0]:
            print("The leftest window is a real duplication, extend the interval to the left one window at a time")
            # Extend the interval to the left by one window
            leftest_interval = sliding_window_interval_list[0]
            leftest_chr, leftest_pos, leftest_end = parse_interval(leftest_interval)
            leftest_interval_hg38_hits = len(align_interval(leftest_interval, self.calling_reference_fasta, self.called_ref_aligner,
                                   self.truth_ref_aligner)[0])
            leftest_interval_hg2_hits = len(align_interval(leftest_interval, self.calling_reference_fasta, self.called_ref_aligner,
                                   self.truth_ref_aligner)[1])
            extended_leftest_interval = ""
            extended_leftest_hg38_hits = 0
            extended_leftest_hg2_hits = 0
            while leftest_interval_hg38_hits*2 < leftest_interval_hg2_hits:
                leftest_pos = leftest_pos - window_size
                leftest_end = leftest_end - window_size
                new_leftest_interval = create_interval(leftest_chr, leftest_pos, leftest_end)
                leftest_interval_hg38_hits = len(align_interval(new_leftest_interval, self.calling_reference_fasta, self.called_ref_aligner,
                                   self.truth_ref_aligner)[0])
                leftest_interval_hg2_hits = len(align_interval(new_leftest_interval, self.calling_reference_fasta, self.called_ref_aligner,
                                   self.truth_ref_aligner)[1])
                print(f"new leftest interval: {new_leftest_interval}, {leftest_interval_hg38_hits}, {leftest_interval_hg2_hits}")
                if leftest_interval_hg2_hits*2 <= leftest_interval_hg2_hits:
                    break
                extended_leftest_interval = new_leftest_interval
                extended_leftest_hg38_hits = leftest_interval_hg38_hits
                extended_leftest_hg2_hits = leftest_interval_hg2_hits
            print(f"extended leftest interval: {extended_leftest_interval}, {extended_leftest_hg38_hits}, {extended_leftest_hg2_hits}")


            # chr, pos, end = parse_interval(self.dup_interval)
            # hg38_hits = len(align_interval(self.dup_interval, self.calling_reference_fasta, self.called_ref_aligner,
            #                    self.truth_ref_aligner)[0])
            # hg2_hits = len(align_interval(self.dup_interval, self.calling_reference_fasta, self.called_ref_aligner,
            #                    self.truth_ref_aligner)[1])
            # print(f"Initial dup interval: {self.dup_interval}, {hg38_hits}, {hg2_hits}")
            #
            # extended_interval = ""
            # extended_hg38_hits = 0
            # extended_hg2_hits = 0
            # while hg38_hits*2 < hg2_hits:
            #     pos = pos - window_size
            #     new_interval = create_interval(chr, pos, end)
            #     hg38_hits = len(align_interval(new_interval, self.calling_reference_fasta, self.called_ref_aligner,
            #                        self.truth_ref_aligner)[0])
            #     hg2_hits = len(align_interval(new_interval, self.calling_reference_fasta, self.called_ref_aligner,
            #                        self.truth_ref_aligner)[1])
            #     if hg38_hits*2 == hg2_hits:
            #         break
            #     extended_interval = new_interval
            #     extended_hg38_hits = hg38_hits
            #     extended_hg2_hits = hg2_hits
            # print(f"extended dup interval: {extended_interval}, {extended_hg38_hits}, {extended_hg2_hits}")
            # extended_chr, extended_pos, extended_end = parse_interval(extended_interval)
            # initial_chr, initial_pos, initial_end = parse_interval(self.dup_interval)
            # # Create a new interval describing the left extending region
            # left_extension_interval = create_interval(extended_chr, extended_pos, initial_pos)
            # left_extension_interval_hg38_hits = len(align_interval(left_extension_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)[0])
            # left_extension_interval_hg2_hits = len(align_interval(left_extension_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)[1])
            # print(f"left extension interval: {left_extension_interval}, size: {interval_size(left_extension_interval)}")
            # print(f"left extension interval hits: {left_extension_interval_hg38_hits}, {left_extension_interval_hg2_hits}")
            #left_extension_pairwise_alignment_dict = check_pairwise_alignment(left_extension_interval, self.calling_reference_fasta, self.truth_ref_fasta, self.called_ref_aligner, self.truth_ref_aligner)
            #print(f"Pairwise alignment dict: {left_extension_pairwise_alignment_dict}")

                #left_extension_interval = create_interval(chr, pos, end)
                #print(f"extended dup interval: {dup_interval}, {hg38_hits}, {hg2_hits}")
            #print(f"Extended dup interval: {dup_interval}")
            #pairwise_alignment_dict = check_pairwise_alignment(dup_interval, self.calling_reference_fasta, self.truth_ref_fasta, self.called_ref_aligner, self.truth_ref_aligner)
            #print(f"Pairwise alignment dict: {pairwise_alignment_dict}")
            #print(f"FINAL extended dup interval: {dup_interval}")
            # new_pos = pos - window_size
            # # Stitch the leftest window with the original leftest window
            # stitched_left_interval = create_interval(chr, new_pos, pos + window_size)
            # print(f"Stitched left interval: {stitched_left_interval}")
            # # Check the duplication status of the new leftest window
            # stitched_left_interval_hg38_hits = len(align_interval(stitched_left_interval, self.calling_reference_fasta, self.called_ref_aligner,
            #                    self.truth_ref_aligner)[0])
            # print('stitched left interval hits:',stitched_left_interval_hg38_hits)
            # if stitched_left_interval_hg38_hits == 1:
            #     stitched_left_interval_dict = check_pairwise_alignment(stitched_left_interval, self.calling_reference_fasta, self.truth_ref_fasta, self.called_ref_aligner, self.truth_ref_aligner)
            #     print('stitched left interval dict:',stitched_left_interval_dict)












