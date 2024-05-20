import numpy as np
from .validator import align_interval
from .interval_parsing import parse_interval, interval_size, create_interval


# DUP basepair correction function of T2T-ACE
# It will extend the DUP interval to the right until no more copies are found
def extend_2_right(dup_interval, calling_reference_fasta, called_ref_aligner,
                                 truth_ref_aligner):
    chr, pos, end = parse_interval(dup_interval)
    if chr == "chrX" or chr == "chrY":
        print("Basepair correction Method doesn't support DUPs on chrX or chrY")
        return dup_interval, 0
    original_end = end
    # The window size is 10% of the interval size
    window_size = int(np.round(interval_size(dup_interval) * 0.1, 0))
    print("extending window size", window_size)

    q_end_list = [q_en for interval, strand, q_st, q_en in
                  align_interval(dup_interval, calling_reference_fasta, called_ref_aligner,
                                 truth_ref_aligner)[1]]
    original_end_count = len(q_end_list)
    new_interval = dup_interval
    copies_alignment_end = ""
    basepair_accuracy = 0
    while min(q_end_list) <= interval_size(new_interval):
        new_end = end + window_size
        new_interval = create_interval(chr, pos, new_end)
        new_alignments = align_interval(new_interval, calling_reference_fasta, called_ref_aligner,
                                        truth_ref_aligner)

        q_end_list = sorted([q_en for interval, strand, q_st, q_en in new_alignments[1]])
        q_end_count = [[q_end, q_end_list.count(q_end)] for q_end in sorted(set(q_end_list))]

        if q_end_count[-1][1] == 2:
            print("new interval", new_interval, interval_size(new_interval))
            print(q_end_count)
            if len(q_end_count) > 1:
                print(f"extend end to the right by {q_end_count[0][0] - interval_size(dup_interval)}bp")
                copies_alignment_end = pos + q_end_count[0][0]
                basepair_accuracy = 1
            elif len(q_end_count) == 1:
                basepair_accuracy = 0
                print(f"the window size is too big")
                if end > original_end:
                    copies_alignment_end = end
                    print(f"Keep the last iteration pos which means extend end to the right {end - original_end}bp")
                else:
                    copies_alignment_end = original_end
                    print(f"Keep the original end")
            break
        elif len(q_end_count) > 1 and q_end_count[-1][0] == interval_size(new_interval):
            print("new interval", new_interval, interval_size(new_interval))
            print(q_end_count)
            print("The sequence isn't homozygous")
            print(f"extend end to the right by {q_end_count[0][0] - interval_size(dup_interval)}bp")
            copies_alignment_end = pos + q_end_count[0][0]
            basepair_accuracy = 0
            break
        elif max(q_end_list) < interval_size(new_interval):
            print("new interval", new_interval, interval_size(new_interval))
            print(q_end_count)
            print("This sequence might hit the unresolved region in hg38")
            print(f"extend end to the right by {q_end_count[0][0] - interval_size(dup_interval)}bp")
            copies_alignment_end = pos + q_end_count[0][0]
            basepair_accuracy = 0
            break
        elif len(q_end_list) > original_end_count:
            print("new interval", new_interval, interval_size(new_interval))
            print(q_end_count)
            print("The sequence has extended to copies that doesn't exist in the original interval alignment")
            print(f"Keep the last iteration pos which means extend end to the right {end - original_end}bp")
            copies_alignment_end = end
            basepair_accuracy = 0
            break
        else:
            print("new interval", new_interval, interval_size(new_interval))
            print(q_end_count)
            min_q_end_percentage = min(q_end_list) / interval_size(new_interval)
            print(min(q_end_list), min_q_end_percentage * 100, "%")
            print("")

        # Update end for the next iteration
        end = new_end
    extended_2_right_interval = create_interval(chr, pos, copies_alignment_end)
    return extended_2_right_interval, basepair_accuracy


# This function is part of the basepair correction function of T2T-ACE
# It will extend the DUP interval to the left until no more copies are found
def extend_2_left(dup_interval, calling_reference_fasta, called_ref_aligner,
                                 truth_ref_aligner):
    chr, pos, end = parse_interval(dup_interval)
    if chr == "chrX" or chr == "chrY":
        print("Basepair correction Method doesn't support DUPs on chrX or chrY")
        return dup_interval, 0
    original_pos = pos
    # The window size is 10% of the interval size
    window_size = int(np.round(interval_size(dup_interval) * 0.1, 0))
    print("extending window size", window_size)

    q_start_list = [q_st for interval, strand, q_st, q_en in
                    align_interval(dup_interval, calling_reference_fasta, called_ref_aligner,
                                   truth_ref_aligner)[1]]
    original_q_start_count = len(q_start_list)
    new_interval = dup_interval
    copies_alignment_start = ""
    basepair_accuracy = 0
    while min(q_start_list) <= 1:
        new_pos = pos - window_size
        new_interval = create_interval(chr, new_pos, end)
        new_alignments = align_interval(new_interval, calling_reference_fasta, called_ref_aligner,
                                        truth_ref_aligner)

        q_start_list = sorted([q_st for interval, strand, q_st, q_en in new_alignments[1]])
        q_start_count = [[q_start, q_start_list.count(q_start)] for q_start in sorted(set(q_start_list))]
        # If only the two matches are at the beginning of the alignment, then the interval pos needs to be moved to the right
        if q_start_count[0][1] == 2 and q_start_count[-1][0] != 1:
            print("new interval", new_interval, interval_size(new_interval))
            print(q_start_count)
            if len(q_start_count) > 1:
                print(f"extend pos to the left by {original_pos-(new_pos+q_start_count[-1][0])}bp")
                copies_alignment_start = new_pos + q_start_count[-1][0]
                basepair_accuracy = 1
            elif len(q_start_count) == 1:
                print(f"the window size is too big")
                basepair_accuracy = 0
                if pos < original_pos:
                    copies_alignment_start = pos
                    print(f"Keep the last iteration pos which means extend pos to the left {original_pos - pos}bp")
                else:
                    copies_alignment_start = original_pos
                    print(f"Keep the original pos")
            break
        # If there is only one non-zero(or non-one) start position for all alignments
        # Then this sequence might hit the unresolved region in hg38
        elif min(q_start_list)>1 and len(q_start_count) >= 1:
            print("new interval", new_interval, interval_size(new_interval))
            print(q_start_count)
            print("This sequence might hit the unresolved region in hg38")
            copies_alignment_start = new_pos + q_start_count[-1][0]
            print(f"extend pos to the left by {original_pos - (new_pos + q_start_count[0][0])}bp")
            basepair_accuracy = 0
            break
        elif len(q_start_count) > 1 and q_start_count[-1][1] ==2 and q_start_count[-1][0] > 3:
            print("new interval", new_interval, interval_size(new_interval))
            print(q_start_count)
            print("The sequence isn't homozygous")
            copies_alignment_start = new_pos + q_start_count[-1][0]
            print(f"extend pos to the left by {original_pos - (new_pos + q_start_count[-1][0])}bp")
            basepair_accuracy = 0
            break
        elif len(q_start_list) > original_q_start_count:
            print("new interval", new_interval, interval_size(new_interval))
            print(q_start_count)
            print("The sequence has extended to copies that doesn't exist in the original interval alignment")
            copies_alignment_start = pos
            print(f"Keep the last iteration pos which means extend pos to the left {original_pos - pos}bp")
            basepair_accuracy = 0
            break
        else:
            print("new interval", new_interval, interval_size(new_interval))
            print(q_start_count)
            max_q_start_percentage = ((end - new_pos) + 1) / interval_size(new_interval)
            print(max(q_start_list), max_q_start_percentage * 100, "%")
            print("")

        # Update end for the next iteration
        pos = new_pos
    extended_2_left_interval = create_interval(chr, copies_alignment_start, end)
    return extended_2_left_interval, basepair_accuracy