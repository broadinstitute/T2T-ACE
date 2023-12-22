import numpy as np
from T2T_ACE.validator import collect_del_flankings, interval_size

class eval_del_interval:
    def __init__(self, del_interval, calling_reference_fasta, called_ref_aligner, truth_ref_aligner):
        self.del_interval = del_interval
        self.calling_reference_fasta = calling_reference_fasta
        self.called_ref_aligner = called_ref_aligner
        self.truth_ref_aligner = truth_ref_aligner
    def analyzeDeletionIntervals(self):
        del_sum_dict = collect_del_flankings(self.del_interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)
        del_interval_size = del_sum_dict['del_interval_size']
        left_flanking_hg38_hits = del_sum_dict['left_flanking_hg38_hits']
        right_flanking_hg38_hits = del_sum_dict['right_flanking_hg38_hits']
        left_flanking_hg2_hits = del_sum_dict['left_flanking_hg2_hits']
        right_flanking_hg2_hits = del_sum_dict['right_flanking_hg2_hits']
        distance_between_flankings = del_sum_dict['distance_between_flankings']
        if left_flanking_hg38_hits > 1 or right_flanking_hg38_hits >1:
            if distance_between_flankings and min(distance_between_flankings) <= del_interval_size*0.5:
                classifcation = 'DEL in DUP'
            else:
                classifcation = 'False DEL'
        elif left_flanking_hg38_hits == 1 and right_flanking_hg38_hits == 1:
            if distance_between_flankings:
                del_flanking_dist = [distance for distance in distance_between_flankings if distance <= del_interval_size*0.5]
                if len(distance_between_flankings) == len(del_flanking_dist):
                    classifcation = 'Homozygous DEL'
                elif len(distance_between_flankings) > len(del_flanking_dist) and len(del_flanking_dist) > 0:
                    classifcation = 'Heterozygous DEL'
                else:
                    classifcation = 'False DEL'
        else:
            classifcation = 'Unclassified DEL'
        return classifcation




