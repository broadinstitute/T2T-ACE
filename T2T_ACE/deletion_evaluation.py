import numpy as np
from T2T_ACE.validator import collect_del_flankings, interval_size

class eval_del_interval:
    def __init__(self, del_interval, calling_reference_fasta, called_ref_aligner, truth_ref_aligner):
        self.del_interval = del_interval
        self.calling_reference_fasta = calling_reference_fasta
        self.called_ref_aligner = called_ref_aligner
        self.truth_ref_aligner = truth_ref_aligner





