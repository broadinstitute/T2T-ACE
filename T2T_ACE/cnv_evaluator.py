import mappy as mp
import T2T_ACE.cnv_alignments as cnv_alignments

import T2T_ACE.alignment_utilities as au
from T2T_ACE.genomic_queries import get_sequence_from_interval, get_flanking_regions


class CNVEvaluator:
    calling_reference: str
    maternal_reference: str
    paternal_reference: str

    calling_aligner: mp.Aligner
    maternal_aligner: mp.Aligner
    paternal_aligner: mp.Aligner

    def __init__(self, calling_reference: str, maternal_reference: str, paternal_reference: str):
        """
        Initializes a CNVEvaluation instance.

        Parameters:
        calling_reference (str): Path to the reference genome used for calling.
        maternal_reference (str): Path to the maternal reference genome.
        paternal_reference (str): Path to the paternal reference genome.
        """
        self.calling_reference = calling_reference
        self.maternal_reference = maternal_reference
        self.paternal_reference = paternal_reference

        self.calling_aligner = au.load_reference(calling_reference)
        self.maternal_aligner = au.load_reference(maternal_reference)
        self.paternal_aligner = au.load_reference(paternal_reference)

    def __str__(self):
        """
        String representation of the CNVEvaluation.
        """
        return (f"Calling Reference: {self.calling_reference}, "
                f"Maternal Reference: {self.maternal_reference}, "
                f"Paternal Reference: {self.paternal_reference}")

    def align_cnv_for_evaluation(self, interval: str, cnv_type: str, genotype: str) -> cnv_alignments.CNVAlignments:
        """
        Evaluates a CNV region.

        Parameters:
        interval (str): Interval to evaluate.
        cnv_type (str): Type of the CNV (e.g., "GAIN", "LOSS").
        genotype (str): Genotype of the CNV.
        """

        raw_sequence = get_sequence_from_interval(self.calling_reference, interval)
        left_flank, right_flank = get_flanking_regions(self.calling_reference, interval, padding=500)

        calling_reference_hits = [_ for _ in self.calling_aligner.map(raw_sequence)]
        maternal_hits = [_ for _ in self.maternal_aligner.map(raw_sequence)]
        paternal_hits = [_ for _ in self.paternal_aligner.map(raw_sequence)]

        left_flank_calling_reference_hits = [_ for _ in self.calling_aligner.map(left_flank)]
        left_flank_maternal_hits = [_ for _ in self.maternal_aligner.map(left_flank)]
        left_flank_paternal_hits = [_ for _ in self.paternal_aligner.map(left_flank)]

        right_flank_calling_reference_hits = [_ for _ in self.calling_aligner.map(right_flank)]
        right_flank_maternal_hits = [_ for _ in self.maternal_aligner.map(right_flank)]
        right_flank_paternal_hits = [_ for _ in self.paternal_aligner.map(right_flank)]

        alignments = cnv_alignments.CNVAlignments(interval,
                                                  cnv_type,
                                                  genotype,
                                                  calling_reference_hits,
                                                  maternal_hits,
                                                  paternal_hits,
                                                  left_flank_calling_reference_hits,
                                                  left_flank_maternal_hits,
                                                  left_flank_paternal_hits,
                                                  right_flank_calling_reference_hits,
                                                  right_flank_maternal_hits,
                                                  right_flank_paternal_hits)

        return alignments

# %%
