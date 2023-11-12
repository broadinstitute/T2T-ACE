import mappy as mp
from typing import List


class CNVEvaluation:

    location: str
    cnv_type: str
    genotype: str

    calling_reference_hits: List[mp.Alignment]
    maternal_hits: List[mp.Alignment]
    paternal_hits: List[mp.Alignment]
    left_flank_calling_reference_hits: List[mp.Alignment]
    left_flank_maternal_hits: List[mp.Alignment]
    left_flank_paternal_hits: List[mp.Alignment]
    right_flank_calling_reference_hits: List[mp.Alignment]
    right_flank_maternal_hits: List[mp.Alignment]
    right_flank_paternal_hits: List[mp.Alignment]

    def __init__(self, location: str, cnv_type: str, genotype: str,
                 calling_reference_hits: List[mp.Alignment],
                 maternal_hits: List[mp.Alignment],
                 paternal_hits: List[mp.Alignment],
                 left_flank_calling_reference_hits: List[mp.Alignment],
                 left_flank_maternal_hits: List[mp.Alignment],
                 left_flank_paternal_hits: List[mp.Alignment],
                 right_flank_calling_reference_hits: List[mp.Alignment],
                 right_flank_maternal_hits: List[mp.Alignment],
                 right_flank_paternal_hits: List[mp.Alignment]):
        """
        Initializes a CNVEvaluation instance.

        Parameters:
        location (str): Location of the CNV.
        cnv_type (str): Type of the CNV.
        genotype (str): Genotype of the CNV.
        """
        self.location = location
        self.cnv_type = cnv_type
        self.genotype = genotype

        self.calling_reference_hits = calling_reference_hits
        self.maternal_hits = maternal_hits
        self.paternal_hits = paternal_hits
        self.left_flank_calling_reference_hits = left_flank_calling_reference_hits
        self.left_flank_maternal_hits = left_flank_maternal_hits
        self.left_flank_paternal_hits = left_flank_paternal_hits
        self.right_flank_calling_reference_hits = right_flank_calling_reference_hits
        self.right_flank_maternal_hits = right_flank_maternal_hits
        self.right_flank_paternal_hits = right_flank_paternal_hits

    def __str__(self):
        """
        String representation of the CNVEvaluation.
        """
        return f"CNV Location: {self.location}, Type: {self.cnv_type}, Genotype: {self.genotype}"
