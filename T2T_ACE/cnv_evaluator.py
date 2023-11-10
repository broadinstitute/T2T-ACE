import mappy as mp
import T2T_ACE.alignment_utilities as au


class CNV_Evaluator:

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
        return f"CNV Location: {self.location}, Type: {self.cnv_type}, Genotype: {self.genotype}"

    def evaluate_deletion(self, cnv_interval: str, gain_or_loss: str, called_reference: str, eval_assembly) -> list:
        """
        Method to evaluate if the CNV is potentially deleterious.
        This can be expanded based on specific criteria or external data integration.
        """
        # Placeholder logic; real implementation requires more detailed criteria
        return self.cnv_type == 'loss'


