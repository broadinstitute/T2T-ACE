# T2T-ACE
 *A*ccurate *C*NV *E*valuation Using Telomere-to-Telomere Assemblies
![T2T-ACE-logo](docs/logo.png)


## About T2T-ACE

T2T-ACE goes beyond conventional CNV assessment methods. It harnesses telomere-to-telomere assemblies to evaluate CNVs with unmatched precision. This innovative tool utilizes alignment-based techniques, enabling researchers to:

- Determine the correctness of CNV events.
- Identify the precise locations of CNVs.
- Genotype CNV events of interest.

## Key Features

- **Unparalleled Accuracy:** T2T-ACE uses comprehensive telomere-to-telomere assemblies to ensure reliable CNV evaluation.

- **Precision Matters:** T2T-ACE can suggest the exact locations of CNV events, providing critical insights 
that can inform mathods or filteting strategies of your CNV caller

- **Genotyping Capability:** Not only can you validate CNV events, but you can also determine their genotypes.

## Design Description
### DEL evaluation
T2T-ACE align the left and right flanking regions of a DEL variant called in reference genome (hg38) to the HG002-T2T reference.
By calculating the distance between the left and right flanking regions are aligned in HG002-T2T reference, we can determine the correctness and genotype of this DEL variant.

- **Correctness:** If the distance between the left and right flanking regions are aligned in HG002-T2T reference is within the range of 0.8 * (length of the DEL variant) and 1.2 * (length of the DEL variant), we consider this DEL variant is FP.
- **Genotype:** 
![DEL](docs/DEL_eval_logic.png)
- **Het DEL Example:** ![DEL](docs/Het_DEL_example.png)
- **Hom DEL Example:** ![DEL](docs/Hom_DEL_example.png)

### DUP evaluation
T2T-ACE align the DUP variant called in reference genome (hg38) to the HG002-T2T reference. T2T-ACE aligns the DNA sequence 
representing a DUP variant called in hg38 to the HG002-T2T reference.  

- **Correctness:** 
A TP DUP event is characterized by a higher copy number 
in HG002-T2T than in hg38, while a FP DUP event shows fewer copies in HG002-T2T than hg38. Since the hg38 assembly is haploid, 
and the HG002-T2T assembly is diploid, we assume a copy neutral event has one copy in hg38 and two copies (one for the maternal 
and paternal haplotypes) in HG002-T2T.  A copy is defined as an alignment made by minimap2 where the alignment length is at least 50%
of the query length. This method allows us to identify not only whether a CNV event is correct, but in the case of duplications 
also allows us to identify the precise number, locations, and genotype of duplication events, even if they occur on 
separate chromosomes from the original call on hg38.
- **Hom DUP Example:** ![DUP](docs/Hom_DUP_example.png)
- **Het DUP Example:** ![DUP](docs/Het_DUP_example.png)
- **Basepair-level Correction:** ![DUP](docs/dup_interval_correction_example.png)

