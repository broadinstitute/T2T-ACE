# T2T-ACE
 *A*ccurate *C*NV *E*valuation Using Telomere-to-Telomere Assemblies
 
## Run T2T-ACE
This tool is designed to evaluate the accuracy of CNV calls using the T2T assembly as a reference. 
The tool will align the CNV calls to the T2T assembly and the hg38 assembly and compare the alignment results. 
```
python3 T2T_ACE/run_T2T-ACE.py --cnv_vcf <cnv_vcf> --t2t_ref <t2t_assembly.fa> --hg38_ref <hg38_assembly.fa>
```

## Download Assembly Files
The T2T assembly and the hg38 assembly can be downloaded from the following links:
* [GENCODE](https://www.gencodegenes.org/human/) hg38 primary assembly: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz
* [HG002](https://github.com/marbl/HG002) T2T assembly: https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.0.1.fasta.gz