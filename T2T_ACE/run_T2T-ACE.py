import ACElib.interval_list_evaluation as IntervalListEvaluation
from ACElib.interval_list_evaluation import read_vcf
from ACElib.alignment_utilities import load_reference
import argparse


argparser = argparse.ArgumentParser(description="T2T-ACE: A tool for evaluating CNV calls for genome with T2T assembly.")
argparser.add_argument('-v', "--cnv_vcf_path", dest="cnv_vcf_path", help="Path to CNV VCF file in hg38.")
argparser.add_argument('-t2t', "--t2t_ref", dest="t2t_ref", help="Path to T2T reference genome.")
argparser.add_argument('-hg38', "--hg38_ref", dest="hg38_ref", help="Path to hg38 reference genome.")
args = argparser.parse_args()

vcf_path = args.cnv_vcf_path
hg002t2t_ref_path = args.t2t_ref
hg38_ref_path = args.hg38_ref

print("CNV VCF Path: ", vcf_path)
print("T2T Reference Path: ", hg002t2t_ref_path)
print("HG38 Reference Path: ", hg38_ref_path)

# Load the minimap2 aligner from reference fasta file
# Load HG002 T2T reference
hg002t2t = load_reference(hg002t2t_ref_path)
# Load hg38 reference
hg38 = load_reference(hg38_ref_path)

# Get the CNV calls from the input VCF file
HG2_cnv_vcf_df = read_vcf(vcf_path)
pass_cnv_list = HG2_cnv_vcf_df[HG2_cnv_vcf_df['FILTER']=='PASS']

# Split the CNV calls into DEL and DUP
del_pass_intervals = []
dup_pass_intervals = []
for index, row in HG2_cnv_vcf_df.iterrows():
    interval = row['CHROM'] + ':' + str(row['POS']) + '-' + str(row['INFO'].split('END=')[1].split(';')[0])
    if row['ALT'] == '<DEL>':
        del_pass_intervals.append(interval)
    elif row['ALT'] == '<DUP>':
        dup_pass_intervals.append(interval)
print('Number of DEL intervals within the Input VCF:', len(del_pass_intervals))
print('Number of DEL intervals within the Input VCF:', len(dup_pass_intervals))

print("Evaluating the CNV calls...")
dup_eval_sum_df = IntervalListEvaluation.eval_interval_list(dup_pass_intervals[:5], hg38_ref_path, hg002t2t_ref_path, hg38, hg002t2t).create_dup_sum()
del_eval_sum_df = IntervalListEvaluation.eval_interval_list(del_pass_intervals[:5], hg38_ref_path, hg002t2t_ref_path, hg38, hg002t2t).classify_list_of_DELs()


dup_eval_sum_df.to_csv('output_dup_eval_sum.csv', index=False)
del_eval_sum_df.to_csv('output_del_eval_sum.csv', index=False)
