import ACElib.interval_list_evaluation as IntervalListEvaluation
from ACElib.interval_list_evaluation import read_vcf
from ACElib.alignment_utilities import load_reference
import argparse


argparser = argparse.ArgumentParser(description="T2T-ACE: A tool for evaluating CNV calls for genome with T2T assembly.")
argparser.add_argument('-v', "--cnv_vcf", dest="cnv_vcf_path", help="Path to CNV VCF file in hg38.", required=False)
argparser.add_argument('-t2t', "--t2t_ref", dest="t2t_ref", help="Path to T2T reference genome.", required=True)
argparser.add_argument('-hg38', "--hg38_ref", dest="hg38_ref", help="Path to hg38 reference genome.", required=True)
argparser.add_argument('-test', "--test", dest="test", help="True, if you want to just run a test on first 5 CNV calls.", type=bool, required=False)
argparser.add_argument('-del', '--del_txt', dest='del_txt_path', help='Path to the txt file containing the list of DEL intervals.', required=False)
argparser.add_argument('-dup', '--dup_txt', dest='dup_txt_path', help='Path to the txt file containing the list of DUP intervals.', required=False)
argparser.add_argument('-o', '--output', dest='output_path', help='Path to the output directory.', required=False)

args = argparser.parse_args()

vcf_path = args.cnv_vcf_path
hg002t2t_ref_path = args.t2t_ref
hg38_ref_path = args.hg38_ref
del_txt_path = args.del_txt_path
dup_txt_path = args.dup_txt_path
output_path = args.output_path

print("CNV VCF Path: ", vcf_path)
print("T2T Reference Path: ", hg002t2t_ref_path)
print("HG38 Reference Path: ", hg38_ref_path)
print("Test: ", args.test)
print("DEL Input Path: ", del_txt_path)
print("DUP Input Path: ", dup_txt_path)
print("Output Path: ", output_path)

# Get CNV DEL and DUP intervals
del_pass_intervals = []
dup_pass_intervals = []
# Get the CNV calls from the input VCF file
if vcf_path:
    HG2_cnv_vcf_df = read_vcf(vcf_path)
    print('Number of CNV calls in the VCF:', HG2_cnv_vcf_df.shape[0])
    passing_cnv_calls = HG2_cnv_vcf_df[HG2_cnv_vcf_df['FILTER'] == 'PASS']
    print('Number of CNV calls passing the filter:', passing_cnv_calls.shape[0])

    for index, row in HG2_cnv_vcf_df.iterrows():
        interval = row['CHROM'] + ':' + str(row['POS']) + '-' + str(row['INFO'].split('END=')[1].split(';')[0])
        if row['ALT'] == '<DEL>' and row['FILTER'] == 'PASS':
            del_pass_intervals.append(interval)
        elif row['ALT'] == '<DUP>' and row['FILTER'] == 'PASS':
            dup_pass_intervals.append(interval)
    print('Number of DEL intervals within the Input VCF:', len(del_pass_intervals))
    print('Number of DEL intervals within the Input VCF:', len(dup_pass_intervals))
elif del_txt_path:
    with open(del_txt_path) as f:
        for line in f:
            del_pass_intervals.append(line.strip())
    print('Number of DEL intervals within the Input TXT:', len(del_pass_intervals))
elif dup_txt_path:
    with open(dup_txt_path) as f:
        for line in f:
            dup_pass_intervals.append(line.strip())
    print('Number of DUP intervals within the Input TXT:', len(dup_pass_intervals))
else:
    raise ValueError("Please provide either the VCF file or the txt file containing the list of DEL and DUP intervals.")

print('Loading the reference genomes...')
# Load the minimap2 aligner from reference fasta file
# Load HG002 T2T reference
hg002t2t = load_reference(hg002t2t_ref_path)
# Load hg38 reference
hg38 = load_reference(hg38_ref_path)

print("Evaluating the CNV calls...")
if dup_pass_intervals:
    if args.test:
        del_pass_intervals = del_pass_intervals[:5]
    dup_eval_sum_df = IntervalListEvaluation.eval_interval_list(dup_pass_intervals, hg38_ref_path, hg002t2t_ref_path, hg38, hg002t2t).create_dup_sum()
    if output_path:
        dup_eval_sum_df.to_csv(output_path + 'DUP_eval_sum.csv', index=False)
    else:
        dup_eval_sum_df.to_csv('output_DUP_eval_sum.csv', index=False)
if del_pass_intervals:
    if args.test:
        del_pass_intervals = del_pass_intervals[:5]
    del_eval_sum_df = IntervalListEvaluation.eval_interval_list(del_pass_intervals, hg38_ref_path, hg002t2t_ref_path, hg38, hg002t2t).classify_list_of_DELs()
    if output_path:
        del_eval_sum_df.to_csv(output_path + 'DEL_eval_sum.csv', index=False)
    else:
        del_eval_sum_df.to_csv('output_DEL_eval_sum.csv', index=False)

# #debugging
# del_eval_sum_df = IntervalListEvaluation.eval_interval_list(["chr2:35751954-35767502"], hg38_ref_path, hg002t2t_ref_path, hg38,
#                                                             hg002t2t).classify_list_of_DELs()
# print(del_eval_sum_df)

