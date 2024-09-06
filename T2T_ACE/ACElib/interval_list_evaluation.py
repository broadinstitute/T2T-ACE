import pandas as pd
import gzip
from .duplication_evaluation import eval_dup_interval
from . import validator as v

def read_vcf(vcf_path):
    # Check if the VCF file is gzipped
    open_func = gzip.open if vcf_path.endswith('.gz') else open
    with open_func(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith("#CHROM"):
                vcf_header_names = line.strip().split('\t')
                break
    new_vcf_header_names = [i.replace("#", "") for i in vcf_header_names]
    df = pd.read_csv(vcf_path, comment='#', sep='\t', header=None, names=new_vcf_header_names)
    return df

class eval_interval_list:
    def __init__(self, interval_list, calling_reference_fasta, truth_reference_fasta, called_ref_aligner, truth_ref_aligner):
        self.interval_list = interval_list
        self.calling_reference_fasta = calling_reference_fasta
        self.called_ref_aligner = called_ref_aligner
        self.truth_reference_fasta = truth_reference_fasta
        self.truth_ref_aligner = truth_ref_aligner

    # This function will create a summary df for all the input DUP intervals
    def create_dup_sum(self):
        hg2_dup_summary_df = pd.DataFrame()
        hg2_dup_summary_df['interval'] = self.interval_list

        # Collect the results for each interval from dup_eval
        dup_interval_size_list = []
        original_hg38_hit_count_list = []
        original_hg2_hit_count_list = []
        original_hg2_mat_hit_count_list = []
        original_hg2_pat_hit_count_list = []
        original_dup_interval_major_classification_list = []
        original_dup_interval_sub_classification_list = []
        original_dup_interval_contain_big_gap_list = []

        for interval in self.interval_list:
            print(interval)
            # If the interval returns True, then it is a valid interval and we can collect the results
            # If the interval returns False, then it is not a valid interval and we will fill the results with 'NA'
            if eval_dup_interval(interval, self.calling_reference_fasta,
                                                                         self.truth_reference_fasta, self.called_ref_aligner,
                                                                         self.truth_ref_aligner).analyzeDuplicationIntervals():
                interval_sum_dict = eval_dup_interval(interval, self.calling_reference_fasta,
                                                                             self.truth_reference_fasta, self.called_ref_aligner,
                                                                             self.truth_ref_aligner).analyzeDuplicationIntervals()
                dup_interval_size_list.append(interval_sum_dict['dup_interval_size'])
                original_hg38_hit_count_list.append(interval_sum_dict['original_hg38_hit_count'])
                original_hg2_hit_count_list.append(interval_sum_dict['original_hg2_hit_count'])
                original_hg2_mat_hit_count_list.append(interval_sum_dict['original_hg2_mat_hit_count'])
                original_hg2_pat_hit_count_list.append(interval_sum_dict['original_hg2_pat_hit_count'])
                original_dup_interval_major_classification_list.append(
                    interval_sum_dict['original_dup_interval_major_classification'])
                original_dup_interval_sub_classification_list.append(
                    interval_sum_dict['original_dup_interval_sub_classification'])
                original_dup_interval_contain_big_gap_list.append(
                    interval_sum_dict['original_dup_interval_contain_big_gap'])
            else:
                dup_interval_size_list.append('NA')
                original_hg38_hit_count_list.append('NA')
                original_hg2_hit_count_list.append('NA')
                original_hg2_mat_hit_count_list.append('NA')
                original_hg2_pat_hit_count_list.append('NA')
                original_dup_interval_major_classification_list.append('NA')
                original_dup_interval_sub_classification_list.append('NA')
                original_dup_interval_contain_big_gap_list.append('NA')

        hg2_dup_summary_df['dup_interval_size'] = dup_interval_size_list
        hg2_dup_summary_df['original_hg38_hit_count'] = original_hg38_hit_count_list
        hg2_dup_summary_df['original_hg2_hit_count'] = original_hg2_hit_count_list
        hg2_dup_summary_df['original_hg2_mat_hit_count'] = original_hg2_mat_hit_count_list
        hg2_dup_summary_df['original_hg2_pat_hit_count'] = original_hg2_pat_hit_count_list
        hg2_dup_summary_df[
            'original_dup_interval_major_classification'] = original_dup_interval_major_classification_list
        hg2_dup_summary_df['original_dup_interval_sub_classification'] = original_dup_interval_sub_classification_list
        hg2_dup_summary_df['original_dup_interval_contain_big_gap'] = original_dup_interval_contain_big_gap_list
        return hg2_dup_summary_df

    # Classify all the DELs by DRAGEN with collect_del_flankings function in validator.py
    def classify_list_of_DELs(self):
        del_sum_df = pd.DataFrame()
        del_interval_list = []
        del_interval_size_list = []
        flanking_size_list = []
        left_flanking_interval_list = []
        right_flanking_interval_list = []
        left_flanking_hg38_hits_list = []
        right_flanking_hg38_hits_list = []
        left_flanking_hg2_hits_list = []
        right_flanking_hg2_hits_list = []
        distance_between_flankings_list = []
        flanking_connection_strand_list = []
        major_classification_list = []
        minor_classification_list = []
        for interval in self.interval_list:
            print(interval)
            sum_dict = v.collect_del_flankings(interval, self.calling_reference_fasta, self.called_ref_aligner, self.truth_ref_aligner)
            del_interval_list.append(sum_dict['del_interval'])
            del_interval_size_list.append(sum_dict['del_interval_size'])
            flanking_size_list.append(sum_dict['flanking_size'])
            left_flanking_interval_list.append(sum_dict['left_flanking_interval'])
            right_flanking_interval_list.append(sum_dict['right_flanking_interval'])
            left_flanking_hg38_hits_list.append(sum_dict['left_flanking_hg38_hits'])
            right_flanking_hg38_hits_list.append(sum_dict['right_flanking_hg38_hits'])
            left_flanking_hg2_hits_list.append(sum_dict['left_flanking_hg2_hits'])
            right_flanking_hg2_hits_list.append(sum_dict['right_flanking_hg2_hits'])
            distance_between_flankings_list.append(sum_dict['distance_between_flankings'])
            flanking_connection_strand_list.append(sum_dict['flanking_connection_strand'])
            major_classification_list.append(sum_dict['classification'])
            minor_classification_list.append(sum_dict['minor_classification'])
        del_sum_df['del_interval'] = del_interval_list
        del_sum_df['del_interval_size'] = del_interval_size_list
        del_sum_df['flanking_size'] = flanking_size_list
        del_sum_df['left_flanking_interval'] = left_flanking_interval_list
        del_sum_df['right_flanking_interval'] = right_flanking_interval_list
        del_sum_df['left_flanking_hg38_hits'] = left_flanking_hg38_hits_list
        del_sum_df['right_flanking_hg38_hits'] = right_flanking_hg38_hits_list
        del_sum_df['left_flanking_hg2_hits'] = left_flanking_hg2_hits_list
        del_sum_df['right_flanking_hg2_hits'] = right_flanking_hg2_hits_list
        del_sum_df['distance_between_flankings'] = distance_between_flankings_list
        del_sum_df['flanking_connection_strand'] = flanking_connection_strand_list
        del_sum_df['major_classification'] = major_classification_list
        del_sum_df['minor_classification'] = minor_classification_list
        return del_sum_df




