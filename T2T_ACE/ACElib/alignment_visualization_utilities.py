import os
import pandas as pd
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict


# Import chrom size data for hg38 and HG002
current_dir = os.path.dirname(os.path.abspath(__file__))
hg38_chrom_size_filepath = os.path.join(current_dir, 'resources/hg38_chrom_size.txt')
hg002_mat_dict_filepath = os.path.join(current_dir, "resources/hg002v0.9.mat_Y_EBV_MT.dict")
hg002_pat_dict_filepath = os.path.join(current_dir, "resources/hg002v0.9.pat_X_EBV_MT.dict")

# Process/Extract chromo size data
hg38_chrom_size_df = pd.read_csv(hg38_chrom_size_filepath, sep='\t', names=['chr','length'])
hg38_chromosome_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
hg38_chrom_size_df.sort_values('chr', key=lambda column: column.map(lambda e: hg38_chromosome_order.index(e)),
                               inplace=True, ignore_index=True)

mat_dict_df = pd.read_csv(hg002_mat_dict_filepath,sep='\t', skiprows=[0], names=['sn_chr','ln_length','m5','source'])
mat_dict_df = mat_dict_df[mat_dict_df['sn_chr']!='SN:chrEBV'][['sn_chr','ln_length']]

mat_chr_size_df = pd.DataFrame()
mat_chr_size_df['chr'] = [i.split(':')[1] for i in mat_dict_df['sn_chr'].tolist()]
mat_chr_size_df['length'] = [int(i.split(':')[1]) for i in mat_dict_df['ln_length'].tolist()]
mat_chromosome_order = [f'chr{i}_MATERNAL' for i in range(1, 23)] + ['chrX_MATERNAL', 'chrY_PATERNAL']
mat_chr_size_df.sort_values('chr', key=lambda column: column.map(lambda e: mat_chromosome_order.index(e)),
                            inplace=True, ignore_index=True)

pat_dict_df = pd.read_csv(hg002_pat_dict_filepath, sep='\t', skiprows=[0], names=['sn_chr','ln_length','m5','source'])
pat_dict_df = pat_dict_df[(pat_dict_df['sn_chr']!='SN:chrM')&(pat_dict_df['sn_chr']!='SN:chrEBV')][['sn_chr','ln_length']]

pat_chr_size_df = pd.DataFrame()
pat_chr_size_df['chr'] = [i.split(':')[1] for i in pat_dict_df['sn_chr'].tolist()]
pat_chr_size_df['length'] = [int(i.split(':')[1]) for i in pat_dict_df['ln_length'].tolist()]
pat_chromosome_order = [f'chr{i}_PATERNAL' for i in range(1, 23)] + ['chrX_MATERNAL', 'chrY_PATERNAL']
pat_chr_size_df.sort_values('chr', key=lambda column: column.map(lambda e: pat_chromosome_order.index(e)),
                            inplace=True, ignore_index=True)

# Set seaborn style
sns.set(style='white')
class PlotIntervals:
    def __init__(self, hg38_interval_list, hg002_interval_list):
        # plot interval class assumes the first hg38 interval is the event interval
        self.hg38_interval_list = hg38_interval_list
        self.hg002_interval_list = hg002_interval_list

    def chromosome_karyotype_plot(self):
        hg002_mat_interval_list = [i for i in self.hg002_interval_list if 'MATERNAL' in i]
        hg002_pat_interval_list = [i for i in self.hg002_interval_list if 'PATERNAL' in i]
        sns.set_style("white")
        # remove top and right axes spines
        fig, axs = plt.subplots(1, 3, figsize=(21, 7))

        sns.barplot(data=hg38_chrom_size_df, x="length", y="chr", linewidth=0.8, edgecolor="0.1", facecolor=(0, 0, 0, 0), ax=axs[0])
        sns.despine()
        axs[0].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, n: int(x/1e6)))
        axs[0].set_xlabel('Chromosome Position (Mbp)', fontsize=14)
        axs[0].set_ylabel('', fontsize=14)
        axs[0].set_title("hg38", fontsize=14)

        query_chr = self.hg38_interval_list[0].split(':')[0]
        query_pos = int(self.hg38_interval_list[0].split(':')[1].split('-')[0])
        query_end = int(self.hg38_interval_list[0].split(':')[1].split('-')[1])
        if query_chr in hg38_chrom_size_df['chr'].values:
            # Get the y position of chr
            y_pos = hg38_chrom_size_df.index[hg38_chrom_size_df['chr'] == query_chr][0]
            # Calculate the x positions for the interval
            max_length = hg38_chrom_size_df['length'].max()
            x_start = query_pos / max_length
            x_end = query_end / max_length
            axs[0].axhspan(y_pos-0.35, y_pos+0.35, xmin=x_start, xmax=x_end, color='red', alpha=1, linewidth=3)
            axs[0].set_xlim([0,max_length])

        if len(self.hg38_interval_list)>1:
            for interval in self.hg38_interval_list[1:]:
                chr = interval.split(':')[0]
                pos = int(interval.split(':')[1].split('-')[0])
                end = int(interval.split(':')[1].split('-')[1])

                if chr in hg38_chrom_size_df['chr'].values:
                    # Get the y position of chr
                    y_pos = hg38_chrom_size_df.index[hg38_chrom_size_df['chr'] == chr][0]
                    # Calculate the x positions for the interval
                    max_length = hg38_chrom_size_df['length'].max()
                    x_start = pos / max_length
                    x_end = end / max_length
                    axs[0].axhspan(y_pos-0.35, y_pos+0.35, xmin=x_start, xmax=x_end, facecolor='black', alpha=1,
                                   linewidth=3, edgecolor=".1")
        # Add number to the end of chromosome
        for chr in set([i.split(':')[0] for i in self.hg38_interval_list]):
            y_pos = hg38_chrom_size_df.index[hg38_chrom_size_df['chr'] == chr][0]
            axs[0].text(hg38_chrom_size_df[hg38_chrom_size_df['chr'] == chr]['length']+5000000,
                        y_pos+0.25, [i.split(':')[0] for i in self.hg38_interval_list].count(chr))

        sns.barplot(data=mat_chr_size_df, x="length", y="chr", linewidth=1, edgecolor=".5",
                    facecolor=(0, 0, 0, 0), ax=axs[1])
        axs[1].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, n: int(x/1e6)))
        axs[1].set_xlabel('Chromosome Position (Mbp)', fontsize=14)
        axs[1].set_ylabel('', fontsize=14)
        axs[1].set_title("HG002 Maternal", fontsize=14)

        for interval in hg002_mat_interval_list:
            chr = interval.split(':')[0]
            pos = int(interval.split(':')[1].split('-')[0])
            end = int(interval.split(':')[1].split('-')[1])

            if chr in mat_chr_size_df['chr'].values:
                # Get the y position of chr
                y_pos = mat_chr_size_df.index[mat_chr_size_df['chr'] == chr][0]
                # Calculate the x positions for the interval
                max_length = mat_chr_size_df['length'].max()
                x_start = pos / max_length
                x_end = end / max_length
                axs[1].axhspan(y_pos-0.35, y_pos+0.35, xmin=x_start, xmax=x_end, facecolor='black', alpha=1,
                               linewidth=3,edgecolor=".1")
                axs[1].set_xlim([0,max_length])
        # Add number to the end of chromosome
        for chr in set([i.split(':')[0] for i in hg002_mat_interval_list]):
            y_pos = mat_chr_size_df.index[mat_chr_size_df['chr'] == chr][0]
            axs[1].text(mat_chr_size_df[mat_chr_size_df['chr'] == chr]['length']+5000000, y_pos+0.25 ,
                        [i.split(':')[0] for i in hg002_mat_interval_list].count(chr))

        sns.barplot(data=pat_chr_size_df, x="length", y="chr", linewidth=1, edgecolor=".5", facecolor=(0, 0, 0, 0),
                    ax=axs[2])
        axs[2].xaxis.set_major_formatter(plt.FuncFormatter(lambda x, n: int(x/1e6)))
        axs[2].set_xlabel('Chromosome Position (Mbp)', fontsize=14)
        axs[2].set_ylabel('', fontsize=14)
        axs[2].set_title("HG002 Paternal", fontsize=14)

        for interval in hg002_pat_interval_list:
            chr = interval.split(':')[0]
            pos = int(interval.split(':')[1].split('-')[0])
            end = int(interval.split(':')[1].split('-')[1])

            if chr in pat_chr_size_df['chr'].values:
                # Get the y position of chr
                y_pos = pat_chr_size_df.index[pat_chr_size_df['chr'] == chr][0]
                # Calculate the x positions for the interval
                max_length = pat_chr_size_df['length'].max()
                x_start = pos / max_length
                x_end = end / max_length
                axs[2].axhspan(y_pos-0.35, y_pos+0.35, xmin=x_start, xmax=x_end, facecolor='black', alpha=1,
                               linewidth=3,edgecolor=".1")
                axs[2].set_xlim([0,max_length])
        # Add number to the end of chromosome
        for chr in set([i.split(':')[0] for i in hg002_pat_interval_list]):
            y_pos = pat_chr_size_df.index[pat_chr_size_df['chr'] == chr][0]
            axs[2].text(pat_chr_size_df[pat_chr_size_df['chr'] == chr]['length']+5000000, y_pos+0.25,
                        [i.split(':')[0] for i in hg002_pat_interval_list].count(chr))
        #plt.show()
        plt.suptitle(self.hg38_interval_list[0], fontsize=18)
        plt.tight_layout()
        plt.savefig(f"{self.hg38_interval_list[0]}_karyoplot.png", dpi=600)

    def parse_interval(self, interval_str):
        parts = interval_str.split(':')
        chromosome = parts[0]
        start, end = map(int, parts[1].split('-'))
        return chromosome, start, end

    def plot_intervals_comparison(self, ratio=6, fig_height=12, flanking=False, save=False, savepath=None):
        hg002_mat_interval_list = [i for i in self.hg002_interval_list if 'MATERNAL' in i]
        hg002_pat_interval_list = [i for i in self.hg002_interval_list if 'PATERNAL' in i]

        # Calculate event length to determine
        event_chromosome, event_start, event_end = self.parse_interval(self.hg38_interval_list[0])
        event_length = abs(event_end-event_start)

        # Ratio is used to determine the event length and x axis length
        x_axis_length = event_length*ratio

        f, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(35, fig_height))

        # Set default x axis if no interval is provided
        default_start = 0
        default_end = 100000  # set a default value that makes sense in your context
        default_length = default_end - default_start

        # Get a colormap for mat
        hg38_colormap = cm.Paired
        # Normalize to mat interval count
        hg38_norm = plt.Normalize(0, len(self.hg38_interval_list))

        # Get hg38 median interval to determine x lim
        median_hg38_chr, median_hg38_start, median_hg38_end = self.parse_interval(sorted(self.hg38_interval_list, key=lambda x: self.parse_interval(x)[1])[int(round(len(self.hg38_interval_list)/2,0))-1])
        hg38_center = (median_hg38_start+median_hg38_end)/2
        hg38_start = hg38_center-x_axis_length/2
        hg38_end = hg38_center+x_axis_length/2
        # Plot Chrom
        ax1.plot([hg38_start,hg38_end],[0,0], linewidth=5, color="gray",alpha=0.5)
        # Plot Each Interval
        for i, interval_str in enumerate(sorted(self.hg38_interval_list, key=lambda x: self.parse_interval(x)[1])):
            chromosome, start, end = self.parse_interval(interval_str)
            if flanking and i == 1:
                # Plot event interval
                ax1.plot([start, end], [0, 0], linewidth=20, solid_capstyle='butt', label=f'{chromosome}:{format(start, ",")}-{format(end, ",")}', color='black')
                # Offset the annotate text to the left
                ax1.annotate(f'{event_length}bp', xy=((start+end)/2, 0), xytext=(((start+end)/2)-event_length*0.15, 0.1), fontsize=30, fontweight='bold')
            # Plot Left Flanking Region
            elif flanking and i ==2:
                ax1.plot([start, end], [0, 0], linewidth=20, solid_capstyle='butt', label=f'{chromosome}:{format(start, ",")}-{format(end, ",")}', color=hg38_colormap(hg38_norm(i-1)))
            # Plot Right Flanking Region or others
            else:
                ax1.plot([start, end], [0, 0], linewidth=20, solid_capstyle='butt', label=f'{chromosome}:{format(start, ",")}-{format(end, ",")}', color=hg38_colormap(hg38_norm(i)))


        ax1.set_yticks([0], ["hg38"], fontweight='bold',fontsize=28)
        ax1.set_ylim([-1,1])
        ax1.set_xlabel('')
        ax1.set_xlim([hg38_start,hg38_end])
        ax1.set_xticks([])

        # Get a colormap for mat
        mat_colormap = cm.Pastel1
        # Normalize to mat interval count
        mat_norm = plt.Normalize(0, len(hg002_mat_interval_list))

        # Handle the case where there is no interval
        if len(hg002_mat_interval_list) == 0:
            mat_start = default_start
            mat_end = default_end
        else:
            # Get the center of hg2 mat intervals to determine x lim
            mat_left_end_interval_pos = self.parse_interval(sorted(hg002_mat_interval_list, key=lambda x: self.parse_interval(x)[1])[0])[1]
            mat_right_end_interval_end = self.parse_interval(sorted(hg002_mat_interval_list, key=lambda x: self.parse_interval(x)[1])[-1])[2]

            mat_center = (mat_left_end_interval_pos+mat_right_end_interval_end)/2
            mat_start = mat_center-x_axis_length/2
            mat_end = mat_center+x_axis_length/2
        # Plot Chrom
        ax2.plot([mat_start,mat_end],[0,0], linewidth=5, color="gray",alpha=0.5)
        # Plot Each Interval
        for i, interval_str in enumerate(sorted(hg002_mat_interval_list, key=lambda x: self.parse_interval(x)[1])):
            chromosome, start, end = self.parse_interval(interval_str)

            # Plot a line for the interval (if flanking use hg38 colormap to match flanking color code)
            if flanking:
                ax2.plot([start, end], [0, 0], linewidth=20, solid_capstyle='butt', label=f'{chromosome}:{format(start, ",")}-{format(end, ",")}', color=hg38_colormap(hg38_norm(i)))
            else:
                ax2.plot([start, end], [0, 0], linewidth=20, solid_capstyle='butt', label=f'{chromosome}:{format(start, ",")}-{format(end, ",")}', color=mat_colormap(mat_norm(i)))


        ax2.set_yticks([0], ["hg002_mat"], fontweight='bold',fontsize=28)
        ax2.set_ylim([-1,1])
        ax2.set_xlabel('')
        ax2.set_xlim([mat_start,mat_end])
        ax2.set_xticks([])
        # If flanking, add distance between the left and right flanking region
        if flanking:
            # Get the end of the LFR
            mat_LRF_end = self.parse_interval(sorted(hg002_mat_interval_list, key=lambda x: self.parse_interval(x)[1])[0])[2]
            mat_RRF_start = self.parse_interval(sorted(hg002_mat_interval_list, key=lambda x: self.parse_interval(x)[1])[1])[1]
            # Calculate the distance between the two
            if mat_LRF_end == mat_RRF_start:
                mat_flanking_distance = 0
            else:
                mat_flanking_distance = mat_RRF_start-mat_LRF_end -1
            # Annotate the distance between the two flanking regions
            ax2.annotate(f'{mat_flanking_distance}bp', xy=((mat_LRF_end + mat_RRF_start) / 2, 0), xytext=(((mat_LRF_end + mat_RRF_start)/2)-mat_flanking_distance*0.15 , 0.1), fontsize=30, fontweight='bold')

        # Normalize to mat interval count
        # Get a colormap
        pat_colormap = cm.Pastel2
        pat_norm = plt.Normalize(0, len(hg002_pat_interval_list))
        # Get hg2 pat median interval to determine x lim
        # Handle the case where there is no interval
        if len(hg002_pat_interval_list) == 0:
            pat_start = default_start
            pat_end = default_end
        else:
            # Get the center of hg2 pat intervals to determine x lim
            pat_left_end_interval_pos = self.parse_interval(sorted(hg002_pat_interval_list, key=lambda x: self.parse_interval(x)[1])[0])[1]
            pat_right_end_interval_end = self.parse_interval(sorted(hg002_pat_interval_list, key=lambda x: self.parse_interval(x)[1])[-1])[2]

            pat_center = (pat_left_end_interval_pos + pat_right_end_interval_end) / 2
            pat_start = pat_center-x_axis_length/2
            pat_end = pat_center+x_axis_length/2
        # Plot Chrom
        ax3.plot([pat_start,pat_end],[0,0], linewidth=5, color="gray",alpha=0.5)
        # Plot Each Interval
        for i, interval_str in enumerate(sorted(hg002_pat_interval_list, key=lambda x: self.parse_interval(x)[1])):
            chromosome, start, end = self.parse_interval(interval_str)

            # Plot a line for the interval (if flanking use hg38 colormap to match flanking color code)
            if flanking:
                ax3.plot([start, end], [0, 0], linewidth=20, solid_capstyle='butt', label=f'{chromosome}:{format(start, ",")}-{format(end, ",")}', color=hg38_colormap(hg38_norm(i)))
            else:
                ax3.plot([start, end], [0, 0], linewidth=20, solid_capstyle='butt', label=f'{chromosome}:{format(start, ",")}-{format(end, ",")}', color=pat_colormap(pat_norm(i)))


        ax3.set_yticks([0], ["hg002_pat"], fontweight='bold',fontsize=28)
        ax3.set_ylim([-1,1])
        ax3.set_xlabel('')
        ax3.set_xlim([pat_start,pat_end])
        ax3.set_xticks([])
        # If flanking, add distance between the left and right flanking region
        if flanking:
            # Get the end of the LFR
            pat_LRF_end = self.parse_interval(sorted(hg002_pat_interval_list, key=lambda x: self.parse_interval(x)[1])[0])[2]
            pat_RRF_start = self.parse_interval(sorted(hg002_pat_interval_list, key=lambda x: self.parse_interval(x)[1])[1])[1]
            # Calculate the distance between the two
            if pat_LRF_end == pat_RRF_start:
                pat_flanking_distance = 0
            else:
                pat_flanking_distance = pat_RRF_start-pat_LRF_end-1
            # Annotate the distance between the two flanking regions
            ax3.annotate(f'{pat_flanking_distance}bp', xy=((pat_LRF_end + pat_RRF_start) / 2, 0), xytext=(((pat_LRF_end + pat_RRF_start)/2)-pat_flanking_distance*0.15 , 0.1), fontsize=30, fontweight='bold')

        # Adding legends for the coordinates
        legend1 = ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=28)
        legend1.set_title("hg38", prop={"size": 25, "weight": "bold"})
        legend2 = ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=28)
        legend2.set_title("HG002 Maternal", prop={"size": 25, "weight": "bold"})
        legend3 = ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=28)
        legend3.set_title("HG002 Paternal", prop={"size": 25, "weight": "bold"})

        # Adjust layout to fit legend
        plt.tight_layout(rect=[0, 0, 0.85, 1])

        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.9)
        if save:
            if savepath is None:
                savepath = os.getcwd()
            else:
                savepath = savepath
            plt.savefig(f"{savepath}/{self.hg38_interval_list[0]}_interval_comparison.png", dpi=600)
        else:
            plt.show()


def compute_dot_matrix(seq1: str, seq2: str, k: int) -> tuple:
    coords_x, coords_y = [], []
    hash_table = defaultdict(list)

    for i in range(len(seq1) - k + 1):
        hash_table[seq1[i:i+k]].append(i)

    for j in range(len(seq2) - k + 1):
        for i in hash_table.get(seq2[j:j+k], []):
            coords_x.append(i)
            coords_y.append(j)

    return np.array(coords_x), np.array(coords_y)


def plot_dot_matrix(coords_x: np.ndarray, coords_y: np.ndarray, seq1_name: str, seq2_name: str) -> None:
    plt.scatter(coords_x, coords_y, s=.1, c='black')
    plt.gca().invert_yaxis()
    plt.xlabel(seq1_name)
    plt.ylabel(seq2_name)
    plt.show()


def dot_plot(seq1: str, seq2: str, seq1_name: str = "Sequence 1", seq2_name: str = "Sequence 2", k: int = 50) -> None:
    if len(seq1) * len(seq2) > 1e8:
        raise ValueError("Sequences are too large; adjust 'k' or trim sequences")

    coords_x, coords_y = compute_dot_matrix(seq1, seq2, k)
    plot_dot_matrix(coords_x, coords_y, seq1_name, seq2_name)

