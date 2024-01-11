import os


class ResourceLocations:
    if os.environ.get('USER') == 'gaoyueya':
        hg002t2t = "/Users/gaoyueya/Documents/Reference_Genome/hg002v1.0.fasta"
        chm13 = "/Users/gaoyueya/Documents/Reference_Genome/chm13v2.0.fa.gz"
        hg38 = "/Users/gaoyueya/Documents/Reference_Genome/Homo_sapiens_assembly38.fasta"
        HG2_DRAGEN_cnv_path = "/Users/gaoyueya/Documents/Projects/TAG-Ticket1639/DRAGEN4_2_4_visualization/DRAGEN_Output/NA24385.cnv_sv.vcf"

    elif os.environ.get('USER') == 'fleharty':
        hg002t2t = "/Users/fleharty/resources/hg002v1.0.fasta.gz"
        chm13 = "/Users/fleharty/resources/chm13v2.0.fa.gz"
        hg38 = "/Users/fleharty/resources/Homo_sapiens_assembly38.fasta"
        HG2_DRAGEN_cnv_path = "/Users/fleharty/resources/NA24385.cnv_sv.vcf.gz"
