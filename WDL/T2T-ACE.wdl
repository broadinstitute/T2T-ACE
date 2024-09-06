version 1.0

workflow T2T_ACE{
    input {
        String SampleName
        File? CNV_VCF
        File? DEL_bed
        File? DUP_bed
        File T2T_Reference
        File hg38_Reference
    }
    call T2T_ACE {
        input:
            SampleName = SampleName,
            CNV_VCF = CNV_VCF,
            T2T_Reference = T2T_Reference,
            hg38_Reference = hg38_Reference,
            DEL_bed = DEL_bed,
            DUP_bed = DUP_bed
    }
    output {
        File? DEL_eval = T2T_ACE.DEL_eval_sum
        File? DUP_eval = T2T_ACE.DUP_eval_sum
    }
    meta {
        author: "Yueyao Gao"
        email: "tag@broadinstitute.org"
        description: "T2T_ACE.wdl is design to run T2T-ACE on a single sample. The input files are a CNV vcf file, a DEL bed file, a DUP bed file, a T2T reference fasta file and a hg38 reference fasta file. The output file is a sheet with T2T-ACE results."
    }
}

task T2T_ACE {
    input {
        String SampleName
        File? CNV_VCF
        File? DEL_bed
        File? DUP_bed
        File T2T_Reference
        File hg38_Reference
        String docker = "us.gcr.io/tag-public/t2t-ace:0.0.0"
        String test = false
        Int memory = 64
        Int cpu = 8
        Int disk_space_gb = 500
        Int preemptible = 3
    }
    command <<<
        set -e

        conda run --no-capture-output -n T2T_ACE_env python3 /BaseImage/T2T-ACE/run_T2T-ACE.py \
        ~{'--cnv_vcf '+ CNV_VCF} \
        ~{'--del_txt '+ DEL_bed} \
        ~{'--dup_txt '+ DUP_bed} \
        --t2t_ref ~{T2T_Reference} \
        --hg38_ref ~{hg38_Reference} \
        --test ~{test}

        if [ -f output_DEL_eval_sum.csv ]; then
            mv output_DEL_eval_sum.csv ~{SampleName}_DEL_eval_sum.csv
        fi
        if [ -f output_DUP_eval_sum.csv ]; then
            mv output_DUP_eval_sum.csv ~{SampleName}_DUP_eval_sum.csv
        fi

    >>>
    output {
        File? DEL_eval_sum = "~{SampleName}_DEL_eval_sum.csv"
        File? DUP_eval_sum = "~{SampleName}_DUP_eval_sum.csv"
    }
    runtime{
        docker: docker
        memory: memory + " GB"
        cpu: cpu
        disk: "local-disk " + disk_space_gb + " SSD"
        preemptible: preemptible
    }
}


