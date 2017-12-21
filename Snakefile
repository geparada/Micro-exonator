configfile : "config.yaml"

############ Params #############################



##################################################

import glob, os

#DATA = [os.path.basename(f).split(".")[0] for f in glob.glob('data/input/*')]
DATA = [os.path.basename(f).split(".")[0] for f in config["samples"]]


rule Splice_Junction_Library:
    input:
        config["Genome_fasta"],
        config["Gene_anontation_fasta"],
        config["Gene_anontation_bed12"]
    params:
        ME_len = config["ME_len"]
    output:
        "Round1/ME_TAGs.fa"
    shell:
        "python2 src/SJ_tags_generator_for_micro_exons.py {input} {params.ME_len} > {output}"


# rule sra_to_fastq:
#     input:
#         config["input_dir"] + "/{sample}.sra"
#     output:
#         temp("data/fastq_paired/{sample}.fastq")
#     shell:
#         "fastq-dump {input} -O data/fastq_paired/"


# rule fastq_gz_to_fastq:
#     input:
#         config["input_dir"] + "/{sample}.fastq.gz"
#     output:
#         temp("data/fastq/{sample}.fastq")
#     shell:
#         "gzip -dc {input} > {output}"
#
# rule fastq_input:
#     input:
#         config["input_dir"] + "/{sample}.fastq"
#     output:
#         "data/fastq/{sample}.fastq"
#     shell:
#         "ln -s {input} {output}"

rule download_to_fastq:
    input:
        "download/{sample}.download.sh"
    output:
        "data/fastq/{sample}.fastq"
    shell:
        "bash {input}"


# rule split_fastq:
#     input:
#         "data/fastq_paired/{sample}.fastq"
#     output:
#         temp("data/fastq/{sample}.fastq")
#     shell:
#         "python2 src/split_paired_end.py {input} > {output}"


rule bwa_index:
    input:
        "Round1/ME_TAGs.fa"
    output:
        "Round1/ME_TAGs.fa.amb"
    shell:
        "bwa index {input}"

rule Round1_bwa_mem_to_tags:
    input:
        "Round1/ME_TAGs.fa",
        "data/fastq/{sample}.fastq",
        "Round1/ME_TAGs.fa.amb"
    output:
        temp("Round1/{sample}.sam")
    threads: 5
    shell:
        "bwa mem -t {threads} -O 2,6 -L 25 {input[0]} {input[1]} > {output}"


rule Round1_alingment_pre_processing:
    input:
        "Round1/{sample}.sam"
    output:
        "Round1/{sample}.sam.pre_processed"
    shell:
        "python2 src/alingment_pre_processing.py {input} F > {output}"



#############################################################################

                     #Round1 - POST-Processing###

#############################################################################




rule row_Micro_Exon_reads:
    input:
        config["Genome_fasta"],
        "Round1/{sample}.sam.pre_processed"
    output:
        "Round1/{sample}.sam.row_ME",
        "Round1/{sample}.sam.row_ME.fastq"
    shell:
        "python2 src/row_ME.py {input} > {output[0]}"


rule hisat2_Genome_index:
    input:
        config["Genome_fasta"]
    output:
        "data/Genome.1.ht2"
    threads: 5
    shell:
        "hisat2-build {input} data/Genome"


rule hisat2_to_Genome:
    input:
        "Round1/{sample}.sam.row_ME.fastq",
        "data/Genome.1.ht2"
    output:
        "Round1/{sample}.sam.row_ME.Genome.Aligned.out.sam"
    threads: 1
    shell:
        "hisat2 -x data/Genome -U {input[0]} > {output}"


rule Round1_filter:
    input:
        config["Genome_fasta"],
        "Round1/{sample}.sam.row_ME",
        "Round1/{sample}.sam.row_ME.Genome.Aligned.out.sam",
        config["GT_AG_U2_5"],
        config["GT_AG_U2_3"],
        config["vertebrates_phylop"]
    params:
        ME_len = config["ME_len"]
    output:
        "Round1/{sample}.sam.row_ME.filter1"
    shell:
        "python2 src/ME_filter1.py {input} {params.ME_len} > {output}"


rule Micro_Exon_table:
    input:
        expand("Round1/{sample}.sam.row_ME.filter1", sample=DATA )
    output:
        "Round1/TOTAL/TOTAL.sam.row_ME.filter1.ME_centric"
    shell:
        "cat {input} > Round1/TOTAL/TOTAL.sam.row_ME.filter1 && python2 src/ME_centric_table.py Round1/TOTAL/TOTAL.sam.row_ME.filter1 > {output}"


rule Micro_Exon_Tags:
    input:
        "Round1/ME_TAGs.fa",
        "Round1/TOTAL/TOTAL.sam.row_ME.filter1.ME_centric"
    output:
        "Round2/ME_canonical_SJ_tags.fa"
    shell:
        "python2 scr/Micro_exons_tags.py  {input} > {output}"



#############################################################################

                     #Round2#

#############################################################################

rule Round2_bowtie_tags_index:
    input:
        "Round2/ME_canonical_SJ_tags.fa"
    output:
        "Round2/ME_canonical_SJ_tags.fa.1.ebwt"
    shell:
        "bowtie-build {input} {input}"


rule Round2_bowtie_to_tags:
    input:
        "Round2/ME_canonical_SJ_tags.fa",
        "data/fastq/{sample}.fastq",
        "Round2/ME_canonical_SJ_tags.fa.1.ebwt"
    output:
        temp("Round2/{sample}.sam")
    threads: 5
    shell:
        "bowtie {input[0]} -p {threads} -q {input[1]} -S -v 2 > {output}"


rule Round2_alingment_pre_processing:
    input:
        "Round2/{sample}.sam"
    output:
        "Round2/{sample}.sam.pre_processed"
    shell:
        "python2 src/alingment_pre_processing_round2_bowtie.py {input} F > {output}"



#############################################################################

                     #Round2 - POST-Processing###

#############################################################################


rule ME_reads:
    input:
        "Round2/{sample}.sam.pre_processed"
    output:
        "Round2/{sample}.sam.pre_processed.fastq"
    shell:
        "python2 src/round2_ME_reads_fastq.py {input} > {output}"


rule bowtie_Genome_index:
    input:
        config["Genome_fasta"]
    output:
        config["Genome_fasta"] + ".1.ebwt"
    shell:
        "bowtie-build {input} {input}"

rule bowtie_to_genome:
    input:
        "Round2/{sample}.sam.pre_processed.fastq",
        config["Genome_fasta"],
        config["Genome_fasta"] + ".1.ebwt"
    output:
        "Round2/{sample}.sam.pre_processed.hg19.sam"
    shell:
        "bowtie {input[1]} -p 1 -q {input[0]} -S -v 2| awk '$2==0 || $2==16'> {output}"


rule Round2_filter:
    input:
        "Round2/{sample}.sam.pre_processed",
        "Round2/{sample}.sam.pre_processed.hg19.sam",
    output:
        "Round2/{sample}.sam.pre_processed.filter1"
    shell:
        "python2 src/Filter1_round2.py {input} > {output}"


rule ME_SJ_coverage:
    input:
        "Round2/ME_canonical_SJ_tags.fa",
        "Round1/TOTAL/TOTAL.sam.row_ME.filter1.ME_centric",
        config["Gene_anontation_bed12"],
        "Round2/{sample}.sam.pre_processed.filter1"
    params:
        ME_len = config["ME_len"]
    output:
        "Round2/{sample}.sam.pre_processed.filter1.ME_SJ_coverage"
    shell:
        "python2 src/ME_SJ_coverage.py {input} {params.ME_len} > {output}"




rule Total_sample_exon_counts:
    input:
        expand("Round2/{sample}.sam.pre_processed.filter1.ME_SJ_coverage", sample=DATA )
    output:
        "Round2/TOTAL.filter1.ME_SJ_coverage"
    shell:
        "cat {input} > {output}"


rule Output:
    input:
        "Round1/TOTAL/TOTAL.sam.row_ME.filter1.ME_centric",
        "Round2/TOTAL.filter1.ME_SJ_coverage"
    output:
        "Report/out_filtered_ME.txt",
        "Report/out_low_scored_ME.txt",
        "Report/out_shorter_than_3_ME.txt",
        "Report/report.html"
    shell:
        "R -e" + ' rmarkdown::render("src/final_filters.Rmd", params = list(ME_table={input[0]}, ME_coverage={input[1]}, out_filtered_ME={output[0]}, out_low_scored_ME={output[1]}, out_shorter_than_3_ME={output[2]}), output_file={output[3]})'


snakemake -j 30000 Round2/TOTAL.filter1.ME_SJ_coverage --cluster-config cluster.json --cluster "bsub -n {cluster.nCPUs} -R {cluster.resources} -c {cluster.tCPU} -G {cluster.Group} -q {cluster.queue} -o {cluster.output} -e {cluster.error} -M {cluster.memory}"

# snakemake --dag Round2/TOTAL.filter1.ME_SJ_coverage | dot -Tsvg > Micro-Exonator.svg

# snakemake --dag Report/report.html | dot -Tsvg > Micro-Exonator.svg
