from os import path

# Load paths for input and output from config
configfile: "config/config.yaml"

dataset = config["dataset"]
fq_dir = config["fq_dir"]
proj_dir = path.join(config["out_dir"], dataset)
ref = config["ref"]
threads = config["threads"]

out_data_dir = path.join(proj_dir, "data")
qc_dir = path.join(out_data_dir, "qc_fq")
reports_dir = path.join(proj_dir, "reports")
results_dir = path.join(proj_dir, "results")

samples, = glob_wildcards(fq_dir + "/{sample}.fastq.gz")

wildcard_constraints:
    sample = "[^\.\/]+"


def get_pe_fq(wc):
    return [path.join(fq_dir, wc.sample + end + ".fastq.gz") for end in ["_R1", "_R2"]]


def get_se_fq(wc):
    return path.join(fq_dir, wc.sample + ".fastq.gz")


rule all:
    input:
        count = expand(
            results_dir + "/count/{dataset}.count", dataset=dataset)

rule fastp:
    input:
        get_se_fq
    output:
        reads = qc_dir + "/{sample}.qc.fq.gz",
        html = reports_dir + "/fastp/{sample}.qc.report.html"
    threads: 8
    conda:
        "config/conda.rnaseq.yaml"
    shell:
        """
        fastp -i {input} -o {output.reads} -5 20 -3 20 -l 30 -n 3 -h {output.html} -w {threads}
        """

rule build_idx:
    input:
        ref = ref
    output:
        fa_idx = ref + ".fai",
        bwa_idx = ref + ".bwt"
    conda:
        "config/conda.rnaseq.yaml"
    shell:
        """
        samtools faidx {input.ref}
        bwa index {input.ref}
        """


rule bwa:
    input:
        reads = rules.fastp.output.reads,
        ref = ref,
        ref_idx = rules.build_idx.output.bwa_idx
    output:
        out_data_dir + "/bam/{sample}.bam"
    conda:
        "config/conda.rnaseq.yaml"
    benchmark:
        reports_dir + "/benchmarks/{sample}.bwa.txt"
    log:
        reports_dir + "/bwa/{sample}.log"
    threads: threads
    shell:
        """
        bwa mem -t {threads} {input.ref} {input.reads} |\
            samtools view -@ {threads} -Shb - 2>> {log} |\
            samtools sort -@ {threads} -m 20G - -o {output} >> {log} 2>&1
        """

rule bam_filter:
    input:
        rules.bwa.output
    output:
        bam = out_data_dir + "/bam/{sample}.filtered.bam",
        bai = out_data_dir + "/bam/{sample}.filtered.bam.bai"
    conda:
        "config/conda.rnaseq.yaml"
    log:
        reports_dir + "/samtools/{sample}.log"
    threads: threads
    shell:
        """
        samtools view -@ {threads} -Shb -F 4 -q 10 {input} 2>> {log} |\
            samtools sort -@ {threads} -m 20G - -o {output.bam} >> {log} 2>&1
        samtools index {output.bam}
        """

rule samcount:
    input:
        bam = rules.bam_filter.output.bam,
        bai = rules.bam_filter.output.bai
    output:
        results_dir + "/count/{sample}.samcount"
    conda:
        "config/conda.rnaseq.yaml"
    threads: threads
    shell:
        """
        samtools idxstats {input.bam} > {output}
        """

rule cat_count:
    input:
        expand(
            results_dir + "/count/{sample}.samcount", sample=samples)
    output:
        results_dir + "/count/{dataset}.count"
    conda:
        "config/conda.rnaseq.yaml"
    params:
        cols = ",".join(["1"] + [str(3 * i)
                                 for i in range(1, len(samples) + 1)])
    threads: threads
    shell:
        """
        cat <(printf "gene\t";sed 's/[^ ]\+\/\|\.samcount//g;s/ \+/\t/g' <(echo {input})) \
            <(csvtk join -TtHf 1 {input}|csvtk cut -TtHf {params.cols}) > {output}
        """
