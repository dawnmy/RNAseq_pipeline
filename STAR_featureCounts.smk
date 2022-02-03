rule all:
    input:
        expand(results_dir + "/featurecounts/{sample}.star_with_xs.genome.featurecounts", sample=samples)


rule fastp:
    input:
        get_pe_fq
    output:
        or1 = qc_dir + "/{sample}.fastp.1.fq.gz",
        or2 = qc_dir + "/{sample}.fastp.2.fq.gz",
        html = reports_dir + "/fastp/{sample}.fastp.report.html"
    threads: 8
    shell:
        """
        fastp -i {input[0]} -I {input[1]} -o {output.or1} -O {output.or2} \
            -5 20 -3 20 -l 50 -h {output.html} -w {threads}
        """

rule star:
    input:
        r1 = rules.fastp.output.or1,
        r2 = rules.fastp.output.or2,
    output:
        bam = results_dir + "/star_with_xs/{sample}.Aligned.sortedByCoord.out.bam"
        # done = touch(results_dir + "/star/{sample}.done")
    params:
        ref = config['star_genome_ref'],
        out_prefix = results_dir + "/star_with_xs/{sample}."
    # conda:
        # "config/conda.rnaseq.yaml"
    benchmark:
        reports_dir + "/benchmarks/{sample}.star_with_xs.txt"
    threads: threads
    shell:
        """
        STAR --genomeDir {params.ref} \
            --runThreadN {threads} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.out_prefix} \
            --quantMode TranscriptomeSAM GeneCounts \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --outSAMstrandField intronMotif \
            --outSAMattributes Standard
        """

rule featurecounts:
    input: 
        bam = rules.star.output.bam,
        gff = config['gff']
    output:
        results_dir + "/featurecounts/{sample}.star_with_xs.genome.featurecounts" 
    threads: threads
    shell:
        """
        featureCounts -p -T {threads} --donotsort -t gene -g ID -a {input.gff} -o {output} {input.bam}
        """
