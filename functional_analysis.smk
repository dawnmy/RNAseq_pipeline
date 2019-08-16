from os import path

# Load paths for input and output from config
configfile: "config/config.yaml"

dataset = config["dataset"]
proj_dir = path.join(config["out_dir"], dataset)
ref = config["ref"]
threads = config["threads"]
kegg = config["kegg"]
gene_ko_map = config["gene_ko_map"]
kofam_dir = config["kofamscan"]

out_data_dir = path.join(proj_dir, "data")
reports_dir = path.join(proj_dir, "reports")
results_dir = path.join(proj_dir, "results")

# samples, = glob_wildcards(fq_dir + "/{sample}.fastq.gz")

# wildcard_constraints:
#    sample = "[^\.\/]+"


rule all:
    input:
        #out_data_dir + "/annotation/gene_family_euk_kegg.blast.out"
        #out_data_dir + "/annotation/gene_family_euk_kegg.diamond.txt",
        #out_data_dir + "/annotation/gene_family_euk_kegg.kofam.out"
        results_dir + \
            "/count/{dataset}.diamond.ko.count".format(dataset=dataset)

# Option 1
# Annotation the CDS with peptides from KEGG euk using blast
rule blast:
    input:
        cds = ref,
        db = kegg
    output:
        out_data_dir + "/annotation/gene_family_euk_kegg.diamond.out"
    threads: threads
    conda:
        "config/conda.blast.yaml"
    shell:
        """
        blastx -num_threads {threads} -max_target_seqs 3 -query {input.cds} \
            -db {input.db} -outfmt "7 std stitle qcovs" -out {output}
        """

# Option 2
# Annotation the CDS with peptides from KEGG euk using blast
rule diamond:
    input:
        cds = ref,
        db = kegg
    output:
        daa = out_data_dir + "/annotation/gene_family_euk_kegg.diamond.daa",
        tab = out_data_dir + "/annotation/gene_family_euk_kegg.diamond.txt"
    threads: threads
    log:
        reports_dir + "/diamond/gene_family_euk_kegg.diamond.log"
    conda:
        "config/conda.blast.yaml"
    shell:
        """
        diamond blastx -d {input.db} --more-sensitive -p {threads} -v --log \
            -q {input.cds} --header --id 50 -e 1e-5 --query-cover 60 \
            -a {output.daa} >> {log}
        diamond view -p {threads} -a {output.daa} -f 6 qseqid sseqid pident evalue \
            length bitscore qcovhsp > {output.tab}
        """

# Option 3
# Annotate cds using kofam
rule kofam:
    input:
        cds = config["ref_pep"],
        ko = config["ko_list"]
    output:
        out_data_dir + "/annotation/gene_family_euk_kegg.kofam.out"
    threads: threads
    params:
        kofam_dir = config["kofamscan"],
        profile = config["profile"]
    shell:
        """
        {params.kofam_dir}/exec_annotation {input.cds} --cpu {threads} -p \
            {params.profile} -k {input.ko} -o {output}
        """


# The annotation option 2 is used here
# Map the CDS ID to KO ID
rule map_ko:
    input:
        diamond_anno = rules.diamond.output.tab,
        gene_ko = gene_ko_map,
        exp_tab = results_dir + "/count/{dataset}.count"
    output:
        gene_ko_exp = results_dir + \
            "/count/{dataset}.diamond.gene.ko.anno.count",
        # ko_exp = results_dir + "/count/{dataset}.diamond.ko.count
    conda:
        "config/conda.rnaseq.yaml"
    shell:
        """
        csvtk join -Ttkf "gene" {input.exp_tab} <(csvtk join -Ttkf "kegg" \
            <(echo "gene\tkegg";cut -f 1-2 {input.diamond_anno}) \
                <(echo "kegg\tko";cat {input.gene_ko})|\
                    awk '$3!="" && !a[$1]++') > {output.gene_ko_exp}
        """

# Group gene expression table to KO expression table
rule make_ko_exp:
    input:
        rules.map_ko.output.gene_ko_exp
    output:
        results_dir + "/count/{dataset}.diamond.ko.count"
    run:
        import pandas as pd
        df = pd.read_csv(str(input), sep="\t", index_col=0,
                         usecols=lambda column: column != "kegg")
        df.fillna("noko", inplace=True)
        ko_exp = df.groupby("ko").sum()
        ko_exp.to_csv(str(output), sep="\t")


# DE analysis based on KO genes


# KEGG pathway enrichment analysis
