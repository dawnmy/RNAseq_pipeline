## RNAseq data processing and transcriptome analysis for P. minimum (co-culture and temperature experiment).

Currently it is for single-end rna-sequencing data, and nanopre long read rna-seq data.
This pipeline can be adapted to other RNAseq analysis.


### Prerequirements

To reproduce the output, you need to use `Bioconda`.

Please follow the instruction [here](https://bioconda.github.io) to install `Bioconda`. 
And then you need to install `snakemake` and Python package `click` and `pandas`:

```shell
conda install snakemake=5.5.4
conda install Click=7.0
conda install pandas=0.25.0
```


After this has been done, download the pipeline onto your system:

```shell
git clone git@github.com:dawnmy/RNAseq_pipeline.git
```

### Modify the config file: `config/config.yaml
All the paths must be either relative path to the parent directory of `config` folder or absolute path.

```yaml
dataset: CPm # name for the dataset
fq_dir: ../data/seq
out_dir: ../outputs
ref: ../ref/Prorocentrum-minimum-CCMP1329.cds.fa
ref_pep: ../ref/Prorocentrum-minimum-CCMP1329.cds.fa
kegg: <path to the family_eukaryotes.pep file> # Please create the corresponding database or index if you use diamond or blastx
gene_ko_map: <path to the KEGG gene and KO ID map file genes_ko.list>  
kofamscan: <dir to the exe of kofamscan>
is_long_read: false # Is it long read rnaseq data

# Path to your KO-HMM database
# A database can be a .hmm file, a .hal file or a directory in which
# .hmm files are. Omit the extension if it is .hal or .hmm file
profile: <dir to kofamscan profiles>

# Path to the KO list file
ko_list: <path to the kofamscan ko_list file>

threads: 20
```


### Run the pipelines

#### Get the expression table for genes

```shell
snakemake -s rnaseq.smk -j 20 --use-conda
```
`-s` to specify the pipeline file, and `-j` to set the number of threads to use and `--use-conda` to 
let the pipeline install required softwares with specified version. The `conda` ENVs will be created under 
the path of the program by default. The program may take ten minutes to create the ENV for the first time.
If you do not wish to create the conda ENV in the working directory, 
please use --conda-prefix parameter to specify the desired path to create the `conda` ENV. 


**If you use SGE for the job submission, you can use the following cmd:**

```shell
snakemake -s rnaseq.smk --latency-wait 30 --use-conda -c "qsub -cwd -q <the job submission queue> \
 -pe multislot {threads} -i /dev/null -e <dir for std error logs> -o <dir for std output logs> \
 -v PATH" -j 2
```

#### Make the KO gene expression table

1. Annotate the genes using KEGG peptide sequences (optional)

You can skip this step if you want to use the gene KO annotation 
file provided in this repo in: `data/annotation/gene_family_euk_kegg.diamond.txt`.
Then you should copy this file into the `<out_dir>/<dataset>/data/annotation/` directory.
If the folder does not exist, please create it.

2. Map the KEGG annotation, KO ID to the gene expression table to make a KO gene expression table

```shell
snakemake -s functional_analysis.smk -j 10 --use-conda
```


**The R scripts for DE, PCA analysis and the KEGG pathway enrichment analysis are under scripts folder**.
Please modify the script (input, output, figure file name, and the group information) to adatpt it to your own case. 
It is recommended to run the R scripts in an interactive way in your local PC for better data understanding. 

### The output structure

```
outputs
└── CPm
    ├── data
    │   ├── annotation
    │   ├── bam
    │   └── qc_fq
    ├── reports
    │   ├── benchmarks
    │   ├── bwa
    │   ├── diamond
    │   ├── fastp
    │   └── samtools
    └── results
        └── count
```



