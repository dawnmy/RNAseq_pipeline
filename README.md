## A pipeline for RNAseq data processing and DE analysis.

Currently it is only for single-end sequencing data, but can be easyly adapted for paired-end


### Prerequirements

To reproduce the output, you need to use `Bioconda`.

Please follow the instruction [here](https://bioconda.github.io) to install `Bioconda`. 
And then you need to install `snakemake` and Python package `click`:

```shell
conda install snakemake=5.5.4
conda install Click=7.0
```

After this has been done, download the pipeline onto your system:

```shell
git clone git@github.com:dawnmy/RNAseq_pipeline.git
```

### Modify the config file: `config/config.yaml
All the paths must be either relative path to the parent directory of `config` folder or absolute path.

```yaml
dataset: CPm # name for the dataset
fq_dir: /vol/projects/dzhiluo/Pminimum/data/seq
out_dir: /vol/projects/dzhiluo/Pminimum/outputs
ref: /vol/projects/dzhiluo/Pminimum/ref/Prorocentrum-minimum-CCMP1329.cds.fa
threads: 20
```


### Run the pipelines

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

**The DE and PCA analysis R script is under scripts folder** 

### The output structure

```
outputs
└── CPm
    ├── data
    │   ├── bam
    │   └── qc_fq
    ├── reports
    │   ├── benchmarks
    │   ├── bwa
    │   ├── fastp
    │   └── samtools
    └── results
        └── count
```



