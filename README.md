# The vdB 16S Pipeline

## Prerequisites

- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [Apptainer/Singularity](https://apptainer.org/): while in many cases we do provide conda envs the only method of execution we support is via containers.
- (optional) [A Snakemake Profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles): this coordinates the execution of jobs on whatever hardware you are using.

> Please set your snakemake profile with `export SNAKEMAKE_PROFILE=/path/to/your/profile`.  You can alternatively use the snakemake cli --profile arg

## Description
This is the workflow used for routine analysis of 16S microbiome sequencing.  Since "routine" is arguably a misnomer, it is a modularized worklfow designed to executed in a few ways.  It consists of 4 main stages:

1. demux:
   - fast demultiplexing via scattered shards using seqkit and cutadapt
2. preprocess:
   - fastqc of pools and samples
   - primer removal
   - quality trimming and filtering with DADA2
   - report generation
3. denoise:
   - calling ASVs with DADA2
4. annotate:
   - taxonomic annotation
   - adapter and host contamination checking

<!-- There is an additional workflows for -->

<!-- - uploading samples to REDcap -->
<!-- - depositing sequences in a database -->

<!-- but due to connectivity issues that can crash the pipeline these are run separate from the main pipeline. -->

## Usage
> **Note**
> All paths here are written for running the commands from the root of this repository.  Where you see relative paths, they are relative to the --directory path.

### Demultiplexing
``` sh
snakemake \
  --directory work \
  --config \
	stage=demux \
    pool=test \
    nshards=2 \
    R1=["$PWD/.test/amplicon/test_input/test_R1_001.fastq.gz"] \
    R2=["$PWD/.test/amplicon/test_input/test_R2_001.fastq.gz"] \
    oligos=$PWD/.test/amplicon/test_input/pool1059.oligos
```

### Preprocessing
``` sh
snakemake \
  --directory work \
  --config \
    stage=preprocess \
    pool=test \
    manifest=demux/test_manifest.tsv \
    primer_F=AYTGGGYDTAAAGNG \
    primer_R=CCGTCAATTYHTTTRAGT
```

### Denoising with DADA2
``` sh
snakemake \
  --directory work \
  --config \
    pool=test \
    manifest=preprocess/test_manifest.tsv
```

### Annotating ASVs with BLAST
``` sh
snakemake \
  --directory work \
  --config \
    pool=test \
	stage=annotate \
    asv_fasta=denoise/test_asvs.fasta \
    adapter_ref=/data/brinkvd/resources/references/synthetic/stephenturner-adapters/93b5f91/adapters_combined_256_unique.fasta \
    human_ref=/data/brinkvd/resources/indexes/human/hg38/hg38/minimap2/hg38.fa.masked.gz.mmi \
    mouse_ref=/data/brinkvd/resources/indexes/mouse/GRCm39/GCA_000001635.9/minimap2/GCA_000001635.9_GRCm39_genomic.fna.gz.mmi \
    anno_db=/data/brinkvd/resources/dbs/ncbi16S/2022/16S_ribosomal_RNA_id_and_taxonomy.txt
```


### Configuration

>we are no longer doing config schema validation becuase we rely on snakedeploy for typical usage.  Snakedeploy only fetches files called config.yaml, not the whole directory of configs and schemas.


Run-specific parameters are set using the `--config`.


<!-- ## REDCap Integration -->
<!-- All samples require a bit of manual review prior to inclusion in wider studies. Additionally, it would be bad to include low-quality samples when running DADA2 in `pooled` mode: hence the two-step execution outlined above.  For production purposes, we run samples without pooling.  At the end of the workflow's execution, samples and the results of the autoexclusion script are uploaded to redcap using credentials found in environment variables.  Users are encouraged to  peruse the multiqc/fastqc reports and then log their review with REDcap.  The redcap review status is a column in the samples table -->

<!-- To upload the sequences to prototype SQLite3 database, run the following command specifying the input directory (the output of a pipeline run), the pool_name, and the database location.  Details about the database are in the scheme file: `scripts/db_admin/db_schema.sql`.  A test db can be created by running `sqlite3 tmp.db < vdb_16S/scripts/db_admin/db_schema.sql`. -->
<!-- ```sh -->
<!-- snakemake \ -->
<!--   --profile ${PWD}/msk-lsf/ \ -->
<!--   --snakefile vdb_16S/rules/add_to_db.smk \ -->
<!--   --config \ -->
<!--     input_directory=${PWD}/vdb_16S/test_output/ \ -->
<!--     pool_name=test_output \ -->
<!--     db=/data/brinkvd/watersn/amplicon.db -->
<!-- # note you can add redcap_action="push/pull" as a config arg to add samples to redcap if the original pipeline run did not upload it for some reason -->
<!-- ``` -->

### Testing

#### Executing on Test Data
The decontamination step requires at least 18GB RAM to load the host minimap index; so, for local execution, use the `skip_contam_search=True` flag in the config

#### test output data
To make development easier, we add a compressed output directory to this repo after running essentially the above command:
```sh
# remove the fastqc zip to save a bit of space
rm vdb_16S/test_output/reports/fastqc/*.zip
tar cf - vdb_16S/test_output/ | xz -9e - > vdb_16S/test_output.tar.xz
# decompress
tar -xf vdb_16S/test_output.tar.xz
```

<!-- #### Generating Unit Tests -->
<!-- Snakemake has the ability to generate unittests based on a set of inputs. After running the pipeline command above, rerun with the `--generate-unit-tests` flag: -->

<!-- ```sh -->
<!-- snakemake \ -->
<!--   --profile ${PWD}/msk-lsf/ \ -->
<!--   --snakefile vdb_16S/Snakefile \ -->
<!--   --directory vdb_16S/test_output/ \ -->
<!--   --notemp \ -->
<!--   --generate-unit-tests \ -->
<!--   --config \ -->
<!--     input_directory=${PWD}/vdb_16S/test_input/ \ -->
<!--     nshards=1 \ -->
<!--     db_directory=${PWD}/vdb_16S/test_resources/ -->
<!-- ``` -->

<!-- this is still a work in progress -->
