#! /bin/bash
set -eux

mode=$1

R1="[$PWD/.test/amplicon/test_input/test_R1_001.fastq.gz]"
R2="[$PWD/.test/amplicon/test_input/test_R2_001.fastq.gz]"
commonargs=" --restart-times 0 -j 1"
echo $commonargs
case $mode in

  full | all)
      snakemake \
	  --directory tmpall/ \
	  --singularity-args "-B ${PWD},/data/brinkvd/,/lila/$PWD/,/lila/data/brinkvd/,/scratch/" \
	  $commonargs \
	  --config \
	  pool="pool1" \
	  stage="all" \
	  sample=473 \
	  oligos=$PWD/.test/amplicon/test_input/pool1059.oligos \
	  R1=$R1 \
	  R2=$R2 \
	  nshards=1
      ;;
  demux )
      snakemake \
	  --directory tmpdemux/ \
	  --singularity-args "-B ${PWD},/data/brinkvd/,/lila/$PWD/,/lila/data/brinkvd/,/scratch/" \
	  $commonargs \
	  --config \
	  stage=demux \
	  sample=473 \
	  pool="pool1" \
	  manifest="$PWD/.test/test_manifest.tsv" \
	  oligos=$PWD/.test/amplicon/test_input/pool1059.oligos \
	  R1=$R1 \
	  R2=$R2 \
	  nshards=1
      ;;
  preprocess )
      snakemake \
	  --singularity-args "-B ${PWD},/data/brinkvd/,/lila/$PWD/,/lila/data/brinkvd/,/scratch/" \
	  $commonargs \
	  --directory tmppre/   \
	  --config \
	  stage=preprocess \
	  sample=473  \
	  pool="pool1" \
	  oligos=$PWD/.test/amplicon/test_input/pool1059.oligos \
	  manifest=$PWD/tmpdemux/demux/pool1_manifest.tsv \
	  R1=$R1 \
	  R2=$R2 \
	  primer_F=AYTGGGYDTAAAGNG \
	  primer_R=CCGTCAATTYHTTTRAGT \
	  nshards=2
      ;;
  preprocess_se )
      snakemake \
	  --singularity-args "-B ${PWD},/data/brinkvd/,/lila/$PWD/,/lila/data/brinkvd/,/scratch/" \
	  $commonargs \
	  --directory tmpprese/   \
	  --config \
	  stage=preprocess \
	  sample=473  \
	  pool="pool1" \
	  oligos=$PWD/.test/amplicon/test_input/pool1059.oligos \
	  manifest=$PWD/.test/manifest_se.tsv \
	  primer_F=AYTGGGYDTAAAGNG \
	  primer_R=CCGTCAATTYHTTTRAGT \
	  lib_layout=single \
	  nshards=2
      ;;
  preprocess_se_noprimers )
      snakemake \
	  --singularity-args "-B ${PWD},/data/brinkvd/,/lila/$PWD/,/lila/data/brinkvd/,/scratch/" \
	  $commonargs \
	  --directory tmpprese_np/   \
	  --config \
	  stage=preprocess \
	  sample=473  \
	  pool="pool1" \
	  oligos=$PWD/.test/amplicon/test_input/pool1059.oligos \
	  manifest=$PWD/.test/manifest_se.tsv \
	  primer_F=AYTGGGYDTAAAGNG \
	  primer_R=CCGTCAATTYHTTTRAGT \
	  lib_layout=single \
	  removeprimers=False \
	  nshards=2
      ;;
  denoise_pooled )
      snakemake \
	  $commonargs \
	  --singularity-args "-B ${PWD},/data/brinkvd/,/lila/$PWD/,/lila/data/brinkvd/,/scratch/" \
	  --directory tmpdenoise/ \
	  --config \
	  stage=denoise \
	  sample=473  \
	  pool="pool1" \
	  manifest=$PWD/tmppre/preprocess/pool1_manifest.tsv \
	  oligos=$PWD/.test/amplicon/test_input/pool1059.oligos \
	  R1=$R1 \
	  R2=$R2 \
	  nshards=2
      ;;
  denoise_pooled_se )
      snakemake \
	  $commonargs \
	  --singularity-args "-B ${PWD},/data/brinkvd/,/lila/$PWD/,/lila/data/brinkvd/,/scratch/" \
	  --directory tmpdenoisese/ \
	  --config \
	  stage=denoise \
	  sample=473  \
	  pool="pool1" \
	  manifest=$PWD/tmppre/preprocess/pool1_manifest.tsv \
	  oligos=$PWD/.test/amplicon/test_input/pool1059.oligos \
	  R1=$R1 \
	  R2=$R2 \
	  lib_layout="single" \
	  pooling="pseudo" \
	  nshards=2
      ;;
  denoise_pseudo )
      snakemake \
	  $commonargs \
	  --singularity-args "-B ${PWD},/data/brinkvd/,/lila/$PWD/,/lila/data/brinkvd/,/scratch/" \
	  --directory tmpdenoise_pseudo/ \
	  --config \
	  stage=denoise \
	  sample=473  \
	  pool="pool1" \
	  pooling="pseudo" \
	  manifest=$PWD/tmppre/preprocess/pool1_manifest.tsv \
	  oligos=$PWD/.test/amplicon/test_input/pool1059.oligos \
	  R1=$R1 \
	  R2=$R2 \
	  nshards=2
      ;;
  denoise )
      snakemake \
	  $commonargs \
	  --singularity-args "-B ${PWD},/data/brinkvd/,/lila/$PWD/,/lila/data/brinkvd/,/scratch/" \
	  --directory tmpdenoise_unpooled/ \
	  --config \
	  stage=denoise \
	  sample=473  \
	  pool="pool1" \
	  pooling=none \
	  manifest=$PWD/tmppre/preprocess/pool1_manifest.tsv \
	  oligos=$PWD/.test/amplicon/test_input/pool1059.oligos \
	  R1=$R1 \
	  R2=$R2 \
	  nshards=2
      ;;
  denoise_se )
      snakemake \
	  $commonargs \
	  --singularity-args "-B ${PWD},/data/brinkvd/,/lila/$PWD/,/lila/data/brinkvd/,/scratch/" \
	  --directory tmpdenoisese_unpooled/ \
	  --config \
	  stage=denoise \
	  sample=473  \
	  pool="pool1" \
	  pooling="none" \
	  manifest=$PWD/tmpprese/preprocess/pool1_manifest.tsv \
	  oligos=$PWD/.test/amplicon/test_input/pool1059.oligos \
	  R1=$R1 \
	  lib_layout="single" \
	  nshards=2
      ;;
  annotate )
      snakemake \
	  $commonargs \
	  --singularity-args "-B ${PWD},/data/brinkvd/,/lila/$PWD/,/lila/data/brinkvd/,/scratch/" \
	  --directory tmpannotate/ \
	  --config \
	  stage=annotate \
	  asv_fasta="${PWD}/tmpdenoise/denoise/pool1_asvs.fasta" \
	  sample=473  \
	  pool="pool1" \
	  pooling="none" \
	  manifest=$PWD/tmppre/preprocess/pool1_manifest.tsv \
	  oligos=$PWD/.test/amplicon/test_input/pool1059.oligos \
	  R1=$R1 \
	  R2=$R2 \
	  nshards=2
      ;;

  *)
    echo -e "unknown mode; please chose from annotate, denoise, denoise_pseudo, demux, preprocess. Exiting\n"
    ;;
esac
