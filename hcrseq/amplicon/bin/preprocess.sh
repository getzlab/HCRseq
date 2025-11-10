#!/bin/bash

R1=$1
R2=$2
ref=$3
primer_config=$4
outstem=$5

sdir=$(dirname "${BASH_SOURCE[0]}")

$sdir/trim_amplicons.sh ${R1} ${R2} $primer_config /dev/stdout $outstem.cutadapt_metrics.json | \
  python $sdir/../paired_consensus.py /dev/stdin /dev/stdout | \
  python $sdir/../extract_umi.py /dev/stdin /dev/stdout | \
  samtools sort -t UR | tee $outstem.tagged.ubam | \
  python $sdir/../dedup.py /dev/stdin /dev/stdout | \
  picard SortSam \
    INPUT=/dev/stdin \
    OUTPUT=$outstem.dedup.ubam \
    SORT_ORDER=queryname \
    USE_JDK_DEFLATER=true \
    USE_JDK_INFLATER=true

  samtools fastq $outstem.dedup.ubam | \
  bwa mem -k10 $ref /dev/stdin | \
  python $sdir/../global_realignment.py /dev/stdin /dev/stdout $ref | \
  picard SortSam \
    INPUT=/dev/stdin \
    OUTPUT=$outstem.aligned.bam \
    SORT_ORDER=queryname \
    USE_JDK_DEFLATER=true \
    USE_JDK_INFLATER=true

  picard MergeBamAlignment \
    REFERENCE_SEQUENCE=$ref \
    UNMAPPED_BAM=$outstem.dedup.ubam \
    ALIGNED_BAM=$outstem.aligned.bam \
    OUTPUT=$outstem.bam \
    CREATE_INDEX=true \
    SORT_ORDER=coordinate \
    USE_JDK_DEFLATER=true \
    USE_JDK_INFLATER=true