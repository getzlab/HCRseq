#!/bin/bash

genome=$1
gtf=$2
plasmid_dir=$3

add_REPORTER_to_reference.py -i $plasmid_reporters

cellranger mkref --genome=${genome} \
        --fasta=${fasta} \
        --genes=${genes} \
        --nthreads=16
