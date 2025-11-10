#!/bin/bash

# Inputs
R1_in=$1
R2_in=$2
primer_config=$3
fq_out=$4
cutadapt_json=$5

# Populates primer sequence variables
source $primer_config

cutadapt -g ^$FWDPRIMER...${RCREVPRIMER}N{$umi_len} \
	-G ^N{$umi_len}$REVPRIMER...$RCFWDPRIMER \
	--discard-untrimmed \
	--no-indels \
	--action retain \
	--interleaved \
	--quiet \
	--json $cutadapt_json \
	-o $fq_out $R1_in $R2_in

