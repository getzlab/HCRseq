#!/bin/bash

hcrseq scrna prepare-reference \
	--hcrseq_ref HCRseq_v0.4.pkl \
	--genome_fasta /home/njharlen/refdata-gex-GRCh38-2024-A/fasta/genome.fa \
	--gtf /home/njharlen/cellranger_genes.gtf \
	--outstem HCRseq_v0_4_scrna
