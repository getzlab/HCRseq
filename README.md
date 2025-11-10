# HCR-seq

This repo contains code for processing and analyzing HCR-seq data

## Installation

To install this package and its prerequisites
```
conda install -f environment.yaml
pip install -e .
```

## HCR-seq Amplicon analysis

### Preparing a reference

The latest version of the reference is [here](reference/reporters/HCRseq_v0.1/amplicon).

This reference was prepared with the following command:

```
hcrseq amplicon prepare-reference --reporter_info=reporter_info.txt \
                                  --mmej_variants=MMEJ_variants.txt \
                                  --forward_primer=ACAACCACTACCTGAG \
                                  --reverse_primer=TCACTTGTACAGCTCGTCCATGC \
                                  --outstem=amplicon/HCRseq_v0.1
```
The files used to build it are [here](reference/reporters/HCRseq_v0.1). The two main inputs are reporter_info.txt which contains information about the sequence and lesion positions of each reporter, and MMEJ_variants.txt which gives expected MMEJ deletion variants to be added to the reference. Once run, the following files will be prepared.

#### Reference fasta

This is a fasta file containing the sequence of all amplicon sequences, including possible MMEJ variants. It will also have associated index files.

#### Primer file

This file specifies the primer sequences used. These are used to identify read pairs conforming to the expected amplicon structure and for adapter trimming. 

#### Lesion info file

This file specifies site specific lesions that should be quantified. 

#### Pathway calculation file

This file specifies how each repair pathway should be calculated.

### Preprocessing

This step performs the following operations
- Creates a consensus between read 1 and read 2
- Extracts the UMI sequence and deduplicates
- Aligns amplicon sequences to reference
  
```
hcrseq amplicon preprocess --fq1=${fq1} \
      --fq2=${fq2} \
      --reference=${ref} \
      --primer_config=${primer_file} \
      --outstem=${outstem}
```

### Quantification

Given an aligned bam, this step counts reporters and calculates pathway activities

```
hcrseq amplicon quantify --bam=${outstem}.bam \
                --lesion_info=${lesion_info} \
                --pathway_info=${pathway_info} \
                --outstem=${outstem}
```

### Running in Terra

See the following workspace as an example of running preprocessing in Terra
https://app.terra.bio/#workspaces/broad-getzlab-fmhcrsparc-terra/Nagel-FM-HCR-Amplicon-experiment7

## Single-cell processing

### Reference preparation

First prepare a cellranger-compatible reference including reporter plasmids as  follows
```
hcrseq scrna prepare-reference --plasmid_fasta=[plasmid_fasta] \
                    --genome_fasta=[genome_fasta] \
                    --gtf=[gtf] \
                    --outstem=[outstem]
```

### Alignment

The data can the be processed using cellranger count. Afterwards, data can be post-processed using hcrseq as follows.

```
hcrseq scrna quantify --h5_file=[h5_file] \
          --bam=[bam] \
          --lesion_info=[lesion_info] \
          --pathway_info=[pathway_info] \
          --outstem=[outstem]
```
This will produce a scanpy h5ad file with repair measurements annotated in adata.obs. 

