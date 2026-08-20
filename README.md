# HCR-seq

This repo contains code for processing and analyzing HCR-seq data

## Installation

To install this package and its prerequisites
```
conda install -f environment.yaml
pip install -e .
```

## Preparing an HCR-seq reference

Both the amplicon and single-cell workflows are built on top of a shared `Reference` object
(pickled to a `.pkl` file) that describes each reporter plasmid, its associated lesions, and how
repair pathways should be calculated from those reporters.

The latest version of the reference is [here](reference/HCRseq_v0.4.1).

It is built with a small script (see [build.py](reference/HCRseq_v0.4.1/build.py)) that calls
`Reference.build`:

```python
from hcrseq.common.ref import Reference

r = Reference.build('plasmid_info.txt', 'reporter_info.txt', 'pathway_calculations.txt')
r.write_pickle('HCRseq_v0.4.1.pkl')
```

The three inputs are:

#### Plasmid info file

Specifies the sequence of each reporter plasmid, along with the position at which its barcode
(the HR donor sequence) is inserted.

#### Reporter info file

Specifies, for each named reporter, which plasmid it derives from, its barcode sequence, and the
position/base-change of the lesion it carries.

#### Pathway calculation file

Specifies how each repair pathway should be calculated from the reporters above (e.g. which
reporter/lesion/metric combination corresponds to which pathway).

This produces a single `HCRseq_v0.4.1.pkl` file, which is the `--hcrseq_ref`/`--ref_path` input
used throughout the amplicon and single-cell commands below.

## HCR-seq Amplicon analysis

### Preparing an amplicon reference

Given the pickled `Reference` above, this derives a primer-trimmed amplicon reference: an indexed
fasta for alignment, plus a matching `Reference` pickle (restricted to reporters containing the
given primers) for quantification.

```
hcrseq amplicon prepare-amplicon-reference --hcrseq_ref=HCRseq_v0.4.1.pkl \
                                            --forward_primer=ACAACCACTACCTGAG \
                                            --reverse_primer=TCACTTGTACAGCTCGTCCATGC \
                                            --umi_len=15 \
                                            --outstem=amplicon/HCRseq_v0.4.1
```

The prebuilt amplicon reference is [here](reference/HCRseq_v0.4.1/amplicon). Once run, the
following files will be prepared.

#### Reference fasta

This is a fasta file containing the sequence of all amplicon sequences. It will also have
associated index files (`.amb`/`.ann`/`.bwt`/`.pac`/`.sa`/`.fai`/`.dict`).

#### Reference pickle

A copy of the `Reference` pickle restricted to the reporters covered by the given primers. This is
the `--ref_path` used by `hcrseq amplicon quantify`.

#### Primer file

This file specifies the primer sequences and UMI length used. These are used to identify read
pairs conforming to the expected amplicon structure and for adapter trimming.

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
                --ref_path=${ref_path} \
                --outstem=${outstem} \
                --min_mapq=5
```

### Sequencing saturation

Estimates a sequencing saturation curve from the UMI deduplication metrics produced during
preprocessing.

```
hcrseq amplicon calculate-sequencing-saturation --bam=${outstem}.bam \
                --cutadapt_json=${outstem}.cutadapt_metrics.json \
                --outstem=${outstem}
```

### Aggregating multiple samples

Combines per-sample cutadapt/count/repair outputs across a set of sample ids into a single set of
aggregated tables.

```
hcrseq amplicon aggregate --cutadapt_files=${cutadapt_json_1},${cutadapt_json_2},... \
                --count_files=${counts_1},${counts_2},... \
                --repair_files=${repair_1},${repair_2},... \
                --ids=${id_1},${id_2},... \
                --outstem=${outstem}
```

### Running in Terra

See the following workspace as an example of running preprocessing in Terra
https://app.terra.bio/#workspaces/broad-getzlab-fmhcrsparc-terra/Nagel-FM-HCR-Amplicon-experiment7

## Single-cell processing

### Reference preparation

First prepare a cellranger-compatible reference including reporter plasmids as follows

```
hcrseq scrna prepare-reference --hcrseq_ref=[hcrseq_ref] \
                    --genome_fasta=[genome_fasta] \
                    --gtf=[gtf] \
                    --outstem=[outstem]
```

### Alignment

The data can then be processed using cellranger count. Afterwards, reads aligning to the reporter
contigs are extracted, realigned, and sorted/indexed using hcrseq as follows.

```
hcrseq scrna postprocess-reporter-bam --bam=[bam] \
          --ref_path=[ref_path] \
          --outstem=[outstem] \
          --threads=8 \
          --mem-per-thread=4G
```

This produces `[outstem]_postprocessed.sorted.bam`, which is the bam passed to `quantify` below.
`--threads`/`--mem-per-thread` control the CPU threads and per-thread memory used when sorting and
indexing the resulting bam.

### Quantification

```
hcrseq scrna quantify --h5_file=[h5_file] \
          --bam=[outstem]_postprocessed.sorted.bam \
          --ref_path=[ref_path] \
          --outstem=[outstem] \
          --min_mapq=5
```
This will produce a scanpy h5ad file with repair measurements annotated in adata.obs. 
