import subprocess

import pysam

from hcrseq.amplicon.quantify import count_umis, quantify_repair
from hcrseq.amplicon.ref import index_reference
from hcrseq.amplicon.aggregate import aggregate_results
from hcrseq.amplicon.downsampling import calculate_sequencing_saturation_curve
from hcrseq.common.ref import Reference
from hcrseq.common.util import rc
import yaml
import glob
import click
import os
import pandas as pd

@click.group()
def amplicon():
    pass

@amplicon.command()
@click.option("--fq1")
@click.option("--fq2")
@click.option("--reference")
@click.option("--primer_config")
@click.option("--outstem")
def preprocess(fq1,fq2,reference,primer_config,outstem):
    """
    Runs preprocessing of hcr-seq amplicon data
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    subprocess.check_output(f'{dir_path}/bin/preprocess.sh {fq1} {fq2} {reference} {primer_config} {outstem}',shell=True)

@amplicon.command()
@click.option("--bam")
@click.option("--cutadapt_json")
@click.option("--outstem")
def calculate_sequencing_saturation(bam,cutadapt_json,outstem):
    """
    Runs preprocessing of hcr-seq amplicon data
    """
    calculate_sequencing_saturation_curve(bam,cutadapt_json,outstem)

@amplicon.command()
@click.option("--bam")
@click.option("--ref_path")
@click.option("--outstem")
@click.option("--min_mapq",default=5,type=int)
def quantify(bam,ref_path,outstem,min_mapq):
    """
    Runs repair quantification on bam file
    """

    print('QUANTIFYING...')
    reference = Reference.load(ref_path)
    counts,del_df,ins_df,mismatch_dist = count_umis(bam,reference,min_mapq=min_mapq)

    with open(f"{outstem}.reporter_metrics.yaml","wt") as f:
        yaml.dump({name: dict(reporter_counts) for name, reporter_counts in counts.items()},f)

    mismatch_df = pd.DataFrame(mismatch_dist).T.explode(['nsnp','ndel','nins','N']).reset_index(names='contig')
    mismatch_df['position'] =  mismatch_df.groupby('contig').cumcount() + 1

    mismatch_df.to_csv(f"{outstem}.mismatch_dist.csv",index=None)

    del_df.to_csv(f"{outstem}.deletions.csv",index=None)
    ins_df.to_csv(f"{outstem}.insertions.csv",index=None)

    q = quantify_repair(counts,reference)
    with open(f"{outstem}.repair_measurements.yaml","wt") as f:
        yaml.dump(q,f)

@amplicon.command()
@click.option("--cutadapt_files")
@click.option("--count_files")
@click.option("--repair_files")
@click.option("--ids")
@click.option("--outstem")
def aggregate(cutadapt_files,count_files,repair_files,ids,outstem):
    aggregate_results(cutadapt_files.split(","),
                      count_files.split(","),
                      repair_files.split(","),
                      ids.split(","),
                      outstem)

@amplicon.command()
@click.option("--hcrseq_ref")
@click.option("--forward_primer")
@click.option("--reverse_primer")
@click.option("--umi_len",default=15,type=int)
@click.option("--outstem")
def prepare_amplicon_reference(hcrseq_ref,forward_primer,reverse_primer,umi_len,outstem):
    """
    Given a pickled HCR-seq Reference (eg. one built for scRNA), derives a primer-trimmed
    amplicon reference: an indexed fasta for alignment plus a matching Reference pickle
    (restricted to reporters containing the given primers) for quantification.
    """
    reference = Reference.load(hcrseq_ref)
    amplicon_reference = reference.to_amplicon((forward_primer,reverse_primer))

    out_fasta = outstem + '.fasta'
    amplicon_reference.write_amplicon_fasta(out_fasta)
    index_reference(out_fasta)
    pysam.dict(out_fasta, '-o', outstem + '.dict')

    amplicon_reference.write_pickle(outstem + '.pkl')

    primer_content = \
        f"""FWDPRIMER={forward_primer}
RCFWDPRIMER={rc(forward_primer)}
RCREVPRIMER={rc(reverse_primer)}
REVPRIMER={reverse_primer}
umi_len={umi_len}
"""
    with open(outstem + '.primers.sh','wt') as pout:
        pout.write(primer_content)




