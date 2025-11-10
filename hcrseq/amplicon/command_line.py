import subprocess

from hcrseq.amplicon.quantify import count_umis, quantify_repair
from hcrseq.amplicon.ref import create_amplicon_reference_v2
from hcrseq.amplicon.aggregate import aggregate_results
from hcrseq.amplicon.downsampling import calculate_sequencing_saturation_curve
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
@click.option("--reference")
@click.option("--lesion_info")
@click.option("--pathway_info")
@click.option("--outstem")
def quantify(bam,reference,lesion_info,pathway_info,outstem):
    """
    Runs repair quantification on bam file
    """

    print('QUANTIFYING...')
    counts,del_df,ins_df,mismatch_dist = count_umis(bam,reference,lesion_info)

    with open(f"{outstem}.reporter_metrics.yaml","wt") as f:
        yaml.dump(counts,f)

    mismatch_df = pd.DataFrame(mismatch_dist).T.explode(['nsnp','ndel','nins','N']).reset_index(names='contig')
    mismatch_df['position'] =  mismatch_df.groupby('contig').cumcount() + 1

    mismatch_df.to_csv(f"{outstem}.mismatch_dist.csv",index=None)

    del_df.to_csv(f"{outstem}.deletions.csv",index=None)
    ins_df.to_csv(f"{outstem}.insertions.csv",index=None)

    q = quantify_repair(counts,pathway_info)
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
@click.option("--reporter_info")
@click.option(("--mmej_variants"))
@click.option("--forward_primer")
@click.option("--reverse_primer")
@click.option("--outstem")
def prepare_reference(reporter_info,mmej_variants,forward_primer,reverse_primer,outstem):
    """
    Given TSV with reporter info, creates fasta and indices, lesion position, and primer files need to run reprocessing
    """
    create_amplicon_reference_v2(reporter_info=reporter_info,
                                 mmej_variants=mmej_variants,
                                 forward_primer=forward_primer,
                                 reverse_primer=reverse_primer,
                                 outstem=outstem)




