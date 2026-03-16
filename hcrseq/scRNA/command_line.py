import subprocess

import click
from hcrseq.scRNA.ref import merge_genome_and_plasmid_ref
from hcrseq.scRNA.quantify import quantify_repair

@click.group()
def scrna():
    pass

@scrna.command()
@click.option("--hcrseq_ref")
@click.option("--genome_fasta")
@click.option("--gtf")
@click.option("--outstem")
def prepare_reference(hcrseq_ref,genome_fasta,gtf,outstem):
    """
    Adds plasmid sequences to reference fasta and gtf, then runs cellranger mkref
    """

    merge_genome_and_plasmid_ref(
        hg_file=genome_fasta,
        hg_gtf=gtf,
        hcrseq_ref=hcrseq_ref,
        out_fasta_file=outstem + '.fasta',
        out_gtf_file=outstem + '.gtf'
    )

    subprocess.check_output(f"cellranger mkref --genome={outstem} --fasta={outstem}.fasta --genes={outstem}.gtf --nthreads=16",shell=True)

@scrna.command()
@click.option("--h5_file")
@click.option("--bam")
@click.option("--lesion_info")
@click.option("--pathway_info")
@click.option("--outstem")
def quantify(h5_file,bam,lesion_info,pathway_info,outstem):
    """
    Takes filtered h5 file from cellranger and extracts repair metrics from bam
    outputs an h5ad file
    """
    adata = quantify_repair(h5_file,bam,
                    lesion_info,
                    pathway_info)
    adata.write(outstem + ".h5ad")