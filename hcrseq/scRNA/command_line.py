import subprocess

import click
import pysam
from hcrseq.common.ref import Reference
from hcrseq.scRNA.ref import merge_genome_and_plasmid_ref
from hcrseq.scRNA.postprocess import realign_reporter_transcripts
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
@click.option("--bam")
@click.option("--ref_path")
@click.option("--outstem")
def postprocess_reporter_bam(bam,ref_path,outstem):
    """
    Extracts reads aligning to the reporter contigs, realigns reporter transcripts,
    then sorts and indexes the resulting bam
    """
    ref = Reference.load(ref_path)
    reporter_contigs = [reporter.name for reporter in ref.reporters]

    reporter_bam = outstem + ".reporters.bam"
    realigned_bam = outstem + ".realigned.bam"
    sorted_bam = outstem + "_postprocessed.sorted.bam"

    # Extract reads aligning to the reporter contigs
    pysam.view("-b","-o",reporter_bam,bam,*reporter_contigs,catch_stdout=False)

    # Realign reporter transcripts
    realign_reporter_transcripts(reporter_bam,realigned_bam,ref_path)

    # Sort and index the final bam
    pysam.sort("-o",sorted_bam,realigned_bam)
    pysam.index(sorted_bam)

@scrna.command()
@click.option("--h5_file")
@click.option("--bam")
@click.option("--ref_path")
@click.option("--outstem")
@click.option("--min_mapq",default=5,type=int)
def quantify(h5_file,bam,ref_path,outstem,min_mapq):
    """
    Takes filtered h5 file from cellranger and extracts repair metrics from bam
    outputs an h5ad file
    """
    adata = quantify_repair(h5_file,bam,
                    ref_path,
                    min_mapq=min_mapq)
    adata.write(outstem + ".h5ad")