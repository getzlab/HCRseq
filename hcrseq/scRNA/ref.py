from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess
from hcrseq.common.ref import Reference


def merge_genome_and_plasmid_ref(hg_file,hcrseq_ref,hg_gtf,out_fasta_file,out_gtf_file):

    ref = Reference.load(hcrseq_ref)

    ## Start by copying/downloading human genome fasta and gtf
    def create_output_file(input_file,output_file):
        if input_file.startswith('gs://'):
            cmd = f'gsutil cp {input_file} {output_file}'
        elif input_file.startswith('ftp://'):
            cmd = f'curl {input_file} > {output_file}'
        else:
            cmd = f'cp {input_file} {output_file}'
        subprocess.check_output(cmd,shell=True)

    create_output_file(hg_file,out_fasta_file)
    create_output_file(hg_gtf,out_gtf_file)

    ## Then append the plasmid sequences
    with open(out_fasta_file, "at") as outfasta, open(out_gtf_file, "at") as outgtf:
        for reporter in ref.reporters:

            SeqIO.write(reporter.to_seqrecord(), outfasta, "fasta")

            outgtf.write(f'{reporter.name}\tCUSTOM\tgene\t{reporter.tx_start}\t{reporter.tx_end}\t.\t+\t.\tgene_id "{reporter.name}";\n')
            outgtf.write(f'{reporter.name}\tCUSTOM\ttranscript\t{reporter.tx_start}\t{reporter.tx_end}\t.\t+\t.\tgene_id "{reporter.name}"; transcript_id "{reporter.name}_t";\n')
            outgtf.write(f'{reporter.name}\tCUSTOM\texon\t{reporter.tx_start}\t{reporter.tx_end}\t.\t+\t.\tgene_id "{reporter.name}"; gene_type "protein_coding"; transcript_id "{reporter.name}_t"; exon_id "{reporter.name}_e";\n')
