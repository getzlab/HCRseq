from Bio import SeqIO
import subprocess

def merge_genome_and_plasmid_ref(hg_file,reporter_fasta,hg_gtf,out_fasta_file,out_gtf_file):

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
    with open(out_fasta_file, "at") as outfasta, open(out_gtf_file, "at") as outgtf, open(reporter_fasta,"rt") as infasta:
        for record in SeqIO.parse(infasta,"fasta"):
            reporter_name = record.id

            SeqIO.write(record, outfasta, "fasta")

            outgtf.write(f'{reporter_name}\tCUSTOM\tgene\t1\t{len(record)}\t.\t+\t.\tgene_id "{reporter_name}";\n')
            outgtf.write(f'{reporter_name}\tCUSTOM\ttranscript\t1\t{len(record)}\t.\t+\t.\tgene_id "{reporter_name}"; transcript_id "{reporter_name}";\n')
            outgtf.write(f'{reporter_name}\tCUSTOM\texon\t1\t{len(record)}\t.\t+\t.\tgene_id "{reporter_name}"; gene_type "protein_coding"; transcript_id "{reporter_name}"; exon_id "{reporter_name}";\n')
