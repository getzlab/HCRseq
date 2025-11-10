import subprocess

from Bio import SeqIO
import pysam

import pandas as pd

from hcrseq.amplicon.util import rc


def create_amplicon_reference_v2(reporter_info,
                                 mmej_variants,
                                 outstem,
                                 forward_primer,
                                 reverse_primer):
    """
    Given list of reporter plasmids, creates full hcr-seq reference
    """

    out_fasta = outstem + '.fasta'
    out_dict = outstem + '.dict'
    out_lesion_file = outstem + '.lesions.txt'
    out_primer_file = outstem + '.primers.sh'

    R = pd.read_csv(reporter_info,sep='\t',index_col=None)
    M = pd.read_csv(mmej_variants,sep='\t')

    # Create entries for MMEJ variants
    for ind,row in M.iterrows():
        base_plasmid = row['reporter']
        base_sequence = R.set_index('reporter').loc[base_plasmid,'sequence']
        variant_sequence = base_sequence[0:(row['start']-1)] + base_sequence[(row['end']):]

        # Add it
        R = pd.concat([R,pd.DataFrame(pd.Series({'reporter' : row['name'],'sequence' : variant_sequence})).T])

    # Because match sequence for DSB is one less than others
    trimmed_forward_primer = forward_primer[1:]
    R['amplicon_start'] = [seq.find(trimmed_forward_primer) + len(trimmed_forward_primer) for seq in R['sequence']]
    R['amplicon_end'] = [seq.find(rc(reverse_primer)) for seq in R['sequence']]


    ## Write amplicon sequences as fasta
    with open(out_fasta,'w') as fout:
        for ind,row in R.iterrows():

            # Write it to fasta
            r = pysam.FastxRecord()
            r.name = row['reporter']
            r.sequence = forward_primer + row['sequence'][row['amplicon_start']:row['amplicon_end']] + rc(reverse_primer)

            fout.write(str(r) + '\n')

    ## Index for bwa
    index_reference(out_fasta)

    ## Dict for picard
    pysam.dict(out_fasta, '-o', out_dict)

    ## update lesions file
    L = R[~R['unrepaired_base'].isna()].copy()
    L['position'] = L['position'] - L['amplicon_start'] + len(forward_primer)
    L[['reporter','position','repaired_base','unrepaired_base']].to_csv(out_lesion_file,sep='\t',index=None)

    ## Write primer file
    primer_content = \
        f"""FWDPRIMER={forward_primer}
RCFWDPRIMER={rc(forward_primer)}
RCREVPRIMER={rc(reverse_primer)}
REVPRIMER={reverse_primer}
umi_len=15
"""
    with open(out_primer_file,'wt') as pout:
        pout.write(primer_content)


def write_ref_file(amps,outfile,format):
    with open(outfile, 'w') as handle:
        SeqIO.write(amps, handle, format=format)



def index_reference(ref):
    subprocess.check_output(f"bwa index {ref}",shell=True)
