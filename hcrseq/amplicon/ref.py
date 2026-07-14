import subprocess

from Bio import SeqIO
import pysam

import pandas as pd

from hcrseq.common.util import rc

def index_reference(ref):
    subprocess.check_output(f"bwa index {ref}",shell=True)
