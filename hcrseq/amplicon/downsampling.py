import matplotlib.pyplot as plt
import numpy as np
import json
from scipy.stats import hypergeom
import pysam
import seaborn as sns
import pandas as pd

def calculate_sequencing_saturation_curve(bam ,cutadapt_file,outstem):

    # Determine total number of reads
    with open(cutadapt_file ,'rt') as f:
        j = json.load(f)
    M = j['read_counts']['input']

    # Points along which to calculate expected #UMI
    N = np.linspace(0 ,M).astype(int)


    # Iterate through each UMI and calculate expectation of detection at given coverage, add to total
    x = np.zeros(len(N))
    i=0
    with pysam.AlignmentFile(bam, 'rb') as bam:
        for r in bam:
            n = r.get_tag('NR')

            x += ( 1 -hypergeom.pmf(0 ,M=M ,n=n ,N=N))

            i+=1

    df = pd.DataFrame({'Input reads': N, 'Expected UMIs recovered': x})

    sns.lineplot(x="Input reads", y="Expected UMIs recovered", data=df)
    plt.savefig(outstem + ".saturation_curve.pdf")

    df.to_csv(outstem + ".saturation_curve.csv")
