import time
import pysam
from hcrseq.common.quantify import HCRseqQuantifier
from hcrseq.common.ref import Reference

if __name__ == '__main__':
    bam_in = "scRNA_test_results/experiment11_postprocessed.sorted.bam"
    ref_path = "../reference/HCRseq_v0.4.1/HCRseq_v0.4.1.pkl"

    ref = Reference.load(ref_path)

    q = HCRseqQuantifier(ref,("CB","UB"))

    q.count_umis(bam_in)
    print(q.get_counts())

