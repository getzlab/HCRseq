import pysam
import sys
from collections import Counter
from umi_tools import UMIClusterer
import tempfile

HEADER={ 'HD': {'VN': '1.0'}}

def extract_umi_to_ubam(fastq_in,
                bam_out,
                umi_len=15):
    """
    Extracts UMI from 3' end of read
    """

    # Save UMI counts for later error correction
    umi_counts = Counter()

    with pysam.FastxFile(fastq_in) as fh, \
        pysam.AlignmentFile(bam_out, "wb", header= HEADER) as outf:
        for read in fh:

            read_len = len(read.sequence) - umi_len
            umi = read.sequence[(-1*umi_len):]
            umi_counts[umi.encode('utf-8')] += 1

            a = pysam.AlignedSegment()
            a.flag = 4
            a.query_name = read.name
            a.query_sequence = read.sequence[0:read_len]
            a.query_qualities = pysam.qualitystring_to_array(read.quality[0:read_len])
            a.tags = (("UR", umi),
                      ("UY", read.quality[(-1*umi_len):]))

            outf.write(a)

    return(umi_counts)

def correct_umis(bam_in,bam_out,umi_counts,threshold=1):
    """
    Counts uncorrected UMIs
    """

    # Use UMI-tools to create mapping of corrected barcodes
    clusterer = UMIClusterer(cluster_method="directional")
    clustered_umis = clusterer(umi_counts, threshold=threshold)
    umi_map = {e.decode("utf-8") : entry[0].decode("utf-8") for entry in clustered_umis for e in entry}

    # Map tags back into bam
    with pysam.AlignmentFile(bam_in, "rb",check_sq=False) as fh_in, \
        pysam.AlignmentFile(bam_out,"wb", template=fh_in) as fh_out:
        for read in fh_in:
            umi = read.get_tag("UR")
            ub = umi_map[umi]
            read.set_tag("UB",ub)
            fh_out.write(read)

if __name__=="__main__":
    fq_in = sys.argv[1]
    bam_out = sys.argv[2]
    use_umi_correction = True

    if use_umi_correction:
        with tempfile.NamedTemporaryFile(mode="w",suffix=".bam") as tmpbam:
            umi_counts = extract_umi_to_ubam(fq_in, tmpbam.name)
            correct_umis(tmpbam.name, bam_out,umi_counts)
    else:
        extract_umi_to_ubam(fq_in,bam_out)