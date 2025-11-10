import pysam
import sys


def extract_umi_to_ubam(fastq_in,
                bam_out,
                umi_len=15):
    """
    Extracts UMI from 3' end of read
    """
    with pysam.FastxFile(fastq_in) as fh, \
        pysam.AlignmentFile(bam_out, "wb", header= { 'HD': {'VN': '1.0'}}) as outf:
        for read in fh:

            read_len = len(read.sequence) - umi_len


            a = pysam.AlignedSegment()
            a.flag = 4
            a.query_name = read.name
            a.query_sequence = read.sequence[0:read_len]
            a.query_qualities = pysam.qualitystring_to_array(read.quality[0:read_len])
            a.tags = (("UR", read.sequence[(-1*umi_len):]),
                      ("UY", read.quality[(-1*umi_len):]))

            outf.write(a)

if __name__=="__main__":
    fq_in = sys.argv[1]
    bam_out = sys.argv[2]

    extract_umi_to_ubam(fq_in,bam_out)