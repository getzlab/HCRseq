import pysam
import gzip
import numpy as np
import sys

from hcrseq.common.util import rc

BP = {'A' : 'T',
      'C' : 'G',
      'G' : 'C',
      'T' : 'A',
      'N' : 'N'}


def create_consensus_from_readpair(R1, R2, offset):
    """
    Takes a single readpair and constructs a consensus read
    """

    min_len = min(len(R1.sequence), len(R2.sequence))

    # Split reads by overlap
    r_pre = R1.sequence[0:offset]
    q_pre = R1.quality[0:offset]

    r1 = np.array(list(R1.sequence[offset:(min_len)]))
    q1 = np.array(list(R1.quality[offset:(min_len)]))

    r2 = np.array(list(rc(R2.sequence)[0:(min_len - offset)]))
    q2 = np.array(list(R2.quality[::-1][0:(min_len - offset)]))

    r_post = rc(R2.sequence)[(min_len - offset):]
    q_post = R2.quality[::-1][(min_len - offset):]

    # Take base with best quality
    idx = q2 > q1
    c = r1
    c[idx] = r2[idx]
    c = "".join(c)
    q = q1
    q[idx] = q2[idx]
    q = "".join(q)

    C = R1
    C.sequence = r_pre + c + r_post
    C.quality = q_pre + q + q_post

    return (C)

def paired_consensus_fastq(fastq_in,outfile,min_agreement=.8):
    """
    Takes overlapping paired end reads in single interleaved fastq file
    R1 : read 1
    R2 : read 2
    offset : the expected offset between the start of R1 and end of R2 (0 if completely overlapping)
    """
    with pysam.FastxFile(fastq_in) as infq, \
        gzip.open(outfile,"wt") as outfq:

        nbad = 0
        read2 = -1
        for read in infq:

            # Assuming interleaved fastq,
            if read2 is None:
                read2 = read
            else:
                read1 = read
                read2 = None
                continue

            offset = align_pair(read1.sequence,rc(read2.sequence))

            if offset==-1:
                nbad +=1
            else:
                C = create_consensus_from_readpair(read1,read2,offset=offset)
                outfq.write(str(C) + '\n')

        sys.stderr.write(f'Filtered {nbad} disagreeing reads\n')





def align_pair(s1, s2, min_overlap=50, max_mismatch=2):
    """
    Given two sequences corresponding to paired reads, will check if there is an offset that will align them
    We assume s1 starts to the left of s2 since we have trimmed sequences outside of adapter sites
    Looks for region of overlap at least min_overlap bases with no more than max_mismatch number of disagreements
    Returns offset as the start position of s2 relative to s1, and returns -1 if no such alignment can be found
    """


    offset = -1
    # Best case, they're already aligned and error-free
    if s1 == s2:
        offset = 0
    else:
        # Next best case, they're not aligned, but we can find an alignment without any mismatches
        overlap_len, s = max_overlap(s1, s2)
        if overlap_len > min_overlap:
            offset = len(s1) - overlap_len

        # Otherwise have to do slower n^2 alignment
        # If min_overlap is set larger than the read length, we will still check the zero offset case
        for i in range(0, max(1, len(s1) - min_overlap)):
            min_len = min(len(s1), len(s2))
            if matches_less_than_k(s1[i:min_len], s2[0:(min_len - i)], max_mismatch=max_mismatch):
                offset = i
                break
    return (offset)


def z_algorithm(s):
    z = [0] * len(s)
    left = right = 0
    for i in range(1, len(s)):
        if i <= right:
            z[i] = min(right - i + 1, z[i - left])
        while i + z[i] < len(s) and s[z[i]] == s[i + z[i]]:
            z[i] += 1
        if i + z[i] - 1 > right:
            left, right = i, i + z[i] - 1
    return z


def max_overlap(s1, s2):
    # s1 suffix with s2 prefix
    combined1 = s2 + "#" + s1
    z1 = z_algorithm(combined1)
    overlap = max([z for i, z in enumerate(z1[len(s2) + 1:]) if z == len(z1) - len(s2) - 1 - i], default=0)

    return (overlap, s1 + s2[overlap:])


def matches_less_than_k(s1, s2, max_mismatch):
    mismatches = 0
    for a, b in zip(s1, s2):
        if a != b:
            mismatches += 1
            if mismatches > max_mismatch:
                return (False)
    return (True)


if __name__ == "__main__":
    fastq_in = sys.argv[1]
    fastq_out = sys.argv[2]
    paired_consensus_fastq(fastq_in=fastq_in,outfile=fastq_out)