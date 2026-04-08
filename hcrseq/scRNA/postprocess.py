from hcrseq.common.ref import Reference
from hcrseq.common.util import check_perfect_match
import pysam
import parasail
import re
import sys

def realign_reporter_transcripts(bam_in,bam_out,
                           reference):
    """
    For each read aligned to HR reporter: checks whether the HR donor sequence is contained. If not
    realigns to the DSB reporter
    """

    ref = Reference.load(reference)
    matrix = parasail.Matrix("dnafull")

    hr_reporter = ref.get_pathway('HR').reporter
    dsb_reporter = ref.get_pathway('NHEJ').reporter

    # For read in bam
    ## if read is aligned to the HR pathway, check whether the sequence at hr_reporter.barcode_position
    ## is aligned to and exactly matches hr_reporter.barcode
    ## if not, realign it to dsb_reporter.sequence
    # Print read to output bam

    print("Starting realignment process...", file=sys.stderr)

    with pysam.AlignmentFile(bam_in, "rb") as input_bam, \
        pysam.AlignmentFile(bam_out, "wb", header=input_bam.header) as output_bam_file:

        for read in input_bam:

            if (not read.is_unmapped) and \
                any(op in read.cigarstring for op in ['S','D','N']) and \
                not check_perfect_match(read,
                    hr_reporter.barcode_position,
                    hr_reporter.barcode_position+len(hr_reporter.barcode)
                                        ):

                if read.reference_id==hr_reporter.name:
                    ref_seq = dsb_reporter.sequence
                else:
                    ref_seq = ref.get_reporter(read.reference_id)

                read = realign_read(read,ref_seq,matrix=matrix)

        output_bam_file.write(read)

def parse_parasail_cigar_to_pysam(parasail_cigar_string):
    """
    Parses a CIGAR string from parasail's trace output into the tuple format
    required by pysam.

    Parasail's trace CIGAR uses 'M', '=', or 'X' for matches/mismatches,
    'I' for insertions, and 'D' for deletions. The SAM/BAM format uses a
    single 'M' operation for both matches and mismatches.

    Args:
        parasail_cigar_string (str): The CIGAR string from parasail (e.g., '15M1I34M').

    Returns:
        list: A list of tuples, where each tuple is (operation, length),
              compatible with pysam's `read.cigartuples`.
    """
    # pysam operation codes: M=0, I=1, D=2, S=4, etc.
    # We will map parasail's match/mismatch ops ('M', '=', 'X') to pysam's 'M' (0).
    OP_MAP = {'M': 0, '=': 0, 'X': 0, 'I': 1, 'D': 2}

    pysam_tuples = []
    # Regex to find all pairs of (length, operation) in the string
    cigar_pairs = re.findall(r'(\d+)([=XMID])', parasail_cigar_string)

    for length, operation in cigar_pairs:
        if operation in OP_MAP:
            pysam_tuples.append((OP_MAP[operation], int(length)))

    return pysam_tuples

def realign_read(read,ref_seq,matrix):
    # Perform a global (Needleman-Wunsch) alignment of the entire read sequence
    # against the entire reference contig sequence.
    result = parasail.sw_stats_strip_32(read.query_sequence,
                               ref_seq,
                               open=10,
                               extend=1,
                               matrix=matrix)

    # Convert parasail's CIGAR string to pysam's tuple format
    new_cigar_tuples = parse_parasail_cigar_to_pysam(result.cigar.decode.decode())

    # Update the read's attributes with the new alignment info
    read.cigartuples = new_cigar_tuples
    # For a global alignment against the entire contig, the start position is 0
    read.reference_start = 0
    # Set the alignment score (AS) tag
    read.set_tag('AS', result.score, 'i')
    # The MD tag is now invalid, so we remove it.
    if read.has_tag('MD'):
        read.set_tag('MD', None)

    return(read)


if __name__ == "__main__":
    bam_in = sys.argv[1]
    bam_out = sys.argv[2]
    ref = sys.argv[3]
    realign_reporter_transcripts(bam_in,bam_out,ref)


