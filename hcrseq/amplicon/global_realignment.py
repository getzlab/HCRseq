import pysam
import parasail
import re
import sys


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


def realign_bam(input_bam_path, output_bam_path, ref_fasta_path):
    """
    Realigns soft-clipped reads from a BAM file using a global aligner and
    writes the results to a new BAM file.

    Args:
        input_bam_path (str): Path to the input BAM file.
        output_bam_path (str): Path for the output BAM file.
        ref_fasta_path (str): Path to the reference FASTA file.
    """
    # Define the scoring matrix for the alignment
    # dnafull provides scores for matches and mismatches between any DNA bases.
    matrix = parasail.Matrix("dnafull")

    # Cache to store reference sequences in memory to avoid redundant file reads
    ref_sequences = {}

    print("Starting realignment process...", file=sys.stderr)

    with pysam.FastaFile(ref_fasta_path) as ref_file, \
            pysam.AlignmentFile(input_bam_path, "rb") as input_bam, \
            pysam.AlignmentFile(output_bam_path, "wb", header=input_bam.header) as output_bam_file:

        for read in input_bam:
            # Check if the read is mapped and has soft clipping in its CIGAR string
            if not read.is_unmapped and 'S' in read.cigarstring:
                try:
                    ref_name = read.reference_name

                    # Fetch the reference sequence if not already in our cache
                    if ref_name not in ref_sequences:
                        ref_sequences[ref_name] = ref_file.fetch(ref_name)

                    ref_seq = ref_sequences[ref_name]

                    # Perform a global (Needleman-Wunsch) alignment of the entire read sequence
                    # against the entire reference contig sequence.
                    result = parasail.nw_trace(read.query_sequence,
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

                    # Write the modified read to the output file
                    output_bam_file.write(read)

                except Exception as e:
                    print(f"Could not process read {read.query_name}: {e}", file=sys.stderr)
                    # Write the original read if an error occurs
                    output_bam_file.write(read)
            else:
                # If the read has no soft clipping, write it to the output file unchanged
                output_bam_file.write(read)
if __name__ == "__main__":
    bam_in = sys.argv[1]
    bam_out = sys.argv[2]
    ref = sys.argv[3]
    realign_bam(bam_in,bam_out,ref)