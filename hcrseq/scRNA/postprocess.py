from hcrseq.common.ref import Reference
from concurrent.futures import ProcessPoolExecutor
import os
import pysam
import parasail
import re
import sys

def realign_reporter_transcripts(bam_in,bam_out,
                           reference,
                                 hr_pathway="HRR",
                                 dsb_pathway="NHEJ"):
    """
    For each read aligned to HR reporter: checks whether the HR donor sequence is contained. If not
    realigns to the DSB reporter
    """

    ref = Reference.load(reference)
    matrix = parasail.Matrix("dnafull")

    hr_reporter = ref.get_pathway(hr_pathway).reporter
    dsb_reporter = ref.get_pathway(dsb_pathway).reporter

    # For read in bam
    ## if read is aligned to the HR pathway, check whether the sequence at hr_reporter.barcode_position
    ## is aligned to and exactly matches hr_reporter.barcode
    ## if not, realign it to dsb_reporter.sequence
    # Print read to output bam

    print("Starting realignment process...", file=sys.stderr)

    with (pysam.AlignmentFile(bam_in, "rb") as input_bam, \
        pysam.AlignmentFile(bam_out, "wb", header=input_bam.header) as output_bam_file):

        for read in input_bam:

            # Determine reference sequence
            ref_seq = ref.get_reporter(read.reference_name).sequence
            ref_name = read.reference_name

            if (read.reference_name in [hr_reporter.name,dsb_reporter.name]) and \
                (not read.is_unmapped) and \
                    (read.mapping_quality == 255) and \
                any(op in read.cigarstring for op in ['S','D','N']):

                if ((read.reference_name==hr_reporter.name) and \
                        not check_barcode_match(read,
                        hr_reporter.sequence,
                        hr_reporter.barcode_position,
                        len(hr_reporter.barcode))):
                    ref_seq = dsb_reporter.sequence
                    ref_name = dsb_reporter.name


                read = realign_read(read,ref_seq,ref_name,matrix=matrix)

                read.reference_id = input_bam.get_tid(read.reference_name)

            if not read.has_tag('MD'):
                try:
                    md_tag, nm_tag = calculate_md_and_nm(read, ref_seq)
                    read.set_tag('MD', md_tag, 'Z')
                    read.set_tag('NM', nm_tag, 'i')
                except Exception as e:
                    # Fallback if reference indexing fails
                    print(f"Warning: Could not calculate MD for {read.query_name}: {e}", file=sys.stderr)

            output_bam_file.write(read)

def check_barcode_match(read,ref_seq,barcode_position,barcode_length):
    aligned_pairs = read.get_aligned_pairs(matches_only=False)

    found = False
    for read_pos,ref_pos in aligned_pairs:
        if ref_pos == barcode_position:
            found = True

        if found:
            # Found an indel
            if (read_pos is None) or (ref_pos is None):
                return False

            # Found a mismatch
            if read.query_sequence[read_pos].lower() != ref_seq[ref_pos].lower():
                return False

            # If we got all the way to the end, return True
            if ref_pos == (barcode_position+barcode_length):
                return True

        # If we finish the loop without returning True, then it wasn't fully covered
    return False


def parse_parasail_cigar_to_pysam(parasail_cigar_string):
    OP_MAP = {'M': 0, '=': 0, 'X': 0, 'I': 1, 'D': 2}
    raw_tuples = []
    cigar_pairs = re.findall(r'(\d+)([=XMID])', parasail_cigar_string)

    for length, operation in cigar_pairs:
        raw_tuples.append((OP_MAP[operation], int(length)))

    if not raw_tuples: return []

    # Merge consecutive identical ops
    merged = []
    curr_op, curr_len = raw_tuples[0]
    for next_op, next_len in raw_tuples[1:]:
        if next_op == curr_op:
            curr_len += next_len
        else:
            merged.append((curr_op, curr_len))
            curr_op, curr_len = next_op, next_len
    merged.append((curr_op, curr_len))
    return merged


def realign_read(read, ref_seq, ref_name, matrix):
    # Perform the alignment
    result = parasail.sw_trace(read.query_sequence, ref_seq, open=10, extend=1, matrix=matrix)

    if result.score <= 0 or result.cigar is None:
        read.is_unmapped = True
        return read

    # 1. Parse and Merge (Turn 1M1M1M into 3M)
    core_cigar = parse_parasail_cigar_to_pysam(result.cigar.decode.decode())

    # 2. Strip leading/trailing Deletions (D=2)
    # Any length removed from the FRONT must be added to reference_start
    actual_ref_start = result.cigar.beg_ref

    while core_cigar and core_cigar[0][0] == 2:
        op_len = core_cigar.pop(0)[1]
        actual_ref_start += op_len

    while core_cigar and core_cigar[-1][0] == 2:
        core_cigar.pop(-1)

    # 3. Handle Soft-Clipping (S=4)
    # The sum of M, I, S, =, X must equal query sequence length
    full_cigar = []

    # Leading Soft Clip: bases in the read before the aligned portion
    if result.cigar.beg_query > 0:
        full_cigar.append((4, result.cigar.beg_query))

    # The cleaned core alignment
    full_cigar.extend(core_cigar)

    # Trailing Soft Clip: remaining bases in the read
    # We count only ops that consume query bases (M, I, S, =, X)
    # pysam codes: M=0, I=1, S=4, =7, X=8. (D=2 and N=3 do NOT consume query)
    query_consumed_in_core = sum(length for op, length in core_cigar if op in [0, 1, 4, 7, 8])
    total_query_accounted_for = result.cigar.beg_query + query_consumed_in_core
    trailing_clip = len(read.query_sequence) - total_query_accounted_for

    if trailing_clip > 0:
        full_cigar.append((4, trailing_clip))

    # Final read updates
    read.cigartuples = full_cigar
    read.reference_start = actual_ref_start
    read.reference_name = ref_name

    # Set alignment score and clean up legacy tags
    read.set_tag('AS', result.score, 'i')
    for tag in ['MD','NM','MC', 'ts']:
        if read.has_tag(tag):
            read.set_tag(tag, None)

    return read


def calculate_md_and_nm(read, ref_seq):
    query_seq = read.query_sequence
    ref_pos = read.reference_start
    query_pos = 0

    md_parts = []
    matches = 0
    nm_count = 0

    for op, length in read.cigartuples:
        if op == 0:  # M (Match or Mismatch)
            for _ in range(length):
                if ref_pos >= len(ref_seq) or query_pos >= len(query_seq):
                    break

                q_base = query_seq[query_pos].upper()
                r_base = ref_seq[ref_pos].upper()

                if q_base == r_base:
                    matches += 1
                else:
                    md_parts.append(str(matches))
                    md_parts.append(r_base)
                    matches = 0
                    nm_count += 1

                query_pos += 1
                ref_pos += 1

        elif op == 1:  # I (Insertion)
            query_pos += length
            nm_count += length

        elif op == 2:  # D (Deletion)
            md_parts.append(str(matches))
            matches = 0
            deleted_bases = ref_seq[ref_pos:ref_pos + length].upper()
            md_parts.append(f"^{deleted_bases}")
            ref_pos += length
            nm_count += length

        elif op == 3:  # N (Skipped region / Intron)
            # IMPORTANT: The MD tag does NOT include the N length.
            # We don't add length to 'matches', we just jump the ref_pos.
            ref_pos += length

        elif op == 4:  # S (Soft clip)
            query_pos += length

    # Final match tally
    md_parts.append(str(matches))

    # Join and clean up
    md_str = "".join(md_parts)

    # Remove leading zeros unless the whole string is 0 (e.g., "0A5" -> "A5")
    # This matches standard SAMtools behavior for consecutive mismatches.
    md_str = re.sub(r'^0+(?=[A-Z\^])', '', md_str)
    # Remove internal 0s between a base and another base/deletion
    md_str = re.sub(r'(?<=[A-Z\^])0+(?=[A-Z\^])', '', md_str)

    return md_str, nm_count

def _subsample_contig(bam,contig,target,seed):
    contig_bam = f"{bam}.{contig}.subsampled.bam"
    total_reads = int(pysam.view("-c",bam,contig))
    fraction = min(target / total_reads, 1.0) if total_reads else 0
    pysam.view("-bh",
               "-o", contig_bam,
               "-s", str(seed + fraction),
               bam,
               contig,
               catch_stdout=False)
    return contig_bam

def create_subsampled_bam(bam,reference,target,seed=0):
    """
    Creates a subsampled bam with up to [target] randomly selected
    reads from each reporter contig in [reference]
    """
    with pysam.AlignmentFile(bam) as f:
        sort_order = f.header.to_dict().get("HD", {}).get("SO")
    if sort_order != "coordinate":
        raise ValueError(f"{bam} must be coordinate-sorted (found SO:{sort_order}) "
                          "for the per-contig merge below to come out sorted")

    ref = Reference.load(reference)
    contigs = [reporter.name for reporter in ref.reporters]

    with ProcessPoolExecutor(max_workers=len(contigs)) as pool:
        contig_bams = list(pool.map(_subsample_contig,
                                     [bam] * len(contigs),
                                     contigs,
                                     [target] * len(contigs),
                                     [seed] * len(contigs)))

    # Each contig_bam is a region-filtered subset of the coordinate-sorted input, so it
    # is itself already coordinate-sorted; merging pre-sorted, single-contig inputs with
    # samtools' default (position) merge mode yields a fully coordinate-sorted bam directly,
    # with no separate sort pass needed.
    pysam.merge("-f", "-o", f"{bam}.subsampled.bam", *contig_bams)
    pysam.index(f"{bam}.subsampled.bam")

    for contig_bam in contig_bams:
        os.remove(contig_bam)


if __name__ == "__main__":
    bam_in = sys.argv[1]
    bam_out = sys.argv[2]
    ref = sys.argv[3]
    realign_reporter_transcripts(bam_in,bam_out,ref)


