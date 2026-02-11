import numpy as np
import pandas as pd
import pysam
from collections import Counter, defaultdict


def count_umis(bam,
               ref_fasta,
               lesion_info,
               min_mapq = 5,
               min_mh=3,
               require_exact_bc=True,
               count_mismatch_distribution=True,
               tags = None):
    """
    Counts UMIs (repaired and otherwise) from bam file
    """

    L = pd.read_csv(lesion_info, sep='\t', index_col='reporter')
    ref = pysam.FastaFile(ref_fasta)

    counts = defaultdict(Counter)
    mismatch_dist = dict() # For each position, counts number of mismatches
    deletions = list()
    insertions=list()

    # First pass to calculate total reads
    with pysam.AlignmentFile(bam) as bam_in:
        for contig in L['contig'].unique():

            if count_mismatch_distribution:
                ref_len = len(ref.fetch(contig))
                mismatch_dist[contig] = {'nsnp': np.zeros(ref_len),
                                     'ndel': np.zeros(ref_len),
                                     'nins': np.zeros(ref_len),
                                     'N': np.zeros(ref_len)}

            for read in bam_in.fetch(contig):

                if read.mapq < min_mapq:
                    continue

                # Check barcode sequence is correct
                if require_exact_bc and (not check_reporter_barcode(read,ref.fetch(contig))):
                    continue

                counts[contig]['total'] += 1

                if count_mismatch_distribution:
                    count_mismatches(read, mismatch_dist)

                if contig in L.index:

                    lesion = L.loc[contig]
                    pos = lesion['position'] - 1

                    ## Check for deletions
                    is_deleted, deletion_info = check_deletion(read,
                                                               pos,
                                                               ref,
                                                               allow_after_base = lesion['unrepaired_base']=='DSB')
                    if is_deleted:
                        counts[contig]['del'] += 1

                        if deletion_info['microhomology_length'] >= min_mh:
                            counts[contig]['del_mh'] += 1
                        else:
                            counts[contig]['del_nomh'] += 1
                        deletions.append(deletion_info)

                    ## Check for insertions
                    has_insertion,insertion_info = check_insertion(read, pos)
                    if has_insertion:
                        counts[contig]['ins'] += 1

                        insertions.append(insertion_info)

                    ## If lesion is a point mutation, then check it
                    if lesion['unrepaired_base'] in 'ACGT':
                        base_aligned,base = check_base(read, pos)

                        if base_aligned:
                            counts[contig]['repaired'] += base == lesion['repaired_base']
                            counts[contig]['unrepaired'] += base == lesion['unrepaired_base']
                            counts[contig][base] += 1

    del_df = pd.DataFrame(deletions)
    ins_df = pd.DataFrame(insertions)

    return(counts,del_df,ins_df,mismatch_dist)

def count_mismatches(read,mismatch_dist):
    contig = read.reference_name

    aligned_pairs = read.get_aligned_pairs(with_seq=True)

    for i,(query_pos, ref_pos, ref_base) in enumerate(aligned_pairs):

        if ref_pos is None or ref_base is None:
            continue

        if ref_base.islower():
            mismatch_dist[contig]['nsnp'][ref_pos] += 1
        if query_pos is None:
            mismatch_dist[contig]['ndel'][ref_pos] += 1
        if (i+1 < len(aligned_pairs)) and (aligned_pairs[i+1][1] is None):
            mismatch_dist[contig]['nins'][ref_pos] += 1
        mismatch_dist[contig]['N'][ref_pos] += 1

def check_reporter_barcode(read,ref_seq,barcode_len=6,barcode_offset = 23):
    """
    TODO : this only works for amplicon version, need to make work for scRNA
    """
    barcode_st = len(ref_seq) - barcode_offset - barcode_len
    return(check_perfect_match(read,barcode_st,barcode_st + barcode_len-1))

def check_perfect_match(read,ref_st,ref_en):
    """
    Returns true if a read perfectly aligns to a particular reference region with no mismatches/indels
    """

    aligned_pairs = read.get_aligned_pairs(with_seq=True,matches_only=False)

    found = False
    for read_pos,ref_pos,base in aligned_pairs:

        if ref_pos ==ref_st:
            found=True

        if found:
            # Found an indel
            if (read_pos is None) or (ref_pos is None):
                return False

            # Found a mismatch
            if base.islower():
                return False

            # If we got all the way to the end, return True
            if ref_pos == ref_en:
                    return True

    # If we finish the loop without returning True, then it wasn't fully covered
    return False


def check_deletion(read,pos,ref,allow_after_base=False):

    contig = read.reference_name
    aligned_pairs = read.get_aligned_pairs(matches_only=False)

    x = np.array(aligned_pairs).astype(float)
    x[x == None] = np.nan

    is_spanning_del = any((np.isnan(x[:, 0])) & (x[:, 1] == pos) )

    # the position represents the base before the breakpoint - we want to include deletions including the base on either side
    if allow_after_base & (not is_spanning_del) and any((np.isnan(x[:, 0])) & (x[:, 1] == (pos+1))):
        is_spanning_del =True
        pos +=1

    deletion_info = {}
    if is_spanning_del:
        aligned_idx = ~np.isnan(x[:, 0])
        st = x[np.where(aligned_idx & (x[:, 1] < pos))[0][-1], 1]
        en = x[np.where(aligned_idx & (x[:, 1] > pos))[0][0], 1]

        deleted_sequence = ref.fetch(contig, st + 1, en)
        flank_sequence = ref.fetch(contig, en, en + len(deleted_sequence))
        max_mh_len = min(len(deleted_sequence),len(flank_sequence))
        mh_len = sum(np.cumprod(np.array(list(deleted_sequence[0:max_mh_len])) == np.array(list(flank_sequence[0:max_mh_len]))))

        # Record entry
        deletion_info = {
            'UMI': read.get_tag('UB'),
            'contig': read.reference_name,
            'start': st+1,
            'end': en+1,
            'deletion_length': len(deleted_sequence),
            'deleted_sequence': deleted_sequence,
            'flank_sequence': flank_sequence,
            'microhomology_length': mh_len}


    return (is_spanning_del, deletion_info)


def check_insertion(read, pos):
    aligned_pairs = read.get_aligned_pairs(matches_only=False)

    for i, (read_idx, ref_idx) in enumerate(aligned_pairs):
        # We are looking for an insertion that starts right after our target position
        has_insertion = ref_idx == pos and (i + 1 < len(aligned_pairs) and aligned_pairs[i + 1][1] is None)

        if has_insertion:
            # We found an insertion. Now find its full length.
            st = aligned_pairs[i + 1][0]
            en = st+1
            for j in range(i + 2, len(aligned_pairs)):
                if aligned_pairs[j][1] is None:
                    en = aligned_pairs[j][0]+1
                else:
                    insertion_sequence = read.query_sequence[st:en]
                    insertion_info = {'UMI': read.get_tag('UB'),
                                      'contig' : read.reference_name,
                                      'pos' : pos+1,
                                      'read_st' : st+1,
                                      'read_en' : en+1,
                                      'insertion_length' : len(insertion_sequence),
                                      'insertion_sequence' : insertion_sequence}
                    return (has_insertion, insertion_info)
    return (False, {})


def check_base(read,pos):
    aligned_pairs = read.get_aligned_pairs(matches_only=True)

    for i, (read_idx, ref_idx) in enumerate(aligned_pairs):
        if ref_idx==pos:
            return(True,read.query_sequence[read_idx])
    return(False,"")
