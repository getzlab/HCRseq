

BP = {'A' : 'T',
      'C' : 'G',
      'G' : 'C',
      'T' : 'A',
      'N' : 'N'}

def rc(seq):
    return("".join([BP[b] for b in reversed(seq)]))


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
