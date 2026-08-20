from hcrseq.common.quantify import HCRseqQuantifier


def count_umis(bam, reference,min_mapq = 5,min_mh=3):
    """
    Counts UMIs (repaired and otherwise) from bam file, using the reporters
    defined in an HCR-seq Reference object
    """

    q = HCRseqQuantifier(reference,tags=None,min_mh=min_mh)
    q.count_umis(bam,min_mapq=min_mapq)

    counts = q.get_counts()
    mismatch_dist = {name: counter.mismatch_dist for name, counter in q.counters.items()}
    del_df = q.get_indel_df('deletions')
    ins_df = q.get_indel_df('insertions')

    return(counts,del_df,ins_df,mismatch_dist)


def quantify_repair(counts,reference):
    """
    Computes repair pathway measurements for each pathway defined in the reference
    """

    repair_measurements = dict()
    for pathway in reference.pathways:
        value, sd = pathway.calculate_repair(counts)
        repair_measurements[pathway.name] = {'value': value, 'sd': sd}

    return(repair_measurements)
