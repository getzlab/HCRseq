import numpy as np
import pandas as pd
import pysam
from collections import Counter, defaultdict

from hcrseq.common.util import check_perfect_match


class UMICounter(object):
    def __init__(self,reporter,tags,min_mh=3):
        """
        :param reporter: Reporter being counted
        :param tags: Read tags being counted, (CB,UB) for scRNA, None for amplicon
        """
        self.reporter = reporter

        ref_len = len(self.reporter.sequence)

        if tags is None:
            self.counts = Counter()
        else:
            # scRNA mode, need two layer dictionary for CB/UMI
            self.counts = defaultdict(lambda: defaultdict(Counter))

        self.mismatch_dist = {'nsnp': np.zeros(ref_len),
                        'ndel': np.zeros(ref_len),
                        'nins': np.zeros(ref_len),
                        'N': np.zeros(ref_len)}  # For each position, counts number of mismatches
        self.deletions = list()
        self.insertions = list()
        self.tags = tags
        self.min_mh = min_mh

    def count(self,read):

        if self.tags is not None:
            if not all(read.has_tag(tag) for tag in self.tags):
                return

        tag_vals = [read.get_tag(tag) for tag in self.tags] if self.tags is not None else None

        self.inc('total',tag_vals)

        self.count_mismatches(read)

        if self.reporter.lesion_position is not None:
            # Report any deletion spanning lesion
            is_deleted, deletion_info = check_deletion(read,
                                                       self.reporter.lesion_position,
                                                       self.reporter.sequence,
                                                       allow_after_base=self.reporter.unrepaired_base == 'DSB')
            if is_deleted:
                self.inc('del',tag_vals)

                if deletion_info['microhomology_length'] >= self.min_mh:
                    self.inc('del_mh',tag_vals)
                else:
                    self.inc('del_nomh',tag_vals)

                if self.tags is not None:
                    for tag in self.tags:
                        deletion_info[tag] = read.get_tag(tag)
                self.deletions.append(deletion_info)

            ## Check for insertions
            has_insertion, insertion_info = check_insertion(read, self.reporter.lesion_position)
            if has_insertion:
                self.inc('ins',tag_vals)

                if self.tags is not None:
                    for tag in self.tags:
                        insertion_info[tag] = read.get_tag(tag)

                self.insertions.append(insertion_info)

            ## Check for point mutations
            if isinstance(self.reporter.unrepaired_base, str) and (self.reporter.unrepaired_base in 'ACGT'):
                base_aligned, base = check_base(read, self.reporter.lesion_position)

                if base_aligned:
                    self.inc('repaired',tag_vals,base == self.reporter.repaired_base)
                    self.inc('unrepaired',tag_vals, base == self.reporter.unrepaired_base)
                    self.inc(base,tag_vals)


    def inc(self,key,tag_vals,val=1):
        if self.tags is None:
            # Amplicon mode - just one set of values
            self.counts[key] += val
        else:
            # scRNA mode, need to keep try of CB/UMI
            self.counts[tag_vals[0]][tag_vals[1]][key] += val

    def count_mismatches(self,read):

        aligned_pairs = read.get_aligned_pairs(with_seq=True)

        for i, (query_pos, ref_pos, ref_base) in enumerate(aligned_pairs):

            if ref_pos is None or ref_base is None:
                continue

            if ref_base.islower():
                self.mismatch_dist['nsnp'][ref_pos] += 1
            if query_pos is None:
                self.mismatch_dist['ndel'][ref_pos] += 1
            if (i + 1 < len(aligned_pairs)) and (aligned_pairs[i + 1][1] is None):
                self.mismatch_dist['nins'][ref_pos] += 1
            self.mismatch_dist['N'][ref_pos] += 1

class HCRseqQuantifier(object):
    def __init__(self,reference,tags,min_mh=3):
        self.reference = reference
        self.tags = tags
        self.counters = dict()

        for reporter in self.reference.reporters:
            self.counters[reporter.name] = UMICounter(reporter,tags,min_mh=min_mh)

    def count_umis(self,bam,min_mapq = 5,require_exact_bc=True):

        with pysam.AlignmentFile(bam) as bam_in:
            for reporter_counter in self.counters.values():
                reporter = reporter_counter.reporter
                for read in bam_in.fetch(reporter.name):

                    if read.mapq < min_mapq:
                        continue

                    # Check barcode sequence is correct
                    if require_exact_bc and (reporter.barcode is not None):
                        if not check_perfect_match(read,
                                                   reporter.barcode_position,
                                               reporter.barcode_position+reporter.barcode_len):
                            continue

                    reporter_counter.count(read)
    def get_indel_df(self,indel_type):
        df = pd.concat([pd.DataFrame(getattr(counter,indel_type)) for counter in self.counters.values()],
                       axis=0)
        return(df)

    def get_counts(self):
        counts = {counter.reporter.name : counter.counts
                  for counter in self.counters.values()}
        return counts

    def get_cell_counts(self):
        """
        Collapses UMI-level counts into a single Counter per (reporter, cell barcode).
        Requires the quantifier to have been constructed with (cell_tag, umi_tag) tags.
        """
        if self.tags is None:
            raise ValueError('get_cell_counts requires per-cell tags, e.g. ("CB", "UB")')

        return {reporter_name: {cb: self._collapse_counts(umi_counts)
                                for cb, umi_counts in reporter_counts.items()}
                for reporter_name, reporter_counts in self.get_counts().items()}

    def quantify_repair_per_cell(self):
        """
        Computes repair pathway measurements independently for each cell barcode.
        Returns dict[pathway_name][cell_barcode] -> {'value': ..., 'sd': ...}
        """
        cell_counts = self.get_cell_counts()
        cell_barcodes = {cb for reporter_counts in cell_counts.values() for cb in reporter_counts}

        repair_measurements = {pathway.name: {} for pathway in self.reference.pathways}
        for cb in cell_barcodes:
            counts = {reporter_name: reporter_counts.get(cb, Counter())
                     for reporter_name, reporter_counts in cell_counts.items()}

            for pathway in self.reference.pathways:
                counts.setdefault(pathway.reporter.name, Counter())
                counts.setdefault(pathway.control.name, Counter())
                value, sd = pathway.calculate_repair(counts)
                repair_measurements[pathway.name][cb] = {'value': value, 'sd': sd}

        return repair_measurements

    def quantify_repair(self):

        counts = self.get_counts()
        aggregated_counts = {}

        for reporter_name, reporter_counts in counts.items():
            aggregated_counts[reporter_name] = self._collapse_counts(reporter_counts)

        repair_measurements = {}
        for pathway in self.reference.pathways:
            if pathway.reporter.name not in aggregated_counts:
                aggregated_counts[pathway.reporter.name] = Counter()
            if pathway.control.name not in aggregated_counts:
                aggregated_counts[pathway.control.name] = Counter()

            value, sd = pathway.calculate_repair(aggregated_counts)
            repair_measurements[pathway.name] = {'value': value, 'sd': sd}

        self.repair_measurements = repair_measurements
        return repair_measurements

    @staticmethod
    def _collapse_counts(counts):
        if isinstance(counts, Counter):
            return counts.copy()

        if not isinstance(counts, dict):
            return Counter()

        if counts and isinstance(next(iter(counts.values())), dict):
            collapsed = Counter()
            for value in counts.values():
                collapsed.update(HCRseqQuantifier._collapse_counts(value))
            return collapsed

        return Counter(counts)


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
        st = int(x[np.where(aligned_idx & (x[:, 1] < pos))[0][-1], 1])
        en = int(x[np.where(aligned_idx & (x[:, 1] > pos))[0][0], 1])

        deleted_sequence = ref[(st + 1):en]
        flank_sequence = ref[en:(en + len(deleted_sequence))]
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
