import pandas as pd
import numpy as np
import pysam
import pickle as pkl
from hcrseq.common.util import rc
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class Reporter(object):

    def __init__(self,name,
                 plasmid,
                 sequence,
                 tx_start,
                 tx_end,
                 barcode_position,
                 barcode,
                 lesion_position,
                 unrepaired_base,
                 repaired_base):
        self.name = name
        self.plasmid = plasmid
        self.sequence = sequence
        self.tx_start = tx_start
        self.tx_end = tx_end
        self.barcode_position = barcode_position
        self.barcode = barcode
        self.barcode_len = 6
        self.lesion_position = lesion_position
        self.unrepaired_base = unrepaired_base
        self.repaired_base = repaired_base

    def to_seqrecord(self):
        record = SeqRecord(
            Seq(self.sequence),
            id=self.name,
        )
        return(record)

    def to_amplicon(self,primers):

        # Find the inner part of the amplicon

        # Trim first primer because DSBR has slightly shorter match
        # TODO make this automatic
        trimmed_forward_primer = primers[0][1:]
        amplicon_start = self.sequence.find(trimmed_forward_primer) + len(trimmed_forward_primer)
        amplicon_end = self.sequence.find(rc(primers[1]))

        # Convert to amplicon sequence
        # Note the full primer sequence will appear regardless of if a partial match in the plasmid sequence
        self.sequence = primers[0] + self.sequence[amplicon_start:amplicon_end] + rc(primers[1])
        self.lesion_position = self.lesion_position - amplicon_start + len(primers[0])
        self.barcode_position = self.barcode_position - amplicon_start + len(primers[0])


class Reference:

    def __init__(self,reporters,pathways,primers = None):
        self.reporters = reporters
        self.pathways = pathways

        if primers is not None:
            self.primers = primers

    def write_amplicon_fasta(self,out_fasta):
        with open(out_fasta, 'w') as fout:
            for _, reporter in self.reporters:
                # Write it to fasta
                r = pysam.FastxRecord()
                r.name = reporter.name
                r.sequence = reporter.sequence

                fout.write(str(r) + '\n')

    def write_gtf(self):
        pass

    def write_pickle(self,out_pickle):
        with open(out_pickle, 'wb') as handle:
            pkl.dump(self, handle)

    @classmethod
    def build(cls,plasmid_info,reporter_info, pathway_info):
        R = pd.read_csv(reporter_info,sep='\t')
        V = pd.read_csv(plasmid_info,sep='\t')
        P = pd.read_csv(pathway_info,sep='\t')

        R = R.join(V.set_index('name'),on='plasmid',how='left')

        # Update barcode sequences
        idx = ~R['barcode_position'].isna()
        for ind,row in R[idx].iterrows():
            st = int(row['barcode_position'])-1 # 0-indexed barcode position
            en = int(row['barcode_position'] + len(row['barcode']))-1 # 0-indexed position after barcode
            seq = row['sequence']
            R.loc[ind,'sequence'] = seq[0:st] + row['barcode'] + seq[en:]


        reporters = {row['name']: Reporter(**row)
                     for _, row in R.iterrows()}

        # Find the control reporter
        ctrl_idx = P['pathway'] == 'control'
        if sum(ctrl_idx) != 1:
            raise ValueError('There must be exactly one control')
        control_name = P.loc[ctrl_idx, 'reporter'].iloc[0]

        pathways = {row['pathway']: Pathway(name=row['pathway'],
                                            metric=row['metric'],
                                            reporter=reporters[row['reporter']],
                                            control=reporters[control_name]) for _, row in P[~ctrl_idx].iterrows()}

        return(cls(reporters=list(reporters.values()),
                   pathways=list(pathways.values())))

    @classmethod
    def load(cls,in_pickle):
        with open(in_pickle, 'rb') as handle:
            return pkl.load(handle)



class Pathway:

    def __init__(self,name,metric,reporter,control):
        self.name = name
        self.metric = metric
        self.reporter = reporter
        self.control = control

        try:
            self.calculate_repair = getattr(self,'_calculate_' + metric)
        except AttributeError:
            print(f'Error! {metric} has no defined function')

    def _calculate_fraction_repaired(self,
                                     counts,
                                     min_count=10):

        repaired = counts[self.reporter.name]['repaired']
        unrepaired = counts[self.reporter.name]['unrepaired']

        n = repaired + unrepaired

        # Return nans if needed plasmid is not present
        if n < min_count:
            return (np.nan, np.nan)

        p = repaired / n

        sd = float(np.sqrt(p * (1 - p) / (n - 1)))

        return (p, sd)

    def _calculate_abundance(self,counts,key='total'):

        reporter = counts[self.reporter.name][key]
        control = counts[self.control.name]['total']

        if control == 0:
            return (np.nan, np.nan)

        r, sd = self._calculate_ratio_and_sd(reporter, control)
        return (r, sd)

    def _calculate_abundance_mmej(self,counts):
        return (self._calculate_abundance(counts,'del_mh'))

    def _calculate_abundance_nommej(self,counts):

        reporter = counts[self.reporter.name]['total'] - counts[self.reporter.name]['del_mh']
        control = counts[self.control.name]['total']

        r, sd = self._calculate_ratio_and_sd(reporter, control)
        return (r, sd)

    @staticmethod
    def _calculate_ratio_and_sd(reporter, control):
        r = reporter / control
        # Assuming poisson variance and using Taylor expansion
        sd = float(np.sqrt(1 / control ** 2 * (reporter + r ** 2 * control)))
        return (r, sd)


if __name__ == '__main__':
    plasmid_info = '../../reference/HCRseq_v0.4/plasmid_info.txt'
    reporter_info = '../../reference/HCRseq_v0.4/reporter_info.txt'
    pathway_info = '../../reference/HCRseq_v0.4/pathway_calculations.txt'
    r = Reference.build(plasmid_info,reporter_info,pathway_info)
    r.write_pickle('../../reference/HCRseq_v0.4/HCRseq_v0.4.pkl')



