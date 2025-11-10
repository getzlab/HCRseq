import pysam
from collections import Counter
import sys


def deduplicate_v2(tagged_bam,output_bam):

    with (pysam.AlignmentFile(tagged_bam,check_sq=False) as bam_in,
          pysam.AlignmentFile(output_bam,"wb",template=bam_in) as bam_out):

        for read_list in UMIIterator(bam_in):

            counter = Counter([read.query_sequence for read in read_list])
            seq, seq_ct = counter.most_common(1)[0]

            for read in read_list:
                if read.query_sequence==seq:
                    representative_read = read
                    break

            # Record number of reads supporting this UMI
            representative_read.set_tag("NR",seq_ct)

            # Write the representative read to the output bam
            bam_out.write(representative_read)

class UMIIterator:
    """
    Iterates over a bam file with reads sorted by UMI sequence.
    At each iteration, returns a list of reads all belonging to a single UMI
    """

    def __init__(self,bam):
        self.read_iter = bam.fetch(until_eof=True)
        self.reads = list()
        self.UMI = ""
        self.finished = False

        # Call once to load in first read
        self.__next__()

    def __iter__(self):
        return self

    def __next__(self):
        """
        Before call, iterator should have single read stored in self.reads
        When __next__() is called, it keeps loading reads until it (1) finds one with a different UMI, or
                                                                    (2) comes to the end
        """

        if self.finished:
            raise StopIteration

        self.finished = True
        for read in self.read_iter:
            umi = read.get_tag('UR')
            if umi==self.UMI:
                self.reads.append(read)
            else:
                self.finished=False
                break

        reads = self.reads

        # If we found a read with non-matching UMI, load that in for the next call
        # Otherwise StopIteration will be raised on the next call
        if not self.finished:
            self.reads = [read]
            self.UMI = umi

        return(reads)

def usage():
    return("USAGE: dedup.py [read_bam] [bam_out]")

def main():
    if len(sys.argv)<3:
        print(usage())
        print('Got:')
        print(sys.argv)
        return

    tagged_bam, bam_out = sys.argv[1:]

    deduplicate_v2(tagged_bam,bam_out)

if __name__ == "__main__":
    main()