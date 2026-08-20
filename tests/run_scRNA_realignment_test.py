import time
import pysam
from hcrseq.scRNA.postprocess import realign_reporter_transcripts

if __name__ == '__main__':
    bam_in = "scRNA_test_data/experiment11_reporters_1pct.bam"
    bam_out = "scRNA_test_results/experiment11_postprocessed.bam"
    ref = "../reference/HCRseq_v0.4.1/HCRseq_v0.4.1.pkl"

    print(f"--- Starting realignment: {bam_in} ---")

    start_time = time.perf_counter()

    realign_reporter_transcripts(bam_in, bam_out, ref)

    end_time = time.perf_counter()
    duration = end_time - start_time

    try:
        with pysam.AlignmentFile(bam_in, "rb") as bam:
            # Note: .mapped might be more accurate if you are only timing mapped reads
            total_reads = bam.mapped + bam.unmapped
    except Exception:
        total_reads = None

    print("-" * 30)
    print(f"Realignment completed in: {duration:.2f} seconds")

    if total_reads:
        reads_per_sec = total_reads / duration
        print(f"Processed ~{total_reads} reads at {reads_per_sec:.2f} reads/sec")
    print(f"Output saved to: {bam_out}")
    print("-" * 30)