version 1.0

workflow hcrseq_scRNA_quantify {
  input {
    File bam
    File bai
    File h5_file
    File ref_path
    String outstem
  }

  call postprocess_reporter_bam {
    input:
      bam = bam,
      bai = bai,
      outstem = outstem,
      ref_path = ref_path
  }

  call quantify {
    input:
      h5_file = h5_file,
      bam = postprocess_reporter_bam.postprocessed_bam,
      bai = postprocess_reporter_bam.postprocessed_bai,
      outstem = outstem,
      ref_path = ref_path
  }

  output {
    File postprocessed_bam = postprocess_reporter_bam.postprocessed_bam
    File postprocessed_bai = postprocess_reporter_bam.postprocessed_bai
    File subsampled_bam = postprocess_reporter_bam.subsampled_bam
    File subsampled_bai = postprocess_reporter_bam.subsampled_bai
    File repair_annotated_h5ad = quantify.repair_annotated_h5ad
  }
}

task postprocess_reporter_bam {
  input {
    File bam
    File bai
    File ref_path
    String outstem
    Int threads = 8
    Int mem_per_thread_gb = 4
    Int subsample_target = 10000
    Int disk_GB = 200
  }

  command <<<
    set -euo pipefail

    hcrseq scrna postprocess-reporter-bam \
      --bam ~{bam} \
      --ref_path ~{ref_path} \
      --outstem ~{outstem} \
      --threads ~{threads} \
      --mem-per-thread ~{mem_per_thread_gb}G \
      --subsample_target ~{subsample_target}
  >>>

  output {
    File postprocessed_bam = "~{outstem}_postprocessed.sorted.bam"
    File postprocessed_bai = "~{outstem}_postprocessed.sorted.bam.bai"
    File subsampled_bam = "~{outstem}_postprocessed.sorted.bam.subsampled.bam"
    File subsampled_bai = "~{outstem}_postprocessed.sorted.bam.subsampled.bam.bai"
  }

  runtime {
    docker: "gcr.io/broad-getzlab-fmhcrsparc/hcrseq:v1.0"
    memory: ceil(1.2 * mem_per_thread_gb * threads) + " GB"
    disks: "local-disk " + disk_GB + " HDD"
    cpu: threads
    preemptible: 2
  }
}

task quantify {
  input {
    File h5_file
    File bam
    File bai
    File ref_path
    String outstem
    Int min_mapq = 5
    Int disk_GB = 200
    Int mem_GB = 32
  }

  command <<<
    set -euo pipefail

    hcrseq scrna quantify --h5_file ~{h5_file} \
      --bam ~{bam} \
      --ref_path ~{ref_path} \
      --outstem ~{outstem} \
      --min_mapq ~{min_mapq}
  >>>

  output {
    File repair_annotated_h5ad = "~{outstem}.h5ad"
  }

  runtime {
    docker: "gcr.io/broad-getzlab-fmhcrsparc/hcrseq:v1.0"
    memory: mem_GB + " GB"
    disks: "local-disk " + disk_GB + " HDD"
    cpu: 4
    preemptible: 2
  }
}
