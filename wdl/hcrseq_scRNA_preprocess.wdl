workflow hcrseq_scRNA_quantify{
  call postprocess_reporter_bam
  call quantify {
    input:
      bam = postprocess_reporter_bam.postprocessed_bam,
      bai = postprocess_reporter_bam.postprocessed_bai
  }
}

task postprocess_reporter_bam{
  File bam
  File bai
  File ref_path
  String outstem
  Int? threads = 8
  String? mem_per_thread = "4G"
  Int? subsample_target = 10000
  Int? disk_GB = 200

  command {
    set -euo pipefail

    hcrseq scrna postprocess-reporter-bam \
			--bam ${bam} \
			--ref_path ${ref_path} \
			--outstem ${outstem} \
			--threads ${threads} \
			--mem-per-thread ${mem_per_thread} \
			--subsample_target ${subsample_target}

  }

  output {
    File postprocessed_bam="${outstem}_postprocessed.sorted.bam"
    File postprocessed_bai="${outstem}_postprocessed.sorted.bam.bai"
    File subsampled_bam="${outstem}_postprocessed.sorted.bam.subsampled.bam"
    File subsampled_bai="${outstem}_postprocessed.sorted.bam.subsampled.bam.bai"
  }

  runtime {
    docker: "gcr.io/broad-getzlab-fmhcrsparc/hcrseq:v1.0"
    memory: "16GB"
    disks: "local-disk " + disk_GB + " HDD"
    cpu: threads
    preemptible: "2"
  }


}

task quantify{
  File h5_file
  File bam
  File bai
  File ref_path
  String outstem
  Int? min_mapq = 5
  Int? disk_GB = 200

  command {
    set -euo pipefail

    hcrseq scrna quantify --h5_file ${h5_file} \
			--bam ${bam} \
			--ref_path ${ref_path} \
			--outstem ${outstem} \
			--min_mapq ${min_mapq}

  }


  output {
    File repair_annotated_h5ad="${outstem}.h5ad"

  }

  runtime {
    docker: "gcr.io/broad-getzlab-fmhcrsparc/hcrseq:v1.0"
    memory: "16GB"
    disks: "local-disk " + disk_GB + " HDD"
    cpu: "4"
    preemptible: "2"
  }


}
