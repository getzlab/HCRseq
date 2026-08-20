workflow hcrseq_aggregate{
  call aggregate

  output {
    File repair_measurements = aggregate.repair_measurements
    File reporter_counts = aggregate.reporter_counts
    File qc_metrics = aggregate.qc_metrics
  }
}

task aggregate{
    Array[File] cutadapt_files
    Array[File] count_files
    Array[File] repair_files
    Array[String] sample_ids
    String outstem

    command {
        set -euo pipefail
	    
        hcrseq amplicon aggregate \
                                --cutadapt_files ${sep=',' cutadapt_files} \
    							--count_files ${sep=',' count_files} \
                                --repair_files ${sep=',' repair_files} \
                                --ids ${sep=',' sample_ids} \
                                --outstem=${outstem}
    

  }
  output {
    File repair_measurements="${outstem}.repair_measurements.csv"
    File reporter_counts="${outstem}.reporter_counts.csv"
    File qc_metrics="${outstem}.qc_metrics.csv"
  }

  runtime {
    docker: "gcr.io/broad-getzlab-fmhcrsparc/hcrseq:v1.0"
    memory: "4GB"
    disks: "local-disk 50 HDD"
    cpu: "1"
    preemptible: "2"
  }


}