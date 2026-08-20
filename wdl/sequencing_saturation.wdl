workflow hcrseq_saturation{
  call calc_saturation

  output {
    File saturation_plot = calc_saturation.saturation_plot
    File saturation_data = calc_saturation.saturation_data
  }
}

task calc_saturation{
  File bam
  File cutadapt_json
  String outstem

  command {
    set -euo pipefail
	hcrseq amplicon calculate-sequencing-saturation \
        --bam=${bam} \
        --cutadapt_json=${cutadapt_json} \
        --outstem=${outstem}
                           
  }


  output {
    File saturation_plot="${outstem}.saturation_curve.pdf"
    File saturation_data="${outstem}.saturation_curve.csv"
  }

  runtime {
    docker: "gcr.io/broad-getzlab-fmhcrsparc/hcrseq:v1.0"
    memory: "4GB"
    disks: "local-disk 50 HDD"
    cpu: "1"
    preemptible: "2"
  }


}