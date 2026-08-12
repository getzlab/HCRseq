workflow hcrseq_preprocess{
  call preprocess
}

task preprocess{
  File fq1
  File fq2
  File ref
  File primer_config
  File ref_path
  String outstem
  Int? min_mapq = 5

  command {
    set -euo pipefail

    bwa index ${ref}
    samtools dict ${ref} > $(sed -E 's/fa(sta)?$/dict/' <<< ${ref})

    hcrseq amplicon preprocess --fq1=${fq1} \
    							--fq2=${fq2} \
                  --reference=${ref} \
                  --primer_config ${primer_config} \
                  --outstem=${outstem}

    hcrseq amplicon quantify --bam=${outstem}.bam \
                --ref_path=${ref_path} \
                --outstem=${outstem} \
                --min_mapq=${min_mapq}

  }


  output {
    File bam="${outstem}.bam"
    File bai="${outstem}.bai"
    File cutadapt_json="${outstem}.cutadapt_metrics.json"
    File reporter_counts="${outstem}.reporter_metrics.yaml"
    File repair_measurements="${outstem}.repair_measurements.yaml"
    File deletions="${outstem}.deletions.csv"
    File insertions="${outstem}.insertions.csv"
    File mismatch_dist="${outstem}.mismatch_dist.csv"
  }

  runtime {
    docker: "gcr.io/broad-getzlab-fmhcrsparc/hcrseq:v1.0"
    memory: "4GB"
    disks: "local-disk 50 HDD"
    cpu: "1"
    preemptible: "2"
  }


}