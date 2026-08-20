workflow hcrseq_preprocess{
  call preprocess

  output {
    File bam = preprocess.bam
    File bai = preprocess.bai
    File cutadapt_json = preprocess.cutadapt_json
    File reporter_counts = preprocess.reporter_counts
    File repair_measurements = preprocess.repair_measurements
    File deletions = preprocess.deletions
    File insertions = preprocess.insertions
    File mismatch_dist = preprocess.mismatch_dist
  }
}

task preprocess{
  File fq1
  File fq2
  File hcrseq_ref
  String? foward_primer = "GACAACCACTACCTGAG"
  String? reverse_primer = "TCACTTGTACAGCTCGTCCATGC"
  String outstem
  Int? min_mapq = 5

  command {
    set -euo pipefail

    hcrseq amplicon prepare-reference \
                  --hcrseq_ref=${hcrseq_ref}
		  --forward_primer=${foward_primer} \
	          --reverse_primer=${reverse_primer} \
	          --outstem=amplicon_reference 

    hcrseq amplicon preprocess --fq1=${fq1} \
    		  --fq2=${fq2} \
                  --reference=amplicon_reference.fasta \
                  --primer_config amplicon_reference.primers.sh \
                  --outstem=${outstem}

    hcrseq amplicon quantify --bam=${outstem}.bam \
                --ref_path=amplicon_reference.pkl \
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
