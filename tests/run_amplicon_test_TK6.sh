
hcrseq amplicon preprocess --fq1 amplicon_test_data/LYG59_1_0462199279_05222025_0013_S3_L001_R1_001.downsample.fastq.gz \
				--fq2 amplicon_test_data/LYG59_1_0462199279_05222025_0013_S3_L001_R2_001.downsample.fastq.gz \
				--reference ../reference/HCRseq_v0.3/amplicon/HCRseq_v0.3.fasta \
				--primer_config ../reference/HCRseq_v0.3/amplicon/HCRseq_v0.3.primers.sh \
				--outstem test_results/amplicon_preprocessTK6_4h 

hcrseq amplicon quantify --bam test_results/amplicon_preprocessTK6_4h.bam \
			--reference ../reference/HCRseq_v0.3/amplicon/HCRseq_v0.3.fasta \
			--lesion_info ../reference/HCRseq_v0.3/amplicon/HCRseq_v0.3.lesions.txt \
			--pathway_info ../reference/HCRseq_v0.3/pathway_calculations.txt \
			--outstem test_results/amplicon_quantifyTK6_4h

hcrseq amplicon calculate-sequencing-saturation \
			--bam test_results/amplicon_preprocessTK6_4h.bam \ 
			--cutadapt_json test_results/amplicon_preprocessTK6_4h.cutadapt_metrics.json \
			--outstem test_results/amplicon_saturationTK6_4h
