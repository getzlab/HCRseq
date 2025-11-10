
hcrseq amplicon preprocess --fq1 amplicon_test_data/Experiment7_TK6_WT_2h_subset_R1.fastq.gz \
				--fq2 amplicon_test_data/Experiment7_TK6_WT_2h_subset_R2.fastq.gz \
				--reference ../reference/HCRseq_v0.3/amplicon/HCRseq_v0.3.fasta \
				--primer_config ../reference/HCRseq_v0.3/amplicon/HCRseq_v0.3.primers.sh \
				--outstem test_results/amplicon_preprocess 

hcrseq amplicon quantify --bam test_results/amplicon_preprocess.bam \
			--reference ../reference/HCRseq_v0.3/amplicon/HCRseq_v0.3.fasta \
			--lesion_info ../reference/HCRseq_v0.3/amplicon/HCRseq_v0.3.lesions.txt \
			--pathway_info ../reference/HCRseq_v0.3/pathway_calculations.txt \
			--outstem test_results/amplicon_quantify
