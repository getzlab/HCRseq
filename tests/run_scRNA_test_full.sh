hcrseq scrna postprocess-reporter-bam \
			--bam scRNA_test_data/P94x_24h.bam \
			--ref_path ../reference/HCRseq_v0.4.1/HCRseq_v0.4.1.pkl\
			--outstem scRNA_test_results/scRNA_test_full


hcrseq scrna quantify --h5_file scRNA_test_data/experiment11_2h_matrix.h5 \
			--bam scRNA_test_results/scRNA_test_full_postprocessed.sorted.bam \
			--ref_path ../reference/HCRseq_v0.4.1/HCRseq_v0.4.1.pkl\
			--outstem scRNA_test_results/scRNA_test_full

