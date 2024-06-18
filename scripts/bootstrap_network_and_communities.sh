NBOOT=10
mkdir -p ./data/gene_network/temp/
mkdir -p ./data/gene_network/bootstrap/

for ITER in {1..10}
do
	echo "Starting iteration ${ITER}"
	UNIQUE_NAME=$(cat /dev/urandom | tr -dc '[:alpha:]' | fold -w ${1:-30} | head -n 1)
	echo "Formatting Counts"
	Rscript ./scripts/format_counts.R \
	  --name single_worms_bootstrap \
	  --counts ./data/annotated_counts/counts_single_worms_bootstrap.csv.gz \
	  --annotations ./data/annotations/sample_annotations_single_worms_bootstrap.csv \
	  --output ./data/gene_network/temp/counts_and_annots_${UNIQUE_NAME}.rds \
	  --bootstrap
#echo "Differential Analysis"
#Rscript ./scripts/differential_analysis.R \
#  --name single_worms_bootstrap \
#  --input ./data/gene_network/temp/counts_and_annots_${UNIQUE_NAME}.rds\
#  --models ./data/differential_analysis/models/single_worms_bootstrap.csv \
#  --output_dir ./data/gene_network/temp/differential_analysis/tables/ \
#  --output_suffix _${UNIQUE_NAME} \
#  --output_dds_dir ./data/gene_network/temp/differential_analysis/deseq
#echo "Batch Correction"
#Rscript ./scripts/batch_correction.R  \
#  --name single_worms_bootstrap \
#  --input_suffix _${UNIQUE_NAME} \
#  --counts ./data/gene_network/temp/counts_and_annots_${UNIQUE_NAME}.rds \
#  --models ./data/differential_analysis/models/single_worms_bootstrap.csv \
#  --deseq_dir ./data/gene_network/temp/differential_analysis/deseq \
#  --output  ./data/gene_network/temp/batch_corrected_${UNIQUE_NAME}.rds

	echo "Tissue_specific Regression"
	Rscript ./scripts/tissue_specific_regression.R \
	  --name single_worm_bootstrap \
	  --input ./data/gene_network/temp/counts_and_annots_${UNIQUE_NAME}.rds \
	  --models ./data/tissue_specific_regression/models/single_worms_bootstrap.csv \
	  --output_dir ./data/gene_network/temp/tissue_specific_regression/tables \
	  --output_suffix _${UNIQUE_NAME} \
	  --specific_model "n2_d8_sw" \
	  --output_rds ./data/gene_network/temp/tissue_specific_regression/rds/rds_${UNIQUE_NAME}.rds \
	  --core_count 1

	echo "Network and Communities"
	Rscript ./scripts/network_and_communities.R \
	 --counts_and_annots ./data/gene_network/temp/counts_and_annots_${UNIQUE_NAME}.rds \
	 --tissue_corrected_counts ./data/gene_network/temp/tissue_specific_regression/rds/rds_${UNIQUE_NAME}.rds \
	 --output_csv ./data/gene_network/bootstrap/${UNIQUE_NAME}.csv.gz \
	 --counts_name n2_d8_sw \
	 --condition1_name n2 \
	 --condition1_subset 'Day == 8 & NumberOfWorms == 1 & Genotype == "N2"' \
	 --min_community_size 20 \
	 --theta .6 \
	 --beta 2 \
	 --min_mu 30 \
	 --trials_infomap 200
	
	echo "Deleting temporary files"
	rm -f ./data/gene_network/temp/counts_and_annots_${UNIQUE_NAME}.rds
	rm -f ./data/gene_network/temp/differential_analysis/tables/*${UNIQUE_NAME}*
	rm -f ./data/gene_network/temp/differential_analysis/deseq/*${UNIQUE_NAME}*
	rm -f ./data/gene_network/temp/tissue_specific_regression/tables/*${UNIQUE_NAME}*
	rm -f ./data/gene_network/temp/tissue_specific_regression/rds/rds_${UNIQUE_NAME}.rds
	echo "Iteration complete"
done
