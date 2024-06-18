NBOOT=10
mkdir -p ./data/gene_network/temp/
mkdir -p ./data/gene_network/bootstrap/

for ITER in {1..$NBOOT}; do

	UNIQUE_NAME=$(cat /dev/urandom | tr -dc '[:alpha:]' | fold -w ${1:-30} | head -n 1)

	Rscript ./scripts/format_counts.R \
	  --name single_worms_bootstrap \
	  --counts ./data/annotated_counts/counts_single_worms_bootstrap.csv.gz \
	  --annotations ./data/annotations/sample_annotations_single_worms_bootstrap.csv \
	  --output ./data/gene_network/temp/counts_and_annots_${UNIQUE_NAME}.rds \
	  --bootstrap

	Rscript ./scripts/differential_analysis.R \
	  --name single_worms_bootstrap \
	  --input ./data/gene_network/temp/counts_and_annots_${UNIQUE_NAME}.rds\
	  --models ./data/differential_analysis/models/single_worms_bootstrap.csv \
	  --output_dds ./data/gene_network/temp/dds_${UNIQUE_NAME}.rds

	Rscript ./scripts/batch_correction.R  \
	  --name single_worms_bootstrap \
	  --input ./data/gene_network/temp/dds_${UNIQUE_NAME}.rds \
	  --models ./data/differential_analysis/models/single_worms_bootstrap.csv \
	  --output  ./data/gene_network/temp/batch_corrected_${UNIQUE_NAME}.rds

	Rscript ./scripts/network_and_communities.R \
	 --counts_and_annots ./data/gene_network/temp/counts_and_annots_${UNIQUE_NAME}.rds \
	 --batch_corrected_counts ./data/gene_network/temp/batch_corrected_${UNIQUE_NAME}.rds \
	 --output_csv ./data/gene_network/bootstrap/${UNIQUE_NAME}.csv.gz \
	 --counts_name d8_sw \
	 --condition1_name n2 \
	 --condition2_name daf2 \
	 --condition1_subset 'Day == 8 & NumberOfWorms == 1 & Genotype == "N2"' \
	 --condition2_subset 'Day == 8 & NumberOfWorms == 1 & Genotype == "daf-2(e1368)"'

	rm -f ./data/gene_network/temp/counts_and_annots_${UNIQUE_NAME}.rds
	rm -f ./data/gene_network/temp/dds_${UNIQUE_NAME}.rds
	rm -f ./data/gene_network/temp/batch_corrected_${UNIQUE_NAME}.rds
done
