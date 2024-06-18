## To replicate results, run the following commands.

```sh
# gene annotations
Rscript ./scripts/format_gene_annotations.R

# merge counts files and writes them to ./data/annotated_counts
Rscript ./scripts/merge_counts.R single_and_pooled_worms
Rscript ./scripts/merge_counts.R single_worms_bootstrap
Rscript ./scripts/merge_counts.R beads
Rscript ./scripts/merge_counts.R tso_only
Rscript ./scripts/merge_counts.R glp1_glp4
Rscript ./scripts/merge_counts.R tissue_specific_rnai
Rscript ./scripts/merge_counts.R net_validation
Rscript ./scripts/merge_counts.R 2022_02_batch

# genome coverage of beads
qsub -N genomecov -l virtual_free=30G -cwd \
  -q short-sl7 -e logs -o logs -t 1-13 \
  ./scripts/genomecov.sh

# format counts; no tissue-specific normalization
Rscript ./scripts/format_counts.R \
  --name beads \
  --counts ./data/annotated_counts/counts_beads.csv.gz \
  --annotations ./data/annotations/sample_annotations_beads.csv \
  --output ./data/formated_counts/counts_and_annots_beads.rds
  
Rscript ./scripts/format_counts.R \
  --name tso_only \
  --counts ./data/annotated_counts/counts_tso_only.csv.gz \
  --annotations ./data/annotations/sample_annotations_tso_only.csv \
  --output ./data/formated_counts/counts_and_annots_tso_only.rds
  
Rscript ./scripts/format_counts.R \
  --name glp1_glp4 \
  --counts ./data/annotated_counts/counts_glp1_glp4.csv.gz \
  --annotations ./data/annotations/sample_annotations_glp1_glp4.csv \
  --output ./data/formated_counts/counts_and_annots_glp1_glp4.rds
  
Rscript ./scripts/format_counts.R \
  --name tissue_specific_rnai \
  --counts ./data/annotated_counts/counts_tissue_specific_rnai.csv.gz \
  --annotations ./data/annotations/sample_annotations_tissue_specific_rnai.csv \
  --output ./data/formated_counts/counts_and_annots_tissue_specific_rnai.rds
  
# differential analysis
Rscript ./scripts/differential_analysis.R \
  --name beads \
  --input ./data/formated_counts/counts_and_annots_beads.rds \
  --models ./data/differential_analysis/models/beads.csv \
  --output_dir ./data/differential_analysis/tables/beads/

Rscript ./scripts/differential_analysis.R \
  --name tso_only \
  --input ./data/formated_counts/counts_and_annots_tso_only.rds \
  --models ./data/differential_analysis/models/tso_only.csv \
  --output_dir ./data/differential_analysis/tables/tso_only/

Rscript ./scripts/differential_analysis.R \
  --name glp1_glp4 \
  --input ./data/formated_counts/counts_and_annots_glp1_glp4.rds \
  --models ./data/differential_analysis/models/glp1_glp4.csv \
  --output_dir ./data/differential_analysis/tables/glp1_glp4/
  
Rscript ./scripts/differential_analysis.R \
  --name tissue_specific_rnai \
  --input ./data/formated_counts/counts_and_annots_tissue_specific_rnai.rds \
  --models ./data/differential_analysis/models/tissue_specific_rnai.csv \
  --output_dir ./data/differential_analysis/tables/tissue_specific_rnai/

# define genes unique to tissues
Rscript ./scripts/define_tissue_unique_genes.R

# format counts with and without tissue-specific normalization
Rscript ./scripts/format_counts.R \
  --name single_and_pooled_worms \
  --counts ./data/annotated_counts/counts_single_and_pooled_worms.csv.gz \
  --annotations ./data/annotations/sample_annotations_single_and_pooled_worms.csv \
  --output ./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds \
#  --tissue_normalization
  
Rscript ./scripts/format_counts.R \
  --name net_validation \
  --counts ./data/annotated_counts/counts_net_validation.csv.gz \
  --annotations ./data/annotations/sample_annotations_net_validation.csv \
  --output ./data/formated_counts/counts_and_annots_net_validation.rds \
#  --tissue_normalization

# differential analysis
  Rscript ./scripts/differential_analysis.R \
  --name net_validation \
  --input ./data/formated_counts/counts_and_annots_net_validation.rds \
  --models ./data/differential_analysis/models/net_validation.csv \
  --output_dir ./data/differential_analysis/tables/net_validation/ \
  --output_dds_dir ./data/differential_analysis/deseq/net_validation
  
Rscript ./scripts/differential_analysis.R \
  --name single_and_pooled_worms \
  --input ./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds \
  --models ./data/differential_analysis/models/single_and_pooled_worms.csv \
  --output_dir ./data/differential_analysis/tables/single_and_pooled_worms/ \
  --output_dds_dir ./data/differential_analysis/deseq/single_and_pooled_worms
  
# batch correction
Rscript ./scripts/batch_correction.R  \
  --name single_and_pooled_worms \
  --deseq_dir ./data/differential_analysis/deseq/single_and_pooled_worms \
  --counts ./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds \
  --models ./data/differential_analysis/models/single_and_pooled_worms.csv \
#  --tissue_normalization \
  --output ./data/batch_corrected_counts/single_and_pooled_worms.rds
  
# tissue-specific regression, normalization and batch correction for day 8 wild-type network
Rscript ./scripts/tissue_specific_regression.R \
  --name single_and_pooled_worms \
  --input ./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds \
  --models ./data/tissue_specific_regression/models/single_and_pooled_worms.csv \
  --output_dir ./data/tissue_specific_regression/tables/single_and_pooled_worms/ \
  --specific_model "n2_d8_sw" 
  --output_rds ./data/tissue_specific_regression/rds/single_and_pooled_worms.rds

# tissue-specific regression, normalization and batch correction for the wild-type aging timeseries
Rscript ./scripts/tissue_specific_regression.R \
  --input "./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds" \
  --name "single_and_pooled_worms_ts" \
  --output_dir "./data/tissue_specific_regression/tables/single_and_pooled_worms_ts/" \
  --models "./data/tissue_specific_regression/models/single_and_pooled_worms.csv" \
  --specific_model "n2_series" \
  --output_rds "./data/tissue_specific_regression/rds/single_and_pooled_worms_n2_series_all.rds" 
  
# tissue-specific regression, normalization and batch correction for the perturbation screen results
Rscript ./scripts/tissue_specific_regression.R \
  --name net_validation \
  --input ./data/formated_counts/counts_and_annots_net_validation.rds \
  --models ./data/tissue_specific_regression/models/net_validation.csv \
  --output_dir ./data/tissue_specific_regression/tables/net_validation/ \
  --output_rds ./data/tissue_specific_regression/rds/net_validation.rds

# gene variability
Rscript ./scripts/gene_variability_single_worms.R -m n2_sw -n 50
#OR
qsub < ./scripts/gene_variability_single_worms_all_models.sh
Rscript ./scripts/gene_variability_single_vs_pooled_worms.R

 # gene network inference 
Rscript ./scripts/network_and_communities.R \
--counts_and_annots "./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds" \
--tissue_corrected_counts "./data/tissue_specific_regression/rds/single_and_pooled_worms.rds" \
--counts_name "n2_d8_sw" \
--condition1_name "n2" \
--condition1_subset 'Day == 8 & NumberOfWorms == 1 & Strain == "QZ0" & Temperature==20 & Food=="NEC937"' \
--output_rds "./data/gene_network/all_genes_gene_network_d8.rds" \
--output_csv "./data/gene_network/all_genes_communities_d8.csv.gz" \
--output_covariance "./data/gene_network/all_genes_communities_d8_compact_file.rds" \
--min_community_size 20 \
--theta .6 \
--beta 2 \
--min_mu 30 \
--trials_infomap 200

# bootstrap gene network inference
qsub -N bootnet -l virtual_free=8G -cwd \
  -q short-sl7 -e logs -o logs -t 1-100 \
  ./scripts/bootstrap_network_and_communities.sh

```
