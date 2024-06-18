#!/bin/bash
#$ -l virtual_free=16G,h_rt=72:00:00
#$ -N net_com   #job name
#$ -q long-sl7  #queue
Rscript /users/nstroustrup/nstroustrup/projects/asynch-seq-2022/scripts/network_and_communities.R --counts_and_annots "./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds" --tissue_corrected_counts "./data/tissue_specific_regression/rds/single_and_pooled_worms.rds" --counts_name "n2_d8_sw" --condition1_name "n2" --condition1_subset 'Day == 8 & NumberOfWorms == 1 & Strain == "QZ0" & Temperature==20 & Food=="NEC937"' --output_rds "./data/gene_network/all_genes_gene_network_d8.rds" --output_csv "./data/gene_network/all_genes_communities_d8.csv.gz" --output_covariance "./data/gene_network/all_genes_communities_d8_compact_file.rds" --min_community_size 20 --theta .6 --beta 2 --min_mu 30 --trials_infomap 200
