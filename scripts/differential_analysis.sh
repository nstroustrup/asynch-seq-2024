#!/bin/bash
#$ -l virtual_free=32G,h_rt=6:00:00
#$ -N diff_a   #job name
#$ -q short-sl7  #queue
#$ -t 1-68
Rscript /users/nstroustrup/nstroustrup/projects/asynch-seq-2022/scripts/differential_analysis.R \
  --name single_and_pooled_worms \
  --input ./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds \
  --models ./data/differential_analysis/models/single_and_pooled_worms.csv \
  --output_dir ./data/differential_analysis/tables/single_and_pooled_worms/ \
  --output_dds_dir ./data/differential_analysis/deseq/single_and_pooled_worms \
  --node_id $SGE_TASK_ID 
