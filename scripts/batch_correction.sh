#!/bin/bash
#$ -l virtual_free=16G,h_rt=24:00:00
#$ -N batch_c   #job name
#$ -q long-sl7  #queue
Rscript /users/nstroustrup/nstroustrup/projects/asynch-seq-2022/scripts/batch_correction.R  \
  --name single_and_pooled_worms \
  --deseq_dir ./data/differential_analysis/deseq/single_and_pooled_worms \
  --counts ./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds \
  --models ./data/differential_analysis/models/single_and_pooled_worms.csv \
  --tissue_normalization \
  --output ./data/batch_corrected_counts/single_and_pooled_worms.rds
