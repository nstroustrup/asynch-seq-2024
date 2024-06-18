#!/bin/bash
#$ -l virtual_free=16G,h_rt=72:00:00
#$ -N TSR_wt   #job name
#$ -q long-sl7  #queue
#$ -pe smp 7
Rscript /users/nstroustrup/nstroustrup/projects/asynch-seq-2022/scripts/tissue_specific_regression.R --input "./data/formated_counts/counts_and_annots_single_and_pooled_worms.rds" --name "single_and_pooled_worms" --output_dir "./data/tissue_specific_regression/tables/single_and_pooled_worms/" --models "./data/tissue_specific_regression/models/single_and_pooled_worms.csv" --specific_model "n2_d8_sw" --output_rds "./data/tissue_specific_regression/rds/single_and_pooled_worms.rds" --core_count 7
