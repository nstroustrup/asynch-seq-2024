#!/bin/bash
#$ -l virtual_free=32G,h_rt=72:00:00
#$ -N TSR_net   #job name
#$ -q long-sl7  #queue
#$ -pe smp 7
Rscript /users/nstroustrup/nstroustrup/projects/asynch-seq-2022/scripts/tissue_specific_regression.R --input "./data/formated_counts/counts_and_annots_net_validation.rds" --name "net_validation" --output_dir "./data/tissue_specific_regression/tables/net_validation_2/" --specific_model "" --output_rds "./data/tissue_specific_regression/rds/net_validation.rds" --models "./data/tissue_specific_regression/models/net_validation.csv" --core_count 1 --exclude_ts T
