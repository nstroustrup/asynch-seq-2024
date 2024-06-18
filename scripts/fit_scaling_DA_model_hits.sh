#!/bin/bash
#$ -l virtual_free=16G,h_rt=10:00:00
#$ -N scal4   #job name
#$ -q long-sl7  #queue
#$ -t 1-240
#Rscript /users/nstroustrup/nstroustrup/projects/asynch-seq-2022/scripts/fit_scaling_DA_model.r -i $SGE_TASK_ID -n 250 --restricted_subset T -g scaling_hits_grid -G "data/scaling_models/genes_that_correlate_with_scale_factor_top.csv"
Rscript /users/nstroustrup/nstroustrup/projects/asynch-seq-2022/scripts/fit_scaling_DA_model.r -i $SGE_TASK_ID -n 250 --restricted_subset T -g all_30min --min_min_gene_count 30 --min_mean_gene_count 0 --print_count F
