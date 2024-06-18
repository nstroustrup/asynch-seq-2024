#!/bin/bash
#$ -l virtual_free=32G,h_rt=6:00:00
#$ -N mars_sep  #job name
#$ -q short-sl7  #queue
#$ -t 1-500
Rscript /users/nstroustrup/nstroustrup/projects/asynch-seq-2022/scripts/gene_variability_single_worms.R -n 250 -i $SGE_TASK_ID -r F -c F --split_bootstraps_across_nodes 500