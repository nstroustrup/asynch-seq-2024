#!/bin/bash
#$ -l virtual_free=32G,h_rt=6:00:00
#$ -N fit_com   #job name
#$ -q short-sl7  #queue
#$ -t 1-28
Rscript /users/nstroustrup/nstroustrup/projects/asynch-seq-2022/scripts/fit_gene_coexpression_groups.R -i $SGE_TASK_ID -t T --da_model_comparison F