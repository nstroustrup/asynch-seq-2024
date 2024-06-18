#!/bin/bash
#$ -l virtual_free=32G,h_rt=6:00:00
#$ -N fit_com2   #job name
#$ -q short-sl7  #queue
Rscript /users/nstroustrup/nstroustrup/projects/asynch-seq-2022/scripts/fit_gene_coexpression_groups.R --compile_and_label T -t F 