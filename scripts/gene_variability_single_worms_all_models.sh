#!/bin/bash
#$ -l virtual_free=32G,h_rt=72:00:00
#$ -N ns_sw_FF   #job name
#$ -q long-sl7  #queue
#$ -t 1-52
Rscript /users/nstroustrup/nstroustrup/projects/asynch-seq-2022/scripts/gene_variability_single_worms.R -n 250 -i $SGE_TASK_ID -r F -c F
