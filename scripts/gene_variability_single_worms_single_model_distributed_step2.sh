#!/bin/bash
#$ -l virtual_free=64G,h_rt=12:00:00
#$ -N n2_sw2   #job name
#$ -q long-sl7  #queue
Rscript /users/nstroustrup/nstroustrup/projects/asynch-seq-2022/scripts/gene_variability_single_worms.R -n 250 -m n2_sw -r T -c T --split_bootstraps_across_nodes 50

