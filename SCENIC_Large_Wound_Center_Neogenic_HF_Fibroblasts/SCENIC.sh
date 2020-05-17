
#! /usr/bin/env bash

#BSUB -J S_LW_Epi
#BSUB -n 56
#BSUB -R "span[hosts=1]"
#BSUB -W 999:59
#BSUB -o Synergy_SCENIC%J.out
#BSUB -e Synergy_SCENIC%J.err
#BSUB -N

source activate scenic
R CMD BATCH SCENIC_vis.R
