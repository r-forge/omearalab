#!/bin/bash
#$ -cwd
#$ -q long*

#$ -M omeara.brian@gmail.com
#$ -m beas
#$ -N fourstoch
#$ -r y
#$ -t 1-10

module load R/2.13.0
/data/apps/R/2.13.0/bin/R CMD BATCH ../../UnifiedApproachScripts/V4_StochasticFourRateCommands.R
