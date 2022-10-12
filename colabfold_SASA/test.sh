#!/bin/bash
#-------------------------------------------------------------------------------
#  SBATCH CONFIG
#-------------------------------------------------------------------------------
## resources
#SBATCH -A PCON0022
#SBATCH --gpus-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --time 3-00:00:00
## labels and outputs
#SBATCH --job-name=test
#SBATCH --output=test_results-%j.log
echo "### Starting at: $(date) ###"
colabfold_batch --templates --use-gpu-relax --num-recycle 3 --model-type AlphaFold2-multimer-v2 ./fasta/QNGGIHRLCQLLSKAHEMVQRRLALAPDA.fa ./result/QNGGIHRLCQLLSKAHEMVQRRLALAPDA
