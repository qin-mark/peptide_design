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
#SBATCH --job-name=lambo
#SBATCH --output=lambo_results-%j.log
echo "### Starting at: $(date) ###"
nvidia-smi
python scripts/black_box_opt.py optimizer=lambo optimizer.encoder_obj=mlm task=proxy_rfp tokenizer=protein surrogate=multi_task_exact_gp acquisition=nehvi