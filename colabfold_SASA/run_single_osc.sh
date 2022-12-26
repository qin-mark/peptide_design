#!/bin/bash
conda activate my_colabfold
fasta=$1
result=$2
mkdir $result
echo $fasta
echo $result
colabfold_batch --use-gpu-relax --num-recycle 3 --model-type AlphaFold2-ptm $fasta $result
