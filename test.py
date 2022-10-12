import torch
# print(torch.cuda.get_device_capability())
import os
import subprocess
#3070 capabilities    8.6
# fasta_path=r'/home/qin/source_code/test.fa'
# result_path=r'/home/qin/source_code/test'
#
# cmd='colabfold_batch --use-gpu-relax --num-recycle 3 --model-type AlphaFold2-multimer-v2 '+fasta_path+' '+result_path
#
# os.system(cmd,shell=True)
# # subprocess.run(cmd,shell=True)

from enum import Enum


class Metrics(Enum):
    SASA=.0
    energy=.0
    ratio_AKT=.0
    ratio_AKT_domain=0.


print(Weekday.energy)