#! /bin/bash
### Submit this Script with: sbatch sbatch_ism.sh ###
 
# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition debug:
#SBATCH --partition=medium #gpu --gres=gpu:1
#  Use one node:
#SBATCH --nodes=1
#  Request 4 cores (hard constraint):
#SBATCH -c 64
#  Request 50GB of memory (hard constraint):
#SBATCH --mem=50GB
#  Request one hour maximal execution time (hard constraint):
#SBATCH --time=0-60:00:00
#SBATCH --job-name=ism


# Initialize the module system:sbat  
source /etc/profile.d/modules.sh
mkdir -p /fast/users/$USER/scratch/tmp
export TMPDIR=/fast/users/$USER/scratch/tmp

echo startInstism_TeloHEAC_CTRL_vs_6h
# mamba activate CNN_TM
python ./ExpSetUpSpecificCNN/ism_spec_TM.py 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo1.txt starrseq-all-final-toorder_oligocomposition.csv allseq-CNN_StarSeq_model_Minna_deepSTAR_lr0.01use_augsetUpSpec('mean_TeloHAEC_CTRL', 'mean_TeloHAEC_IL1b_6h') 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo2.txt TeloHEAC_CTRL_vs_6h

echo finished
# allseq-CNN_StarSeq_model_Minna_deepSTAR_lr0.01use_augsetUpSpec('mean_TeloHAEC_CTRL', 'mean_TeloHAEC_IL1b_6h')