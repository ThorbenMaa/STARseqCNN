#! /bin/bash
### Submit this Script with: sbatch sbatch_Train_CNN_TM.sh ###
 
# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition debug:
#SBATCH --partition=gpu --gres=gpu:1
#  Use one node:
#SBATCH --nodes=1
#  Request 4 cores (hard constraint):
#SBATCH -c 4
#  Request 50GB of memory (hard constraint):
#SBATCH --mem=50GB
#  Request one hour maximal execution time (hard constraint):
#SBATCH --time=0-4:00:00
#SBATCH --job-name=TrainAndTest


# Initialize the module system:sbat  
source /etc/profile.d/modules.sh
mkdir -p /fast/users/$USER/scratch/tmp
export TMPDIR=/fast/users/$USER/scratch/tmp


echo startInst
# mamba activate CNN_TM
python model_train_eval_interpretation/trainCNN.py \
--activityFile 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo1.txt \
--activityFile 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo2.txt \
--seqFile starrseq-all-final-toorder_oligocomposition.csv \
--holdOutChr chr8 \
--earlyStop True
--modelName TM_earlyStop

   

echo finished
