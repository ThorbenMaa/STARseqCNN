#! /bin/bash
### Submit this Script with: sbatch sbatch_Train_CNN_TM.sh ###
 
# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition debug:
#SBATCH --partition=medium #gpu --gres=gpu:1
#  Use one node:
#SBATCH --nodes=1
#  Request 4 cores (hard constraint):
#SBATCH -c 16
#  Request 50GB of memory (hard constraint):
#SBATCH --mem=50GB
#  Request one hour maximal execution time (hard constraint):
#SBATCH --time=0-05:00:00
#SBATCH --job-name=TrainAndTest


# Initialize the module system:sbat  
source /etc/profile.d/modules.sh
mkdir -p /fast/users/$USER/scratch/tmp
export TMPDIR=/fast/users/$USER/scratch/tmp

echo hepG2VStelo_ctrl_noaug_normStattlog
echo startInst
# mamba activate CNN_TM
python ExpSetUpSpecificCNN/train_or_eval_CNNs_setUpSpec_norm.py \
--counts1 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo1.txt  \
--counts2 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo2.txt  \
--seqs starrseq-all-final-toorder_oligocomposition.csv  \
--mode train  \
--holdOut chr8  \
--useAug no_aug  \
--model noModel \
--compare mean_HepG2_untreatedPilot \
--compare mean_TeloHAEC_CTRL

python ExpSetUpSpecificCNN/train_or_eval_CNNs_setUpSpec_norm.py \
--counts1 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo1.txt  \
--counts2 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo2.txt  \
--seqs starrseq-all-final-toorder_oligocomposition.csv  \
--mode load  \
--holdOut chr8  \
--useAug no_aug  \
--model allseq-CNN_StarSeq_model_Minna_deepSTAR_lr0.01no_augsetUpSpecNorm\(\'mean_HepG2_untreatedPilot\',\ \'mean_TeloHAEC_CTRL\'\)   \
--compare mean_HepG2_untreatedPilot \
--compare mean_TeloHAEC_CTRL
   

echo finished
