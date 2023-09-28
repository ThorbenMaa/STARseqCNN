#! /bin/bash
### Submit this Script with: sbatch script.sh ###
 
# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition debug:
#SBATCH --partition=shortterm
#  Use one node:
#SBATCH --nodes=1
#  Request 4 cores (hard constraint):
#SBATCH -c 4
#  Request GPU SBATCH --gres=gpu:1
#  Request 50GB of memory (hard constraint):
#SBATCH --mem=50GB
#  Request one hour maximal execution time (hard constraint):
#SBATCH --time=0-04:00:00
#  Request 100 GB of local scratch disk (hard constraint):
#SBATCH --tmp=1GB
#  Notify me at job start and end:
#SBATCH --mail-type=ALL
#  Send the notifications to (change to your email address!):
#SBATCH --mail-user=tho.maass@uni-luebeck.de
#  Find your job easier with a name:
#SBATCH --job-name=sanityCheckSNV


# Initialize the module system:
source /etc/profile.d/modules.sh




echo startInst
# mamba activate CNN_TM
python sanityCheck_variantEffects.py 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo1.txt starrseq-all-final-toorder_oligocomposition.csv 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo2.txt allseq-CNN_StarSeq_model_Minna_deepSTAR_lr0.01no_aug chr8 all
                    

echo finished
