#! /bin/bash
### Submit this Script with: sbatch abatch_ism_and_modisco.sh ###
 
# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition debug:
#SBATCH --partition=medium
#  Use one node:
#SBATCH --nodes=1
#  Request 4 cores (hard constraint):
#SBATCH -c 16
#  Request 50GB of memory (hard constraint):
#SBATCH --mem=50GB
#  Request one hour maximal execution time (hard constraint):
#SBATCH --time=0-10:00:00
#SBATCH --job-name=modisco

# Initialize the module system:sbat  
source /etc/profile.d/modules.sh
mkdir -p /fast/users/$USER/scratch/tmp
export TMPDIR=/fast/users/$USER/scratch/tmp

echo startInst
# mamba activate modisco-lite
python modisco_TM.py hypothetical_contribution_scores_mean_diffTeloHEAC_6h_vs_24h.npz Sequences.npz



modisco report -i modisco_resultshypothetical_contribution_scores_mean_diffTeloHEAC_6h_vs_24h.npz.h5 -o report_diffTeloHEAC_6h_vs_24h/ -s report_diffTeloHEAC_6h_vs_24h/ -m JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme_nice.txt


echo finished
