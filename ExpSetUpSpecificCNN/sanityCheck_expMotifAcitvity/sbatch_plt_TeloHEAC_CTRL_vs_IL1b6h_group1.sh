#! /bin/bash
### Submit this Script with: sbatch sbatch_Train_CNN_TM.sh ###
 
# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition debug:
#SBATCH --partition=short
#  Use one node:
#SBATCH --nodes=1
#  Request 4 cores (hard constraint):
#SBATCH -c 4
#  Request 50GB of memory (hard constraint):
#SBATCH --mem=20GB
#  Request one hour maximal execution time (hard constraint):
#SBATCH --time=0-00:20:00
#SBATCH --job-name=plot


# Initialize the module system:sbat  
source /etc/profile.d/modules.sh
mkdir -p /fast/users/$USER/scratch/tmp
export TMPDIR=/fast/users/$USER/scratch/tmp
#command: bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/
# prepare seqs
python ./ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/sanityCheck_extract_Sequences.py --seqFile starrseq-all-final-toorder_oligocomposition.csv --output all_seqs.fa

# TelloHEAC_CTRL vs 6h IL1b
## prepare motif PWM
python ./ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/sanityCheck_extract_motifPWM.py \
--motif MA1130.1_FOSL2::JUN \
--motif MA1141.1_FOS::JUND \
--motif MA1128.1_FOSL1::JUN \
--output test.txt

## run fimo
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/FIMO.sh test.txt all_seqs.fa

## plot
python ./ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plot_experimental_activities.py --setUps TeloHAEC_CTRL --setUps TeloHAEC_IL1b_6h
