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
python model_train_eval_interpretation/modisco_TM.py \
--seq_file TM_s_model_early_stopSequences.npz \
--contrib_file TM_s_model_early_stophypothetical_contribution_scores_mean_HASMC_Chol.npz \
--contrib_file TM_s_model_early_stophypothetical_contribution_scores_mean_HepG2_untreatedPilot.npz \
--contrib_file TM_s_model_early_stophypothetical_contribution_scores_mean_HASMC_untreatedPilot.npz \
--contrib_file TM_s_model_early_stophypothetical_contribution_scores_mean_TeloHAEC_IL1b_24h.npz \
--contrib_file TM_s_model_early_stophypothetical_contribution_scores_mean_TeloHAEC_IL1b_6h.npz \
--contrib_file TM_s_model_early_stophypothetical_contribution_scores_mean_RAW_TGFB.npz \
--contrib_file TM_s_model_early_stophypothetical_contribution_scores_mean_TeloHAEC_CTRL.npz \
--contrib_file TM_s_model_early_stophypothetical_contribution_scores_mean_RAW_IL1B.npz \
--contrib_file TM_s_model_early_stophypothetical_contribution_scores_mean_RAW_CTRL.npz \
--contrib_file TM_s_model_early_stophypothetical_contribution_scores_mean_ccell_3T3_undiff_CTRL.npz \
--contrib_file TM_s_model_early_stophypothetical_contribution_scores_mean_cell_3T3_undiff_TGFB.npz \
--contrib_file TM_s_model_early_stophypothetical_contribution_scores_mean_cell_3T3_diff_CTRL.npz \


modisco report -i modisco_TM_s_model_early_stophypothetical_contribution_scores_mean_HASMC_Chol.npz.h5 -o report_TM_early_stop_HASMC_chol/ -s report_TM_early_stop_HASMC_Chol/ -m JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme_nice.txt
modisco report -i modisco_TM_s_model_early_stophypothetical_contribution_scores_mean_HepG2_untreatedPilot.npz.h5 -o report_TM_early_stop_HepG2_untreatedPilot/ -s report_TM_early_stop_HepG2_untreatedPilot/ -m JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme_nice.txt
modisco report -i modisco_TM_s_model_early_stophypothetical_contribution_scores_mean_HASMC_untreatedPilot.npz.h5 -o report_TM_early_stop_HASMC_untreatedPilot/ -s report_TM_early_stop_HASMC_untreatedPilot/ -m JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme_nice.txt
modisco report -i modisco_TM_s_model_early_stophypothetical_contribution_scores_mean_TeloHAEC_IL1b_24h.npz.h5 -o report_TM_early_stop_TeloHAEC_IL1b_24h/ -s report_TM_early_stop_TeloHAEC_IL1b_24h/ -m JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme_nice.txt
modisco report -i modisco_TM_s_model_early_stophypothetical_contribution_scores_mean_TeloHAEC_IL1b_6h.npz.h5 -o report_TM_early_stop_TeloHAEC_IL1b_6h/ -s report_TM_early_stop_TeloHAEC_IL1b_6h/ -m JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme_nice.txt
modisco report -i modisco_TM_s_model_early_stophypothetical_contribution_scores_mean_RAW_TGFB.npz.h5 -o report_TM_early_stop_RAW_TGFB/ -s report_TM_early_stop_RAW_TGFB/ -m JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme_nice.txt
modisco report -i modisco_TM_s_model_early_stophypothetical_contribution_scores_mean_TeloHAEC_CTRL.npz.h5 -o report_TM_early_stop_TeloHAEC_CTRL/ -s report_TM_early_stop_TeloHAEC_CTRL/ -m JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme_nice.txt
modisco report -i modisco_TM_s_model_early_stophypothetical_contribution_scores_mean_RAW_IL1B.npz.h5 -o report_TM_early_stop_RAW_IL1B/ -s report_TM_early_stop_RAW_IL1B/ -m JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme_nice.txt
modisco report -i modisco_TM_s_model_early_stophypothetical_contribution_scores_mean_RAW_CTRL.npz.h5 -o report_TM_early_stop_RAW_CTRL/ -s report_TM_early_stop_RAW_CTRL/ -m JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme_nice.txt
modisco report -i modisco_TM_s_model_early_stophypothetical_contribution_scores_mean_ccell_3T3_undiff_CTRL.npz.h5 -o report_TM_early_stop_ccell_3T3_undiff_CTRL/ -s report_TM_early_stop_ccell_3T3_undiff_CTRL/ -m JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme_nice.txt
modisco report -i modisco_TM_s_model_early_stophypothetical_contribution_scores_mean_cell_3T3_undiff_TGFB.npz.h5 -o report_TM_early_stop_cell_3T3_undiff_TGFB/ -s report_TM_early_stop_cell_3T3_undiff_TGFB/ -m JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme_nice.txt
modisco report -i modisco_TM_s_model_early_stophypothetical_contribution_scores_mean_cell_3T3_diff_CTRL.npz.h5 -o report_TM_early_stop_cell_3T3_diff_CTRL/ -s report_TM_early_stop_cell_3T3_diff_CTRL/ -m JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme_nice.txt

echo finished
