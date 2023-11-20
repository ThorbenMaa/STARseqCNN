#!/bin/bash

#cd /data/gpfs-1/work/groups/ag_kircher/Kaikkonen_2023/filesFromMinnaKaikkonen/github_repo/STARseqCNN/

# train CNN
#mamba activate CNN_TM

#sbatch_Train_CNN_setUpSpec_TM.sh




# model interpretation
#sbatch sbatch_ism_spec.sh

#mamba activate modisco_lite_v3

#sbatch sbatch_tfmodisco_spec.sh


variants="MPRAlm_Alec"
setUpComp="diffTeloHEAC_CTRL_vs_6h"
setUpOne="TeloHAEC_CTRL"
setUpTwo="TeloHAEC_IL1b_6h"

# analyse motifs
export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.4:$PATH

python ./allMotifsWithSignificantEffects/trimMotifs.py \
--outputPWM ./results/${setUpComp}/${setUpComp}_PWMs.txt \
--reportFile modisco_resultshypothetical_contribution_scores_mean_${setUpComp}.npz.h5


python ./allMotifsWithSignificantEffects/removeSimilarMotifs.py \
--PWM_file ./results/${setUpComp}/${setUpComp}_PWMs.txt \
--outputPWM ./results/${setUpComp}/condensed_${setUpComp}_PWMs.txt \
--tomtomSig 0.01
 
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/FIMO.sh \
./results/${setUpComp}/condensed_${setUpComp}_PWMs.txt \
all_seqs.fa \


python ./allMotifsWithSignificantEffects/plot_experimental_activities_perSetUp_comp.py \
--output ./results/${setUpComp}/${setUpComp} \
--outputPWMs ./results/${setUpComp}/${setUpComp}_sigPWMs.txt \
--setUps ${setUpOne} \
--setUps ${setUpTwo} \
--PWMfile ./results/${setUpComp}/condensed_${setUpComp}_PWMs.txt \
--tomtomSig 0.01 \
--number-motifs_to_plot 1 \
--fimoQual 0.2




# check varaint effects
python VariantEffects/extractVariantEffectsFromMatchFile.py \
--output ./results/${setUpComp}/${variants}_${setUpComp} \
--p_value 0.01 \
--matchFile TeloHAEC_CTRL-vs-TeloHAEC_IL1b_6h.csv #HARD CODED!

#run fimo sh
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/FIMO.sh \
./results/${setUpComp}/${setUpComp}_sigPWMs.txt \
./results/${setUpComp}/${variants}_${setUpComp}_seqs.fa \

#gives Number of variants that can be explained by overlapping with a motif (not same as of number of motifs that align with a variant!!!! This would be the output of plot_alignment... script)
python VariantEffects/check_variant_within_motif.py \
--variantsFile ./results/${setUpComp}/${variants}_${setUpComp}.csv \
--output ./results/${setUpComp}/${variants}_${setUpComp}_withMotifs \
--q_thres 0.2

python VariantEffects/plot_alignment_withPWM.py \
--aligMotifFile ./results/${setUpComp}/${variants}_${setUpComp}_withMotifs.csv \
--PWMsFile ./results/${setUpComp}/${setUpComp}_sigPWMs.txt \
--output ./results/${setUpComp}/plot_${variants}_${setUpComp} \
--chunk_size 2




#tbd: plotting script