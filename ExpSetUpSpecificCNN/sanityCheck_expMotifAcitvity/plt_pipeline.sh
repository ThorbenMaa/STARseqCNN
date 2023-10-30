#!/bin/bash
#command: bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline.sh
source /etc/profile.d/modules.sh

motifs=$(echo "\
--motif MA0107.1_RELA \
--motif MA0101.1_REL \
")

setUps=$(echo "\
--setUps TeloHAEC_CTRL \
--setUps TeloHAEC_IL1b_6h \
")

group=$(echo "--output group2")




# prepare seqs
python ./ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/sanityCheck_extract_Sequences.py --seqFile starrseq-all-final-toorder_oligocomposition.csv --output all_seqs.fa


# prepare motif PWM
python ./ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/sanityCheck_extract_motifPWM.py $1 --output test.txt

# run fimo
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/FIMO.sh test.txt all_seqs.fa

# plot
python ./ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plot_experimental_activities.py $2 $3
