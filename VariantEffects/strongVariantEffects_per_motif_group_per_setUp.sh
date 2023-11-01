#!/bin/bash
#command: bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline_on_all.sh
setUp="TeloHAEC_IL1b_6h"
TFgroup="group1"

echo ${TFgroup}
# generate motifs in fimo format
python ./ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/sanityCheck_extract_motifPWM.py \
--motif MA1130.1_FOSL2::JUN \
--motif MA1141.1_FOS::JUND \
--motif MA1128.1_FOSL1::JUN \
--output VariantEffects/${setUp}_${TFgroup}_motifs.txt

# fimo on TeloHEAC_IL1b_6h top variant effects
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/FIMO.sh VariantEffects/${setUp}_${TFgroup}_motifs.txt VariantEffects/10percent__df_highestVariantEffects_mean_${setUp}.fa

# how many of the top variants effects can I explain with them lying in one of the TF motifs?
python VariantEffects/check_variant_within_motif.py \
--activityFile VariantEffects/df_highestVariantEffects_mean_${setUp}.csv \
--expSetUp ${setUp} \
--output VariantEffects/${setUp}_${TFgroup}





TFgroup="group2" #--motif auch ändern

echo ${TFgroup}
# generate motifs in fimo format
python ./ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/sanityCheck_extract_motifPWM.py \
--motif MA0107.1_RELA \
--motif MA0101.1_REL \
--output VariantEffects/${setUp}_${TFgroup}_motifs.txt

# fimo on TeloHEAC_IL1b_6h top variant effects
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/FIMO.sh VariantEffects/${setUp}_${TFgroup}_motifs.txt VariantEffects/10percent__df_highestVariantEffects_mean_${setUp}.fa

# how many of the top variants effects can I explain with them lying in one of the TF motifs?
python VariantEffects/check_variant_within_motif.py \
--activityFile VariantEffects/df_highestVariantEffects_mean_${setUp}.csv \
--expSetUp ${setUp} \
--output VariantEffects/${setUp}_${TFgroup}






TFgroup="group3" #--motif auch ändern

echo ${TFgroup}
# generate motifs in fimo format
python ./ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/sanityCheck_extract_motifPWM.py \
--motif MA0488.1_JUN \
--motif MA0018.4_CREB1 \
--motif MA0492.1_JUND \
--output VariantEffects/${setUp}_${TFgroup}_motifs.txt

# fimo on TeloHEAC_IL1b_6h top variant effects
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/FIMO.sh VariantEffects/${setUp}_${TFgroup}_motifs.txt VariantEffects/10percent__df_highestVariantEffects_mean_${setUp}.fa

# how many of the top variants effects can I explain with them lying in one of the TF motifs?
python VariantEffects/check_variant_within_motif.py \
--activityFile VariantEffects/df_highestVariantEffects_mean_${setUp}.csv \
--expSetUp ${setUp} \
--output VariantEffects/${setUp}_${TFgroup}



TFgroup="group4" #--motif auch ändern

echo ${TFgroup}
# generate motifs in fimo format
python ./ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/sanityCheck_extract_motifPWM.py \
--motif MA1633.2_BACH1 \
--motif MA0655.1_JDP2 \
--motif MA0841.1_NFE2 \
--output VariantEffects/${setUp}_${TFgroup}_motifs.txt

# fimo on sequences with TeloHEAC_IL1b_6h top variant effects using TF motifs
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/FIMO.sh VariantEffects/${setUp}_${TFgroup}_motifs.txt VariantEffects/10percent__df_highestVariantEffects_mean_${setUp}.fa

# how many of the top variants effects can I explain with them lying in one of the TF motifs?
python VariantEffects/check_variant_within_motif.py \
--activityFile VariantEffects/df_highestVariantEffects_mean_${setUp}.csv \
--expSetUp ${setUp} \
--output VariantEffects/${setUp}_${TFgroup}




TFgroup="group5" #--motif auch ändern

echo ${TFgroup}
# generate motifs in fimo format
python ./ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/sanityCheck_extract_motifPWM.py \
--motif MA0466.3_CEBPB \
--motif MA0838.1_CEBPG \
--motif MA0837.2_CEBPE \
--output VariantEffects/${setUp}_${TFgroup}_motifs.txt

# fimo on TeloHEAC_IL1b_6h top variant effects
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/FIMO.sh VariantEffects/${setUp}_${TFgroup}_motifs.txt VariantEffects/10percent__df_highestVariantEffects_mean_${setUp}.fa

# how many of the top variants effects can I explain with them lying in one of the TF motifs?
python VariantEffects/check_variant_within_motif.py \
--activityFile VariantEffects/df_highestVariantEffects_mean_${setUp}.csv \
--expSetUp ${setUp} \
--output VariantEffects/${setUp}_${TFgroup}

