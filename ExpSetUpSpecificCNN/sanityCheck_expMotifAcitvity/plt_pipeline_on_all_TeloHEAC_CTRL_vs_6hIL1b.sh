#!/bin/bash
#command: bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline_on_all.sh
source /etc/profile.d/modules.sh

echo group1
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline.sh \
"\
--motif MA1130.1_FOSL2::JUN \
--motif MA1141.1_FOS::JUND \
--motif MA1128.1_FOSL1::JUN \
" \
"\
--setUps TeloHAEC_CTRL \
--setUps TeloHAEC_IL1b_6h \
" \
"--output group1"

echo group2
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline.sh \
"\
--motif MA0107.1_RELA \
--motif MA0101.1_REL \
" \
"\
--setUps TeloHAEC_CTRL \
--setUps TeloHAEC_IL1b_6h \
" \
"--output group2"

echo group3
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline.sh \
"\
--motif MA0488.1_JUN \
--motif MA0018.4_CREB1 \
--motif MA0492.1_JUND \
" \
"\
--setUps TeloHAEC_CTRL \
--setUps TeloHAEC_IL1b_6h \
" \
"--output group3"

echo group4
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline.sh \
"\
--motif MA1633.2_BACH1 \
--motif MA0655.1_JDP2 \
--motif MA0841.1_NFE2 \
" \
"\
--setUps TeloHAEC_CTRL \
--setUps TeloHAEC_IL1b_6h \
" \
"--output group4"

echo group5
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline.sh \
"\
--motif MA0466.3_CEBPB \
--motif MA0838.1_CEBPG \
--motif MA0837.2_CEBPE \
" \
"\
--setUps TeloHAEC_CTRL \
--setUps TeloHAEC_IL1b_6h \
" \
"--output group5"
