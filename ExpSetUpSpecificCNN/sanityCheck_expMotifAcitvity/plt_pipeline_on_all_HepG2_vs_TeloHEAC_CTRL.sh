#!/bin/bash
#command: bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline_on_all.sh
source /etc/profile.d/modules.sh

echo group1
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline.sh \
"\
--motif MA1950.1_FLI1::FOXI1 \
--motif MA1942.1_ETV2::FOXI1 \
--motif MA1956.1_FOXO1::FLI1 \
" \
"\
--setUps HepG2_untreatedPilot \
--setUps TeloHAEC_CTRL \
" \
"--output group1"

echo group2
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline.sh \
"\
--motif MA0655.1_JDP2 \
--motif MA0841.1_NFE2 \
--motif MA1135.1_FOSB::JUNB \
" \
"\
--setUps HepG2_untreatedPilot \
--setUps TeloHAEC_CTRL \
" \
"--output group2"

echo group3a
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline.sh \
"\
--motif MA1653.1_ZNF148 \
--motif MA0516.3_SP2 \
--motif MA0742.2_KLF12 \
--motif MA1513.1_KLF15 \
--motif MA1511.2_KLF10 \
--motif MA0753.2_ZNF740 \
" \
"\
--setUps HepG2_untreatedPilot \
--setUps TeloHAEC_CTRL \
" \
"--output group3a"

echo group3b
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline.sh \
"\
--motif MA1522.1_MAZ \
--motif MA1723.1_PRDM9 \
--motif MA1965.1_SP5 \
--motif MA0493.2_KLF1 \
--motif MA0740.2_KLF14 \
--motif MA0599.1_KLF5 \
" \
"\
--setUps HepG2_untreatedPilot \
--setUps TeloHAEC_CTRL \
" \
"--output group3b"

echo group4
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline.sh \
"\
--motif MA1633.2_BACH1 \
--motif MA0099.3_FOS::JUN \
--motif MA1988.1_Atf3 \
" \
"\
--setUps HepG2_untreatedPilot \
--setUps TeloHAEC_CTRL \
" \
"--output group4"

echo group5
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline.sh \
"\
--motif MA1120.1_SOX13 \
--motif MA0077.1_SOX9 \
--motif MA1152.1_SOX15 \
" \
"\
--setUps HepG2_untreatedPilot \
--setUps TeloHAEC_CTRL \
" \
"--output group5"


echo group6
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline.sh \
"\
--motif MA0679.2_ONECUT1 \
" \
"\
--setUps HepG2_untreatedPilot \
--setUps TeloHAEC_CTRL \
" \
"--output group6"

echo group7
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline.sh \
"\
--motif MA0036.3_GATA2 \
--motif MA0037.4_Gata3 \
" \
"\
--setUps HepG2_untreatedPilot \
--setUps TeloHAEC_CTRL \
" \
"--output group7"

echo group8
bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plt_pipeline.sh \
"\
--motif MA1525.2_NFATC4 \
" \
"\
--setUps HepG2_untreatedPilot \
--setUps TeloHAEC_CTRL \
" \
"--output group8"
