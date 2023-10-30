#!/bin/bash
#Description: the script needs a vep annotated vcf file as input. It calls scripts that add ESM score annotations in the info column of the vcf file. The ESM annotated vcf file is the output of the script.
#Author: Thorben Maass
#Contact: tho.maass@uni-luebeck.de
#Year: 2023
#example command: bash ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/FIMO.sh test.txt test.fa
work_dir=$(pwd)
export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.4:$PATH
fimo $1 $2