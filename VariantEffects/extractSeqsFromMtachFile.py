"""
Description:       
                    

Example commands:   python VariantEffects/extractVariantEffectsFromMatchFile.py --quantile 0.1 --output test
Outputs:            fasta sequence file with sequences cooresponding to variants with top x % of all variant effects. CSV with haplo seqs and variant effect

"""

#set seed and load dependencies

import click
from tqdm import tqdm #progress bar
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from scipy import stats
import os

@click.command()
@click.option(
    "--input",
    "infile",
    required=True,
    multiple=False,
    type=str,
    default="MPRAlm_Alec_all_haplo_SNPs_diffTeloHEAC_CTRL_vs_6h.csv",
    help="e.g. MPRAlm_Alec_all_haplo_SNPs_diffTeloHEAC_CTRL_vs_6h.csv",
)
def cli(infile):

    

    #import sequences
    seqs=open(str(infile)+".csv", "r")
    seq_entries=seqs.readlines()
    seqs.close()
    print (len(seq_entries))
    print(len(seq_entries[1].split(",")))
    print(seq_entries[1].split(","))
        
    # write to otfile
    outfile=open(str(infile)+"_seqs"+".fa", "w")
    for i in tqdm(range (1, len(seq_entries), 1)):
        print("hallo", seq_entries[i])
        outfile.write(">"+str(seq_entries[i].split(",")[10]) + "\n" + str(seq_entries[i].split(",")[12]) + "\n" + ">"+str(seq_entries[i].split(",")[11]) + "\n" + str(seq_entries[i].split(",")[25]) + "\n")
    outfile.close()
            
if __name__ == "__main__":
    cli()












