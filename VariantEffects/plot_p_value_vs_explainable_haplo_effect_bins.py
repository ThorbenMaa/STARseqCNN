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
    "--input_explainable",
    "infile_explainable",
    required=True,
    multiple=False,
    type=str,
    default="MPRAlm_Alec_diff_HASMC_ctrl_vs_chol_withMotifs.csv",
    help="e.g. MPRAlm_Alec_diff_HASMC_ctrl_vs_chol_withMotifs.csv",
)
@click.option(
    "--input_all",
    "infile_all",
    required=True,
    multiple=False,
    type=str,
    default="MPRAlm_Alec_all_haplo_SNPs_diffTeloHEAC_CTRL_vs_6h.csv",
    help="e.g. MPRAlm_Alec_all_haplo_SNPs_diffTeloHEAC_CTRL_vs_6h.csv",
)
@click.option(
    "--bin_size",
    "bin_size",
    required=True,
    multiple=False,
    type=int,
    default=100,
    help="p value bin size",
)
@click.option(
    "--p_val_type",
    "p_val_type",
    required=True,
    multiple=False,
    type=str,
    default="P.Value",
    help="P.Value or adj.P.Val",
)
@click.option(
    "--out_path",
    "out_path",
    required=True,
    multiple=False,
    type=str,
    help="path where output file should end up",
)
def cli(infile_explainable, infile_all, bin_size, p_val_type, out_path):

    df_exp=pd.read_csv(infile_explainable, sep=",", low_memory=False)
    df_exp=df_exp.drop_duplicates(subset=["ID1","ID2"]) #to remove cases where multiple motifs match or multiple variants
    #df_exp = df_exp.sort_values (by=["p-value haplotype effect"])
    #print(df_exp)
    
    

    df_all=pd.read_csv(infile_all, sep=",", low_memory=False)
    print(df_all)
    df_all=df_all.drop_duplicates(subset=["ID1","ID2"]) 
    df_all=df_all.sort_values(by=[p_val_type])
    df_all=df_all.reset_index(drop=True)
    #p_val_array=df_all["P.Value"].to_list() #sorted already
    #print (len(p_val_array))
    #print(df_all["ID1"])



    number_of_explainable_effects=[]
    for i in range(bin_size, df_all.shape[0], bin_size):
        df_temp_all=df_all.iloc[i-bin_size:i]
        #print(df_temp_all["P.Value"])
        df_temp_exp=pd.merge(df_temp_all, df_exp, on=["ID1", "ID2"], how="inner")
        #print(df_temp_exp)
        number_of_explainable_effects.append(df_temp_exp.shape[0]/bin_size)

    np_bins=np.arange(1, len(number_of_explainable_effects)+1)
    print(len(number_of_explainable_effects))
    #print(len(np_bins.shape))
    plt.bar(np_bins, number_of_explainable_effects)
    plt.ylabel("fraction of explainable effects")
    plt.xlabel("p-value bin n = "+str(bin_size))
    
    
    #clacl std dev of explainable fraction
    print ("mean")
    print (np.mean(np.array(number_of_explainable_effects)))
    print ("std dev")
    print (np.std(np.array(number_of_explainable_effects)))

    plt.savefig(str(out_path)+"/exp_effects_pvalue_bins_"+str(p_val_type)+".svg")



if __name__ == "__main__":
    cli()