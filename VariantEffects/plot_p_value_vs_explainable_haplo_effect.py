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
def cli(infile_explainable, infile_all, p_val_type, out_path):

    df_exp=pd.read_csv(infile_explainable, sep=",", low_memory=False)
    df_exp=df_exp.drop_duplicates(subset=["ID1","ID2"]) #to remove cases where multiple motifs match or multiple variants
    if p_val_type=="P.Value":
        df_exp = df_exp.sort_values (by=["p-value haplotype effect"])
    elif p_val_type=="adj.P.Val":
        df_exp = df_exp.sort_values (by=["adj. p-value haplotype effect"])
    
    

    df_all=pd.read_csv(infile_all, sep=",", low_memory=False)
    df_all=df_all.drop_duplicates(subset=["ID1","ID2"]) 
    df_all=df_all.sort_values(by=[p_val_type])
    p_val_array=df_all[p_val_type].to_list() #sorted already



    number_of_explainable_effects=[]
    for i in range(0, len(p_val_array), 1):
        """"
        if p_val_type=="P.Value":
            df_temp_exp=df_exp.loc[df_exp['p-value haplotype effect'] <= p_val_array[i]]
        elif p_val_type=="adj.P.Val":
            df_temp_exp=df_exp.loc[df_exp['adj. p-value haplotype effect'] <= p_val_array[i]]
        """
        df_temp_all=df_all.loc[df_all[p_val_type] <= p_val_array[i]]
        df_temp_exp=pd.merge(df_temp_all, df_exp, on=["ID1", "ID2"], how="inner")
        number_of_explainable_effects.append(df_temp_exp.shape[0]/df_temp_all.shape[0])

    plt.vlines(x=np.log10(p_val_array)[100::100], ymin=0, ymax=1, colors='purple', ls='--', lw=0.5, label="100 p-val interval")
    plt.scatter(np.log10(p_val_array), number_of_explainable_effects, s=10)
    plt.ylabel("fraction of explainable effects")
    plt.xlabel("log_10(p-value threshold for haplotype effect)")
    plt.legend()

    plt.savefig(str(out_path)+"/exp_effects_pvalue_"+str(p_val_type)+".svg")


if __name__ == "__main__":
    cli()