"""
Description:            This scripts plots Boxplots without outlyers of experimental STARseq activity and differences for Sequences containing a motif of interest or not. Note that the motifs should be provided by the user
                        using the "motivsToCheck" array. The array can have as many motifs as wanted.

Input:                  Input 1 are labels (experimental activities). Input 2 are sequences. Note that the motif of interest has to be provided using the "motivsToCheck" array in this file and are not parsed together with the
                        bash command! Input 3 is like input 1

Output:                  Boxplots for activities with and without motifs of interest for TeloHEAC_CTRL and TeloHEAC treated with IL1b after 6h

example bash command:   python ./ExpSetUpSpecificCNN/sanityCheck_modisco_results_diff_TeloHEAC_CTRL_vs_IL1b_6h.py 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo1.txt starrseq-all-final-toorder_oligocomposition.csv 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo2.txt

"""


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os
os.environ["TF_CPP_MIN_LOG_LEVEL"]="2" #to supress warnings with loading tf
import tensorflow as tf
from tensorflow import keras
import numpy as np
import sys
from scipy import stats
import random


motivsToCheck=[ "TGAGTCA", "GGAATTTCCC", "ATGATGTCA", "TGACTCA", "TTGCCCAA"]

#define cell type
cellType=[ "TeloHAEC_CTRL", "TeloHAEC_IL1b_6h"]#"cell_3T3_diff_CTRL", "ccell_3T3_undiff_CTRL", "cell_3T3_undiff_TGFB", "RAW_CTRL" , "RAW_IL1B", "RAW_TGFB", "TeloHAEC_CTRL", "TeloHAEC_IL1b_24h", "TeloHAEC_IL1b_6h", "HASMC_untreatedPilot", "HASMC_Chol", "HepG2_untreatedPilot"]

#import labels
df_IDs_reg_labels=pd.read_csv(sys.argv[1], sep="\t", decimal=',', low_memory=False)
df_IDs_reg_labels2=pd.read_csv(sys.argv[3], sep="\t", decimal=',', low_memory=False)

df_IDs_reg_labels = pd.concat([df_IDs_reg_labels, df_IDs_reg_labels2], axis=0)

df_IDs_reg_labels=df_IDs_reg_labels.drop_duplicates()

#import sequences
df_IDs_Sequences=pd.read_csv(sys.argv[2], sep=",",low_memory=False)

#remove ">" as first character from ID column
df_IDs_Sequences["name"]=df_IDs_Sequences["name"].str[1:]

#drop duplicates
df_IDs_Sequences=df_IDs_Sequences.dropna()

#merge data frames on name/oligo columns
df_IDs_Sequences=df_IDs_Sequences.rename(columns={'name': 'ID'})
df_IDs_reg_labels=df_IDs_reg_labels.rename(columns={'Oligo': 'ID'})
df_IDs_seqs_reg_labels=pd.merge(df_IDs_Sequences, df_IDs_reg_labels, on=["ID"])

#select sequences from data set (here, no filtering is done)
df_IDs_seqs_reg_labels_all=df_IDs_seqs_reg_labels

#average data over replica to generate labels
#3T3
#df_IDs_seqs_reg_labels['mean_cell_3T3_diff_CTRL'] = df_IDs_seqs_reg_labels.loc[:, ["cell_3T3_diff_CTRL_rep1_2022_12_14", "cell_3T3_diff_CTRL_rep2_2022_12_14", "cell_3T3_diff_CTRL_rep3_2022_12_14", "cell_3T3_diff_CTRL_rep4_2022_12_14"]].mean(axis=1)
#df_IDs_seqs_reg_labels['mean_ccell_3T3_undiff_CTRL'] = df_IDs_seqs_reg_labels.loc[:, ["cell_3T3_undiff_CTRL_rep1_2022_12_14", "cell_3T3_undiff_CTRL_rep2_2022_12_14", "cell_3T3_undiff_CTRL_rep3_2022_12_14", "cell_3T3_undiff_CTRL_rep4_2022_12_14"]].mean(axis=1)
#df_IDs_seqs_reg_labels['mean_cell_3T3_undiff_TGFB'] = df_IDs_seqs_reg_labels.loc[:, ["cell_3T3_undiff_TGFB_rep1_2022_12_14", "cell_3T3_undiff_TGFB_rep2_2022_12_14", "cell_3T3_undiff_TGFB_rep3_2022_12_14", "cell_3T3_undiff_TGFB_rep4_2022_12_14"]].mean(axis=1)
#RAW
#df_IDs_seqs_reg_labels['mean_RAW_CTRL'] = df_IDs_seqs_reg_labels.loc[:, ["RAW_CTRL_rep1_2022_12_14", "RAW_CTRL_rep2_2022_12_14", "RAW_CTRL_rep3_2022_12_14", "RAW_CTRL_rep4_2022_12_14"]].mean(axis=1)
#df_IDs_seqs_reg_labels['mean_RAW_IL1B'] = df_IDs_seqs_reg_labels.loc[:, ["RAW_IL1B_rep1_2022_12_14", "RAW_IL1B_rep2_2022_12_14", "RAW_IL1B_rep3_2022_12_14", "RAW_IL1B_rep4_2022_12_14"]].mean(axis=1)
#df_IDs_seqs_reg_labels['mean_RAW_TGFB'] = df_IDs_seqs_reg_labels.loc[:, ["RAW_TGFB_rep1_2022_12_14", "RAW_TGFB_rep2_2022_12_14", "RAW_TGFB_rep3_2022_12_14", "RAW_TGFB_rep4_2022_12_14"]].mean(axis=1)
#TeloHAEC
df_IDs_seqs_reg_labels['mean_TeloHAEC_CTRL'] = df_IDs_seqs_reg_labels.loc[:, ["TeloHAEC_CTRL_rep1", "TeloHAEC_CTRL_rep2", "TeloHAEC_CTRL_rep3"]].mean(axis=1)
#df_IDs_seqs_reg_labels['mean_TeloHAEC_IL1b_24h'] = df_IDs_seqs_reg_labels.loc[:, ["TeloHAEC_IL1b_24h_rep1", "TeloHAEC_IL1b_24h_rep2", "TeloHAEC_IL1b_24h_rep3"]].mean(axis=1)
df_IDs_seqs_reg_labels['mean_TeloHAEC_IL1b_6h'] = df_IDs_seqs_reg_labels.loc[:, ["TeloHAEC_IL1b_6h_rep1", "TeloHAEC_IL1b_6h_rep2", "TeloHAEC_IL1b_6h_rep3"]].mean(axis=1)
#HASMC
#df_IDs_seqs_reg_labels['mean_HASMC_untreatedPilot'] = df_IDs_seqs_reg_labels.loc[:, ["HASMC_untreatedPilot_rep1", "HASMC_untreatedPilot_rep2", "HASMC_untreatedPilot_rep3"]].mean(axis=1)
#df_IDs_seqs_reg_labels['mean_HASMC_Chol'] = df_IDs_seqs_reg_labels.loc[:, ["HASMC_Chol_rep1", "HASMC_Chol_rep2", "HASMC_Chol_rep3"]].mean(axis=1)
#HepG2
#df_IDs_seqs_reg_labels['mean_HepG2_untreatedPilot'] = df_IDs_seqs_reg_labels.loc[:, ["HepG2_untreatedPilot_rep1", "HepG2_untreatedPilot_rep2", "HepG2_untreatedPilot_rep3"]].mean(axis=1)

#split data to dataframe corresponding to sequences having or not having the motif of interest.
for i in range (0, len(motivsToCheck), 1):
    data = []
    xTicks= []
    for j in range (0, len(cellType), 1):
    
    
        df_IDs_seqs_reg_labels_hasMotiv=df_IDs_seqs_reg_labels.loc[(df_IDs_seqs_reg_labels['enhancer'].str.contains(str(motivsToCheck[i])))==True]["mean_"+str(cellType[j])].to_numpy()
        data.append(df_IDs_seqs_reg_labels_hasMotiv)
        xTicks.append(str(cellType[j])+" +")
        df_IDs_seqs_reg_labels_hasntMotiv=df_IDs_seqs_reg_labels.loc[(df_IDs_seqs_reg_labels['enhancer'].str.contains(str(motivsToCheck[i])))==False]["mean_"+str(cellType[j])].to_numpy()
        data.append(df_IDs_seqs_reg_labels_hasntMotiv)
        xTicks.append(str(cellType[j])+" -")
    
    # differences
    df_IDs_seqs_reg_labels_hasMotiv_diff = (df_IDs_seqs_reg_labels.loc[(df_IDs_seqs_reg_labels['enhancer'].str.contains(str(motivsToCheck[i])))==True]["mean_"+str(cellType[0])].to_numpy()
                                            -df_IDs_seqs_reg_labels.loc[(df_IDs_seqs_reg_labels['enhancer'].str.contains(str(motivsToCheck[i])))==True]["mean_"+str(cellType[1])].to_numpy()
                                            )
    data.append(df_IDs_seqs_reg_labels_hasMotiv_diff)
    xTicks.append("Difference +")

    df_IDs_seqs_reg_labels_hasNotMotiv_diff = (df_IDs_seqs_reg_labels.loc[(df_IDs_seqs_reg_labels['enhancer'].str.contains(str(motivsToCheck[i])))==False]["mean_"+str(cellType[0])].to_numpy()
                                            -df_IDs_seqs_reg_labels.loc[(df_IDs_seqs_reg_labels['enhancer'].str.contains(str(motivsToCheck[i])))==False]["mean_"+str(cellType[1])].to_numpy()
                                            )
    data.append(df_IDs_seqs_reg_labels_hasNotMotiv_diff)
    xTicks.append("Difference -")

        
        
    print(data)
    box=plt.boxplot(data, showfliers=False, patch_artist=True)
    colors = ['cyan', 'r', "cyan", "r", "cyan", "r"]
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
    #plt.set_ylim=(-2, 100)
    plt.ylabel("experimental activity")
    plt.xticks([1,2,3,4,5,6], xTicks, rotation=90)
    plt.title(" n("+str(motivsToCheck[i])+")="+str(data[0].shape[0])+"; "+'n(other)='+str(data[1].shape[0]))    
    plt.tight_layout()       
    plt.savefig("boxplot_diff_"+str(motivsToCheck[i])+".svg")
    plt.close()

