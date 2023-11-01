"""
Description:       
                    

Example commands:   python VariantEffects/extractSeqsWithHighVariantEffect.py --quantile 0.1 --output test
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

@click.command()
@click.option(
    "--seqFile",
    "seq_file",
    required=True,
    multiple=False,
    type=str,
    default="starrseq-all-final-toorder_oligocomposition.csv",
    help="e.g. starrseq-all-final-toorder_oligocomposition.csv",
)
@click.option(
    "--activityFile",
    "activity_file",
    required=True,
    multiple=True,
    type=str,
    default=["2023-01-10_22-29-33 myCounts.minDNAfilt.depthNorm.keepHaps - starr.haplotypes.oligo1.txt", "2023-01-10_22-29-33 myCounts.minDNAfilt.depthNorm.keepHaps - starr.haplotypes.oligo2.txt" ],
    help="e.g. 2023-01-10_22-29-33 myCounts.minDNAfilt.depthNorm.keepHaps - starr.haplotypes.oligo1.txt",
)
@click.option(
    "--quantile",
    "quantile",
    required=True,
    multiple=False,
    type=float,
    default=0.1,
    help="bla",
)
@click.option(
    "--output",
    "out",
    required=True,
    multiple=False,
    type=str,
    default="bla",
    help="bla",
)
def cli(seq_file, activity_file, quantile, out):

    #import labels
    df_IDs_reg_labels=pd.read_csv(activity_file[0], sep="\t", decimal=',', low_memory=False)
    df_IDs_reg_labels2=pd.read_csv(activity_file[1], sep="\t", decimal=',', low_memory=False)

    df_IDs_reg_labels = pd.concat([df_IDs_reg_labels, df_IDs_reg_labels2], axis=0)

    df_IDs_reg_labels=df_IDs_reg_labels.drop_duplicates()


    #import sequences
    df_IDs_Sequences=pd.read_csv(seq_file, sep=",",low_memory=False)

    #remove ">" as first character from ID column
    df_IDs_Sequences["name"]=df_IDs_Sequences["name"].str[1:]

    #drop duplicates
    df_IDs_Sequences=df_IDs_Sequences.dropna()

    #merge data frames on name/oligo columns
    df_IDs_Sequences=df_IDs_Sequences.rename(columns={'name': 'ID'})
    df_IDs_reg_labels=df_IDs_reg_labels.rename(columns={'Oligo': 'ID'})
    df_IDs_seqs_reg_labels=pd.merge(df_IDs_Sequences, df_IDs_reg_labels, on=["ID"])

    #average data over replica to generate labels, systematic +1 to be able to log transform, normalize by dividing by mean_input2022Dec

    #input for normalization
    df_IDs_seqs_reg_labels['mean_input2022Dec'] = df_IDs_seqs_reg_labels.loc[:, ["input2022Dec_50ng_rep1_2022_12_14", "input2022Dec_50ng_rep2_2022_12_14", "input2022Dec_50ng_rep3_2022_12_14", "input2022Dec_50ng_rep4_2022_12_14"]].mean(axis=1) + 1

    #3T3
    df_IDs_seqs_reg_labels['mean_cell_3T3_diff_CTRL'] = df_IDs_seqs_reg_labels.loc[:, ["cell_3T3_diff_CTRL_rep1_2022_12_14", "cell_3T3_diff_CTRL_rep2_2022_12_14", "cell_3T3_diff_CTRL_rep3_2022_12_14", "cell_3T3_diff_CTRL_rep4_2022_12_14"]].mean(axis=1) + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    df_IDs_seqs_reg_labels['mean_ccell_3T3_undiff_CTRL'] = df_IDs_seqs_reg_labels.loc[:, ["cell_3T3_undiff_CTRL_rep1_2022_12_14", "cell_3T3_undiff_CTRL_rep2_2022_12_14", "cell_3T3_undiff_CTRL_rep3_2022_12_14", "cell_3T3_undiff_CTRL_rep4_2022_12_14"]].mean(axis=1) + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    df_IDs_seqs_reg_labels['mean_cell_3T3_undiff_TGFB'] = df_IDs_seqs_reg_labels.loc[:, ["cell_3T3_undiff_TGFB_rep1_2022_12_14", "cell_3T3_undiff_TGFB_rep2_2022_12_14", "cell_3T3_undiff_TGFB_rep3_2022_12_14", "cell_3T3_undiff_TGFB_rep4_2022_12_14"]].mean(axis=1) + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    #RAW
    df_IDs_seqs_reg_labels['mean_RAW_CTRL'] = df_IDs_seqs_reg_labels.loc[:, ["RAW_CTRL_rep1_2022_12_14", "RAW_CTRL_rep2_2022_12_14", "RAW_CTRL_rep3_2022_12_14", "RAW_CTRL_rep4_2022_12_14"]].mean(axis=1) + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    df_IDs_seqs_reg_labels['mean_RAW_IL1B'] = df_IDs_seqs_reg_labels.loc[:, ["RAW_IL1B_rep1_2022_12_14", "RAW_IL1B_rep2_2022_12_14", "RAW_IL1B_rep3_2022_12_14", "RAW_IL1B_rep4_2022_12_14"]].mean(axis=1) + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    df_IDs_seqs_reg_labels['mean_RAW_TGFB'] = df_IDs_seqs_reg_labels.loc[:, ["RAW_TGFB_rep1_2022_12_14", "RAW_TGFB_rep2_2022_12_14", "RAW_TGFB_rep3_2022_12_14", "RAW_TGFB_rep4_2022_12_14"]].mean(axis=1) + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    #TeloHAEC
    df_IDs_seqs_reg_labels['mean_TeloHAEC_CTRL'] = df_IDs_seqs_reg_labels.loc[:, ["TeloHAEC_CTRL_rep1", "TeloHAEC_CTRL_rep2", "TeloHAEC_CTRL_rep3"]].mean(axis=1) + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    df_IDs_seqs_reg_labels['mean_TeloHAEC_IL1b_24h'] = df_IDs_seqs_reg_labels.loc[:, ["TeloHAEC_IL1b_24h_rep1", "TeloHAEC_IL1b_24h_rep2", "TeloHAEC_IL1b_24h_rep3"]].mean(axis=1) + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    df_IDs_seqs_reg_labels['mean_TeloHAEC_IL1b_6h'] = df_IDs_seqs_reg_labels.loc[:, ["TeloHAEC_IL1b_6h_rep1", "TeloHAEC_IL1b_6h_rep2", "TeloHAEC_IL1b_6h_rep3"]].mean(axis=1) + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    #HASMC
    df_IDs_seqs_reg_labels['mean_HASMC_untreatedPilot'] = df_IDs_seqs_reg_labels.loc[:, ["HASMC_untreatedPilot_rep1", "HASMC_untreatedPilot_rep2", "HASMC_untreatedPilot_rep3"]].mean(axis=1) + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    df_IDs_seqs_reg_labels['mean_HASMC_Chol'] = df_IDs_seqs_reg_labels.loc[:, ["HASMC_Chol_rep1", "HASMC_Chol_rep2", "HASMC_Chol_rep3"]].mean(axis=1) + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    #HepG2
    df_IDs_seqs_reg_labels['mean_HepG2_untreatedPilot'] = df_IDs_seqs_reg_labels.loc[:, ["HepG2_untreatedPilot_rep1", "HepG2_untreatedPilot_rep2", "HepG2_untreatedPilot_rep3"]].mean(axis=1) + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']

    #print (df_IDs_seqs_reg_labels)
    df_IDs_seqs_reg_labels=df_IDs_seqs_reg_labels.sort_values(by=["ID"])
    #print (df_IDs_seqs_reg_labels)


    #initalize output df
    df_diff=pd.DataFrame(columns=["Seq1",
                                "ID1",
                                "mean_cell_3T3_diff_CTRL1",
                                "mean_ccell_3T3_undiff_CTRL1",
                                "mean_cell_3T3_undiff_TGFB1",
                                "mean_RAW_CTRL1",
                                "mean_RAW_IL1B1",
                                "mean_RAW_TGFB1",
                                "mean_TeloHAEC_CTRL1",
                                "mean_TeloHAEC_IL1b_24h1",
                                "mean_TeloHAEC_IL1b_6h1",
                                "mean_HASMC_untreatedPilot1",
                                "mean_HASMC_Chol1",
                                "mean_HepG2_untreatedPilot1",
                                "Seq2",
                                "ID2",
                                "mean_cell_3T3_diff_CTRL2",
                                "mean_ccell_3T3_undiff_CTRL2",
                                "mean_cell_3T3_undiff_TGFB2",
                                "mean_RAW_CTRL2",
                                "mean_RAW_IL1B2",
                                "mean_RAW_TGFB2",
                                "mean_TeloHAEC_CTRL2",
                                "mean_TeloHAEC_IL1b_24h2",
                                "mean_TeloHAEC_IL1b_6h2",
                                "mean_HASMC_untreatedPilot2",
                                "mean_HASMC_Chol2",
                                "mean_HepG2_untreatedPilot2",
                                "variantPos"
                                ])

    #idetify haplotype pairs
    df_IDs_seqs_reg_labels_list=df_IDs_seqs_reg_labels.values.tolist()
    #print(df_IDs_seqs_reg_labels.iloc[[1]])
    #print (str(df_IDs_seqs_reg_labels_list[1][0]))#)[:-1])

    for i in tqdm(range (0, len(df_IDs_seqs_reg_labels_list), 1)):
        for j in range (0, len(df_IDs_seqs_reg_labels_list), 1):
            if str(df_IDs_seqs_reg_labels_list[i][0])[:-1] == str(df_IDs_seqs_reg_labels_list[j][0])[:-1]:
                if str(df_IDs_seqs_reg_labels_list[i][0]) != str(df_IDs_seqs_reg_labels_list[j][0]): #also "gleiche" seq aber anderer haplotype
                    
                    # calulcate variantPos
                    variantPos=False
                    seq1=df_IDs_seqs_reg_labels.iloc[[i]]["enhancer"].values[0]
                    seq2=df_IDs_seqs_reg_labels.iloc[[j]]["enhancer"].values[0]
                    for k in range (0, len(seq1), 1):
                        if seq1[k]!=seq2[k]:
                            variantPos=k+1
                            break

                    # fill df
                    df_diff.loc[len(df_diff)]=[df_IDs_seqs_reg_labels.iloc[[i]]["enhancer"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[i]]["ID"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[i]]["mean_cell_3T3_diff_CTRL"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[i]]["mean_ccell_3T3_undiff_CTRL"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[i]]["mean_cell_3T3_undiff_TGFB"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[i]]["mean_RAW_CTRL"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[i]]["mean_RAW_IL1B"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[i]]["mean_RAW_TGFB"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[i]]["mean_TeloHAEC_CTRL"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[i]]["mean_TeloHAEC_IL1b_24h"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[i]]["mean_TeloHAEC_IL1b_6h"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[i]]["mean_HASMC_untreatedPilot"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[i]]["mean_HASMC_Chol"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[i]]["mean_HepG2_untreatedPilot"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[j]]["enhancer"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[j]]["ID"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[j]]["mean_cell_3T3_diff_CTRL"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[j]]["mean_ccell_3T3_undiff_CTRL"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[j]]["mean_cell_3T3_undiff_TGFB"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[j]]["mean_RAW_CTRL"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[j]]["mean_RAW_IL1B"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[j]]["mean_RAW_TGFB"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[j]]["mean_TeloHAEC_CTRL"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[j]]["mean_TeloHAEC_IL1b_24h"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[j]]["mean_TeloHAEC_IL1b_6h"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[j]]["mean_HASMC_untreatedPilot"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[j]]["mean_HASMC_Chol"].values[0],
                            df_IDs_seqs_reg_labels.iloc[[j]]["mean_HepG2_untreatedPilot"].values[0],
                            variantPos,


                            ]#,
                            #columns= df_diff.columns), ignore_index=True)

                    

                                
                    
                    #df_diff.loc[len(df_diff)]=df_list
        """
        if i==2:
            break
        """
    #calc diff of activities 
    df_diff["diff_activity mean_cell_3T3_diff_CTRL"]=(df_diff["mean_cell_3T3_diff_CTRL1"]) - (df_diff["mean_cell_3T3_diff_CTRL2"])
    df_diff["diff_activity mean_ccell_3T3_undiff_CTRL"]=(df_diff["mean_ccell_3T3_undiff_CTRL1"]) - (df_diff["mean_ccell_3T3_undiff_CTRL2"])
    df_diff["diff_activity mean_cell_3T3_undiff_TGFB"]=(df_diff["mean_cell_3T3_undiff_TGFB1"]) - (df_diff["mean_cell_3T3_undiff_TGFB2"])
    df_diff["diff_activity mean_RAW_CTRL"]=(df_diff["mean_RAW_CTRL1"]) - (df_diff["mean_RAW_CTRL2"])
    df_diff["diff_activity mean_RAW_IL1B"]=(df_diff["mean_RAW_IL1B1"]) - (df_diff["mean_RAW_IL1B2"])
    df_diff["diff_activity mean_RAW_TGFB"]=(df_diff["mean_RAW_TGFB1"]) - (df_diff["mean_RAW_TGFB2"])
    df_diff["diff_activity mean_TeloHAEC_CTRL"]=(df_diff["mean_TeloHAEC_CTRL1"]) - (df_diff["mean_TeloHAEC_CTRL2"])
    df_diff["diff_activity mean_TeloHAEC_IL1b_24h"]=(df_diff["mean_TeloHAEC_IL1b_24h1"]) - (df_diff["mean_TeloHAEC_IL1b_24h2"])
    df_diff["diff_activity mean_TeloHAEC_IL1b_6h"]=(df_diff["mean_TeloHAEC_IL1b_6h1"]) - (df_diff["mean_TeloHAEC_IL1b_6h2"])
    df_diff["diff_activity mean_HASMC_untreatedPilot"]=(df_diff["mean_HASMC_untreatedPilot1"]) - (df_diff["mean_HASMC_untreatedPilot2"])
    df_diff["diff_activity mean_HASMC_Chol"]=(df_diff["mean_HASMC_Chol1"]) - (df_diff["mean_HASMC_Chol2"])
    df_diff["diff_activity mean_HepG2_untreatedPilot"]=(df_diff["mean_HepG2_untreatedPilot1"]) - (df_diff["mean_HepG2_untreatedPilot2"])

    #use abs value (important to be able to drop duplicates)
    df_diff["diff_activity mean_cell_3T3_diff_CTRL"]=df_diff["diff_activity mean_cell_3T3_diff_CTRL"].abs()
    df_diff["diff_activity mean_ccell_3T3_undiff_CTRL"]=df_diff["diff_activity mean_ccell_3T3_undiff_CTRL"].abs()
    df_diff["diff_activity mean_cell_3T3_undiff_TGFB"]=df_diff["diff_activity mean_cell_3T3_undiff_TGFB"].abs()
    df_diff["diff_activity mean_RAW_CTRL"]=df_diff["diff_activity mean_RAW_CTRL"].abs()
    df_diff["diff_activity mean_RAW_IL1B"]=df_diff["diff_activity mean_RAW_IL1B"].abs()
    df_diff["diff_activity mean_RAW_TGFB"]=df_diff["diff_activity mean_RAW_TGFB"].abs()
    df_diff["diff_activity mean_TeloHAEC_CTRL"]=df_diff["diff_activity mean_TeloHAEC_CTRL"].abs()
    df_diff["diff_activity mean_TeloHAEC_IL1b_24h"]=df_diff["diff_activity mean_TeloHAEC_IL1b_24h"].abs()
    df_diff["diff_activity mean_TeloHAEC_IL1b_6h"]=df_diff["diff_activity mean_TeloHAEC_IL1b_6h"].abs()
    df_diff["diff_activity mean_HASMC_untreatedPilot"]=df_diff["diff_activity mean_HASMC_untreatedPilot"].abs()
    df_diff["diff_activity mean_HASMC_Chol"]=df_diff["diff_activity mean_HASMC_Chol"].abs()
    df_diff["diff_activity mean_HepG2_untreatedPilot"]=df_diff["diff_activity mean_HepG2_untreatedPilot"].abs()


    print("before and afte drop nan")
    print (df_diff.shape[0])
    df_diff=df_diff.dropna()
    print (df_diff.shape[0])


    # drop duplicates (weil sonst werte dopplet vorkommen, da a-b == b-a bei abs values). Hier koennte man auch alle machen bei subset, aber wahrscheinlichkeit bei 5 zufaellig gleich zu sein ist schon sehr klein
    print("before and after drop duplicates")
    print (df_diff)
    df_diff=df_diff.drop_duplicates(subset=["diff_activity mean_cell_3T3_diff_CTRL", "diff_activity mean_ccell_3T3_undiff_CTRL", "diff_activity mean_cell_3T3_undiff_TGFB", "diff_activity mean_RAW_CTRL", "diff_activity mean_RAW_IL1B"])
    print (df_diff)

    # select rows with hightest variant effects and write into set up specific data frame
    totalNumberforQuantile=int(quantile*df_diff.shape[0])
    df_highestVariantEffects_mean_cell_3T3_diff_CTRL=df_diff.sort_values('diff_activity mean_cell_3T3_diff_CTRL',ascending=False).iloc[:totalNumberforQuantile,]#["ID1", "Seq1", "ID2", "Seq2", "diff_activity mean_cell_3T3_diff_CTRL", "mean_cell_3T3_diff_CTRL1", "mean_cell_3T3_diff_CTRL2" , "variantPos"]
    df_highestVariantEffects_mean_ccell_3T3_undiff_CTRL=df_diff.sort_values('diff_activity mean_ccell_3T3_undiff_CTRL',ascending=False).iloc[:totalNumberforQuantile,]#["ID1", "Seq1", "ID2", "Seq2", "diff_activity mean_ccell_3T3_undiff_CTRL", "mean_ccell_3T3_undiff_CTRL1", "mean_ccell_3T3_undiff_CTRL2", "variantPos"]
    df_highestVariantEffects_mean_cell_3T3_undiff_TGFB=df_diff.sort_values('diff_activity mean_cell_3T3_undiff_TGFB',ascending=False).iloc[:totalNumberforQuantile,]#["ID1", "Seq1", "ID2", "Seq2", "diff_activity mean_cell_3T3_undiff_TGFB", "mean_cell_3T3_undiff_TGFB1", "mean_cell_3T3_undiff_TGFB2", "variantPos"]
    df_highestVariantEffects_mean_RAW_CTRL=df_diff.sort_values('diff_activity mean_RAW_CTRL',ascending=False).iloc[:totalNumberforQuantile,]#["ID1", "Seq1", "ID2", "Seq2", "diff_activity mean_RAW_CTRL", "mean_RAW_CTRL1", "mean_RAW_CTRL2", "variantPos"]
    df_highestVariantEffects_mean_RAW_IL1B=df_diff.sort_values('diff_activity mean_RAW_IL1B',ascending=False).iloc[:totalNumberforQuantile,]#["ID1", "Seq1", "ID2", "Seq2", "diff_activity mean_RAW_IL1B", "mean_RAW_IL1B1", "mean_RAW_IL1B2", "variantPos"]
    df_highestVariantEffects_mean_RAW_TGFB=df_diff.sort_values('diff_activity mean_RAW_TGFB',ascending=False).iloc[:totalNumberforQuantile,]#["ID1", "Seq1", "ID2", "Seq2", "diff_activity mean_RAW_TGFB", "mean_RAW_TGFB1", "mean_RAW_TGFB2", "variantPos"]
    df_highestVariantEffects_mean_TeloHAEC_CTRL=df_diff.sort_values('diff_activity mean_TeloHAEC_CTRL',ascending=False).iloc[:totalNumberforQuantile,]#["ID1", "Seq1", "ID2", "Seq2", "diff_activity mean_TeloHAEC_CTRL", "mean_TeloHAEC_CTRL1", "mean_TeloHAEC_CTRL2", "variantPos"]
    df_highestVariantEffects_mean_TeloHAEC_IL1b_24h=df_diff.sort_values('diff_activity mean_TeloHAEC_IL1b_24h',ascending=False).iloc[:totalNumberforQuantile,]#["ID1", "Seq1", "ID2", "Seq2", "diff_activity mean_TeloHAEC_IL1b_24h", "mean_TeloHAEC_IL1b_24h1", "mean_TeloHAEC_IL1b_24h2", "variantPos"]
    df_highestVariantEffects_mean_TeloHAEC_IL1b_6h=df_diff.sort_values('diff_activity mean_TeloHAEC_IL1b_6h',ascending=False).iloc[:totalNumberforQuantile,]#["ID1", "Seq1", "ID2", "Seq2", "diff_activity mean_TeloHAEC_IL1b_6h", "mean_TeloHAEC_IL1b_6h1", "mean_TeloHAEC_IL1b_6h2", "variantPos"]
    df_highestVariantEffects_mean_HASMC_untreatedPilot=df_diff.sort_values('diff_activity mean_HASMC_untreatedPilot',ascending=False).iloc[:totalNumberforQuantile,]#["ID1", "Seq1", "ID2", "Seq2", "diff_activity mean_HASMC_untreatedPilot", "mean_HASMC_untreatedPilot1", "mean_HASMC_untreatedPilot2", "variantPos"]
    df_highestVariantEffects_mean_HASMC_Chol=df_diff.sort_values('diff_activity mean_HASMC_Chol',ascending=False).iloc[:totalNumberforQuantile,]#["ID1", "Seq1", "ID2", "Seq2", "diff_activity mean_HASMC_Chol", "mean_HASMC_Chol1", "mean_HASMC_Chol2", "variantPos"]
    df_highestVariantEffects_mean_HepG2_untreatedPilot=df_diff.sort_values('diff_activity mean_HepG2_untreatedPilot',ascending=False).iloc[:totalNumberforQuantile,]#["ID1", "Seq1", "ID2", "Seq2", "diff_activity mean_HepG2_untreatedPilot", "mean_HepG2_untreatedPilot1", "mean_HepG2_untreatedPilot2", "variantPos"]

    df_list=[df_highestVariantEffects_mean_cell_3T3_diff_CTRL,
             df_highestVariantEffects_mean_ccell_3T3_undiff_CTRL,
             df_highestVariantEffects_mean_cell_3T3_undiff_TGFB,
             df_highestVariantEffects_mean_RAW_CTRL,
             df_highestVariantEffects_mean_RAW_IL1B,
             df_highestVariantEffects_mean_RAW_TGFB,
             df_highestVariantEffects_mean_TeloHAEC_CTRL,
             df_highestVariantEffects_mean_TeloHAEC_IL1b_24h,
             df_highestVariantEffects_mean_TeloHAEC_IL1b_6h,
             df_highestVariantEffects_mean_HASMC_untreatedPilot,
             df_highestVariantEffects_mean_HASMC_Chol,
             df_highestVariantEffects_mean_HepG2_untreatedPilot,
             ]
    
    df_list_names=["df_highestVariantEffects_mean_cell_3T3_diff_CTRL",
             "df_highestVariantEffects_mean_ccell_3T3_undiff_CTRL",
             "df_highestVariantEffects_mean_cell_3T3_undiff_TGFB",
             "df_highestVariantEffects_mean_RAW_CTRL",
             "df_highestVariantEffects_mean_RAW_IL1B",
             "df_highestVariantEffects_mean_RAW_TGFB",
             "df_highestVariantEffects_mean_TeloHAEC_CTRL",
             "df_highestVariantEffects_mean_TeloHAEC_IL1b_24h",
             "df_highestVariantEffects_mean_TeloHAEC_IL1b_6h",
             "df_highestVariantEffects_mean_HASMC_untreatedPilot",
             "df_highestVariantEffects_mean_HASMC_Chol",
             "df_highestVariantEffects_mean_HepG2_untreatedPilot",
             ]    

    for dfIndex in range (0, len(df_list), 1):
        
        # save dataframe
        df_list[dfIndex].to_csv("VariantEffects/"+str(df_list_names[dfIndex])+".csv")

        #import sequences
        seqs=open("VariantEffects/"+str(df_list_names[dfIndex])+".csv", "r")
        seq_entries=seqs.readlines()
        seqs.close()
        
        # write to otfile
        outfile=open("VariantEffects/"+str(out)+"_"+str(df_list_names[dfIndex])+".fa", "w")
        for i in tqdm(range (1, len(seq_entries), 1)):
            outfile.write(">"+str(seq_entries[i].split(",")[2]) + "\n" + str(seq_entries[i].split(",")[1]) + "\n" + ">"+str(seq_entries[i].split(",")[16]) + "\n" + str(seq_entries[i].split(",")[15]) + "\n")
        outfile.close()


            
if __name__ == "__main__":
    cli()












