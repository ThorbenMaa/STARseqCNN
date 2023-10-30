"""
Description:            This scripts plots Boxplots without outlyers of experimental STARseq activity and differences for Sequences containing a motif of interest or not. Note that the motifs should be provided by the user
                        using the "motivsToCheck" array. The array can have as many motifs as wanted.


Output:                  Boxplots for activities with and without motifs of interest for TeloHEAC_CTRL and TeloHEAC treated with IL1b after 6h

example bash command:   python ./ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/plot_experimental_activities.py --setUps TeloHAEC_CTRL --setUps TeloHAEC_IL1b_6h
"""


from email.policy import default
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from scipy import stats
import click
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
    "--fimoFile",
    "fimo_file",
    required=True,
    multiple=False,
    type=str,
    default="fimo_out/fimo.tsv",
    help="file generated by FIMO",
)
@click.option(
    "--output",
    "out",
    required=True,
    multiple=False,
    type=str,
    default="diffBoxplot",
    help="e.g. diffBoxplot. How output files will start",
)
@click.option(
    "--setUps",
    "exp_set_ups",
    required=True,
    multiple=True,
    type=str,
    default= [ "TeloHAEC_CTRL", "TeloHAEC_IL1b_6h"],
    help="has to be two",
)
def cli(seq_file, fimo_file, activity_file, out, exp_set_ups):

    #import labels
    df_IDs_reg_labels=pd.DataFrame()
    for i in range(0, len(activity_file), 1):
        df_IDs_reg_labels=pd.concat([df_IDs_reg_labels, pd.read_csv(activity_file[i], sep="\t", decimal=',', low_memory=False)], axis=0, )

    #drop duplicates
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

    #average data over replica to generate labels
    #3T3
    df_IDs_seqs_reg_labels['mean_cell_3T3_diff_CTRL'] = df_IDs_seqs_reg_labels.loc[:, ["cell_3T3_diff_CTRL_rep1_2022_12_14", "cell_3T3_diff_CTRL_rep2_2022_12_14", "cell_3T3_diff_CTRL_rep3_2022_12_14", "cell_3T3_diff_CTRL_rep4_2022_12_14"]].mean(axis=1)
    df_IDs_seqs_reg_labels['mean_ccell_3T3_undiff_CTRL'] = df_IDs_seqs_reg_labels.loc[:, ["cell_3T3_undiff_CTRL_rep1_2022_12_14", "cell_3T3_undiff_CTRL_rep2_2022_12_14", "cell_3T3_undiff_CTRL_rep3_2022_12_14", "cell_3T3_undiff_CTRL_rep4_2022_12_14"]].mean(axis=1)
    df_IDs_seqs_reg_labels['mean_cell_3T3_undiff_TGFB'] = df_IDs_seqs_reg_labels.loc[:, ["cell_3T3_undiff_TGFB_rep1_2022_12_14", "cell_3T3_undiff_TGFB_rep2_2022_12_14", "cell_3T3_undiff_TGFB_rep3_2022_12_14", "cell_3T3_undiff_TGFB_rep4_2022_12_14"]].mean(axis=1)
    #RAW
    df_IDs_seqs_reg_labels['mean_RAW_CTRL'] = df_IDs_seqs_reg_labels.loc[:, ["RAW_CTRL_rep1_2022_12_14", "RAW_CTRL_rep2_2022_12_14", "RAW_CTRL_rep3_2022_12_14", "RAW_CTRL_rep4_2022_12_14"]].mean(axis=1)
    df_IDs_seqs_reg_labels['mean_RAW_IL1B'] = df_IDs_seqs_reg_labels.loc[:, ["RAW_IL1B_rep1_2022_12_14", "RAW_IL1B_rep2_2022_12_14", "RAW_IL1B_rep3_2022_12_14", "RAW_IL1B_rep4_2022_12_14"]].mean(axis=1)
    df_IDs_seqs_reg_labels['mean_RAW_TGFB'] = df_IDs_seqs_reg_labels.loc[:, ["RAW_TGFB_rep1_2022_12_14", "RAW_TGFB_rep2_2022_12_14", "RAW_TGFB_rep3_2022_12_14", "RAW_TGFB_rep4_2022_12_14"]].mean(axis=1)
    #TeloHAEC
    df_IDs_seqs_reg_labels['mean_TeloHAEC_CTRL'] = df_IDs_seqs_reg_labels.loc[:, ["TeloHAEC_CTRL_rep1", "TeloHAEC_CTRL_rep2", "TeloHAEC_CTRL_rep3"]].mean(axis=1)
    df_IDs_seqs_reg_labels['mean_TeloHAEC_IL1b_24h'] = df_IDs_seqs_reg_labels.loc[:, ["TeloHAEC_IL1b_24h_rep1", "TeloHAEC_IL1b_24h_rep2", "TeloHAEC_IL1b_24h_rep3"]].mean(axis=1)
    df_IDs_seqs_reg_labels['mean_TeloHAEC_IL1b_6h'] = df_IDs_seqs_reg_labels.loc[:, ["TeloHAEC_IL1b_6h_rep1", "TeloHAEC_IL1b_6h_rep2", "TeloHAEC_IL1b_6h_rep3"]].mean(axis=1)
    #HASMC
    df_IDs_seqs_reg_labels['mean_HASMC_untreatedPilot'] = df_IDs_seqs_reg_labels.loc[:, ["HASMC_untreatedPilot_rep1", "HASMC_untreatedPilot_rep2", "HASMC_untreatedPilot_rep3"]].mean(axis=1)
    df_IDs_seqs_reg_labels['mean_HASMC_Chol'] = df_IDs_seqs_reg_labels.loc[:, ["HASMC_Chol_rep1", "HASMC_Chol_rep2", "HASMC_Chol_rep3"]].mean(axis=1)
    #HepG2
    df_IDs_seqs_reg_labels['mean_HepG2_untreatedPilot'] = df_IDs_seqs_reg_labels.loc[:, ["HepG2_untreatedPilot_rep1", "HepG2_untreatedPilot_rep2", "HepG2_untreatedPilot_rep3"]].mean(axis=1)

    # import FIMO file with sequences containing motifs of interest
    df_fimo=pd.read_csv(fimo_file, sep="\t", skipfooter=3)
    #print(df_fimo)
    df_fimo=df_fimo.rename(columns={'sequence_name': 'ID'})
    #print(df_fimo)

    
    # generate list with motifs scrutinized in FIMO file
    motivsToCheck=[]
    df_fimo_temp = df_fimo.sort_values(by=["motif_id"])
    df_fimo_temp = df_fimo_temp.drop_duplicates(subset=["motif_id"])
    motivsToCheck=df_fimo_temp["motif_id"].to_list()

    print(motivsToCheck)
    


    #split data to dataframe corresponding to sequences having or not having the motif of interest and extract corresponding experimental activities
    
    fig, axs = plt.subplots(nrows=len(motivsToCheck), ncols=1, sharex=True, sharey=True, figsize=(8,4*len(motivsToCheck)))
    if len(motivsToCheck)==1: # bug repair, da sonst bei nur einem motif ax[i] mit i=0 nicht gefunden wird. Dann ist axs keine liste sondern nur ein ding
        axs=[axs]
    fig.subplots_adjust(hspace=0.5)
    for i in range (0, len(motivsToCheck), 1):
        data = []
        xTicks= []
        for j in range (0, len(exp_set_ups), 1):   
            df_fimo_temp = df_fimo.loc[(df_fimo['motif_id'].str.contains(str(motivsToCheck[i])))==True]
            #print("Number of seqs with motif", df_fimo_temp.shape)
            df_fimo_temp=df_fimo_temp.drop_duplicates(subset="ID") #drop duplicates damit sequencen in dem motif mehrfach vorkommt nur einmal gezeahlt werden
            #print("Number of seqs with motif without duplicates", df_fimo_temp.shape)
            temp_df_IDs_seqs_reg_labels_hasMotiv=df_IDs_seqs_reg_labels[df_IDs_seqs_reg_labels["ID"].isin(df_fimo_temp["ID"].to_list())] #to_list() damit reihenfolge egal ist
            #print("Overlap fimo labels", temp_df_IDs_seqs_reg_labels_hasMotiv.shape)
            np_IDs_seqs_reg_labels_hasMotiv=temp_df_IDs_seqs_reg_labels_hasMotiv["mean_"+str(exp_set_ups[j])].to_numpy()
            data.append(np_IDs_seqs_reg_labels_hasMotiv)
            xTicks.append(str(exp_set_ups[j])+" +")


            df_fimo_temp = df_fimo.loc[(df_fimo['motif_id'].str.contains(str(motivsToCheck[i])))==True]
            df_fimo_temp=df_fimo_temp.drop_duplicates(subset="ID")
            temp_df_IDs_seqs_reg_labels_hasNotMotiv=df_IDs_seqs_reg_labels[~df_IDs_seqs_reg_labels["ID"].isin(df_fimo_temp["ID"].to_list())]
            np_IDs_seqs_reg_labels_hasNotMotiv=temp_df_IDs_seqs_reg_labels_hasNotMotiv["mean_"+str(exp_set_ups[j])].to_numpy()
            data.append(np_IDs_seqs_reg_labels_hasNotMotiv)
            xTicks.append(str(exp_set_ups[j])+" -")

        
        # differences
        np_IDs_seqs_reg_labels_hasMotiv_diff = data[0] - data[2]
        data.append(np_IDs_seqs_reg_labels_hasMotiv_diff)
        xTicks.append("Difference +")

        np_IDs_seqs_reg_labels_hasNotMotiv_diff = data[1] - data[3]
        data.append(np_IDs_seqs_reg_labels_hasNotMotiv_diff)
        xTicks.append("Difference -")

            
            
        #print(data)
        box=axs[i].boxplot(data, showfliers=False, patch_artist=True)
        colors = ['cyan', 'r', "cyan", "r", "cyan", "r"]
        for patch, color in zip(box['boxes'], colors):
            patch.set_facecolor(color)
        #plt.set_ylim=(-2, 100)

                #plt siginificance bars
        if stats.kruskal(data[0], data[1])[1]<0.01:
            level=1
            # Get the y-axis limits
            bottom, top = axs[i].get_ylim()
            y_range = top - bottom

            # Plot the bar
            bar_height = (y_range * 0.002 * level) + top
            bar_tips = bar_height - (y_range * 0.02)
            axs[i].plot(
                [1, 1, 2, 2],
                [bar_tips, bar_height, bar_height, bar_tips], lw=1, c='k'
            )

            text_height = bar_height + (y_range * 0.01)
            axs[i].text((1 + 2) * 0.5, text_height, "*", ha='center', va='bottom', c='k')   


        
        #plt siginificance bars
        if stats.kruskal(data[2], data[3])[1]<0.01:
            level=1
            # Get the y-axis limits
            bottom, top = axs[i].get_ylim()
            y_range = top - bottom

            # Plot the bar
            bar_height = (y_range * 0.002 * level) + top
            bar_tips = bar_height - (y_range * 0.02)
            axs[i].plot(
                [3, 3 , 4, 4],
                [bar_tips, bar_height, bar_height, bar_tips], lw=1, c='k'
            )

            text_height = bar_height + (y_range * 0.01)
            axs[i].text((3 + 4) * 0.5, text_height, "*", ha='center', va='bottom', c='k')   


        #plt siginificance bars
        if stats.kruskal(data[4], data[5])[1]<0.01:
            level=1
            # Get the y-axis limits
            bottom, top = axs[i].get_ylim()
            y_range = top - bottom

            # Plot the bar
            bar_height = (y_range * 0.002 * level) + top
            bar_tips = bar_height - (y_range * 0.02)
            axs[i].plot(
                [5, 5, 6, 6],
                [bar_tips, bar_height, bar_height, bar_tips], lw=1, c='k'
            )

            text_height = bar_height + (y_range * 0.01)
            axs[i].text((5 + 6) * 0.5, text_height, "*", ha='center', va='bottom', c='k')   
        

        
         
        axs[i].set_xticks([1,2,3,4,5,6], xTicks, rotation=90)
        axs[i].set_ylabel("experimental activity")
        axs[i].set_title(" n("+str(motivsToCheck[i])+")="+str(data[0].shape[0])+"; "+'n(other)='+str(data[1].shape[0]))    

        fig.tight_layout()         
        plt.savefig(str(out)+str(exp_set_ups)+".svg")
        #plt.close()

        



        
if __name__ == "__main__":
    cli()
