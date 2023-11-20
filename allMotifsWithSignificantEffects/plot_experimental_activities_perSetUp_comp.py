"""
Description:            This scripts plots Boxplots without outlyers of experimental STARseq activity and differences for Sequences containing a motif of interest or not. Note that the motifs should be provided by the user
                        using the "motivsToCheck" array. The array can have as many motifs as wanted.


Output:                  Boxplots for activities with and without motifs of interest for TeloHEAC_CTRL and TeloHEAC treated with IL1b after 6h

example bash command:   python ./allMotifsWithSignificantEffects/plot_experimental_activities_perSetUp_comp.py 
"""
from wordcloud import WordCloud
from email.policy import default
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from scipy import stats
import click
import logomaker
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
    default="all",
    help=" How output figure file will start",
)
@click.option(
    "--outputPWMs",
    "outputPWM",
    required=True,
    multiple=False,
    type=str,
    default="TeloHAEC_CTRL_sigPWMs.txt",
    help=" How output figure file will start",
)
@click.option(
    "--setUps",
    "exp_set_ups",
    required=True,
    multiple=True,
    type=str,
    default= [ "TeloHAEC_CTRL", "TeloHAEC_IL1b_6h" ],
    help="set ups to compare. choose only 2",
)
@click.option(
    "--PWMfile",
    "PWM_file",
    required=True,
    multiple=False,
    type=str,
    default= "condensed_test.txt",
    help="contains all PWMs to be tested",
)
@click.option(
    "--dataBase",
    "data_base_file",
    required=True,
    multiple=False,
    type=str,
    default="JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme_nice.txt",
    help="e.g. JASPAR",
)
@click.option(
    "--tomtomSig",
    "tomtom_qual_level",
    required=True,
    multiple=False,
    default=0.01,
    help="alpha for sig testing for tomtom matches",
)
@click.option(
    "--number-motifs_to_plot",
    "number_motifs_to_plot",
    required=True,
    multiple=False,
    default=1,
    help="alpha for sig testing for tomtom matches",
)
@click.option(
    "--fimoQual",
    "fimo_qual_level",
    required=True,
    multiple=False,
    default=0.05,
    help="alpha for sig testing for tomtom matches",
)
def cli(seq_file, fimo_file, activity_file, out, exp_set_ups, PWM_file, data_base_file, tomtom_qual_level, number_motifs_to_plot, outputPWM, fimo_qual_level ):




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
    #plot experimental activities

    #initialize outfile with significant PWMs in outPWM  set up
    outPWM=open(outputPWM, "w")
    background = 0.25
    
    outPWM.write('MEME version 4\n\n')
    outPWM.write('ALPHABET= ACGT\n\n')
    outPWM.write('strands: + -\n\n')
    outPWM.write('Background letter frequencies (from unknown source):\n')
    outPWM.write('A %.3f C %.3f G %.3f T %.3f\n\n' % (background, background, background, background))
    
    
    fig, axs = plt.subplots(nrows=len(motivsToCheck), ncols=3+number_motifs_to_plot*2+1, sharex=False, sharey=False, figsize=((3+number_motifs_to_plot*2+1)*8 ,4*len(motivsToCheck)))
    if len(motivsToCheck)==1: # bug repair, da sonst bei nur einem motif ax[i] mit i=0 nicht gefunden wird. Dann ist axs keine liste sondern nur ein ding
        axs=[axs]
    fig.subplots_adjust(hspace=0.5)
    for i in range (0, len(motivsToCheck), 1):
        # motif logo
            
        # read PWM file
        data_base=open(PWM_file, "r")
        data_base_entries=data_base.read().split("MOTIF")
        data_base.close()

        # look for particular motif and write them to output file
        for p in range (0, len(data_base_entries), 1):
            #print(data_base_entries[p].split("\n")[0].split()[0])
            if data_base_entries[p].split("\n")[0].split()[0] == str(motivsToCheck[i]):
                
                plot_logo(data_base_entries[p], axs, i, 1, "CNN motif " +data_base_entries[p].split("\n")[0].split()[0])
                
                # look for JASPAR matches

                # write PWM in meme file
                tmp_file=open("tmp.txt", "w")
                background = 0.25
    
                tmp_file.write('MEME version 4\n\n')
                tmp_file.write('ALPHABET= ACGT\n\n')
                tmp_file.write('strands: + -\n\n')
                tmp_file.write('Background letter frequencies (from unknown source):\n')
                tmp_file.write('A %.3f C %.3f G %.3f T %.3f\n\n' % (background, background, background, background))
                append_to_meme_file(data_base_entries[p], tmp_file)
                tmp_file.close()

                
                # run tomtom gaiannst JASPAR
                cmd = 'export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.4:$PATH'
                os.system(cmd)
                cmd = 'tomtom -no-ssc -oc . --verbosity 1 -text -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10.0 tmp.txt %s > %s' % (data_base_file , "tmp_tomtom")
                os.system(cmd)
                
                # select most siginifcant tomtom matches with pvalue < 0.01
                if os.path.getsize("tmp_tomtom")!=0:
                    df_tomtom=pd.read_csv("tmp_tomtom", sep="\t", skipfooter=4)
                    df_tomtom_sig=df_tomtom.loc[df_tomtom['q-value'] < tomtom_qual_level]
                    df_tomtom_sig=df_tomtom_sig.sort_values(by=['q-value'])
                    df_tomtom_sig=df_tomtom_sig.reset_index(drop=True)

                    # look for corresponding PWSs in JASPAR
                        # read data base entries
                    JASPAR=open(data_base_file, "r")
                    JASPAR_entries=JASPAR.read().split("MOTIF")
                    JASPAR.close()
                    #print(len(JASPAR_entries))
                    #print(JASPAR_entries[3])

                    # look for particular motif and plot PWM logo
                    otherMotifs=[]
                    for index, row in df_tomtom_sig.iterrows():
                        #additional_match=0
                        tmp_motif=row['Target_ID']
                        if int(index)<number_motifs_to_plot:
                            for g in range (0, len(JASPAR_entries), 1):
                                if JASPAR_entries[g].split("\n")[0].split()[0] == tmp_motif:
                                    plot_logo(JASPAR_entries[g], axs, i, 3+index*2, "data base match: "+JASPAR_entries[g].split("\n")[0].split()[0])
                                    break
                        else:
                            for g in range (0, len(JASPAR_entries), 1):
                                if JASPAR_entries[g].split("\n")[0].split()[0] == tmp_motif:
                                    otherMotifs.append(JASPAR_entries[g].split("\n")[0].split()[0][9:])
                                    #axs[i,2+motifs_to_plot*2+1].text(0.2 , -0.3+(index-1)*0.05, JASPAR_entries[g].split("\n")[0].split()[0], size=6) #((index-1)//15)*0.4
                                    #additional_match=additional_match+1
                                    break
                    if len(otherMotifs)>0:
                        axs[i,2+number_motifs_to_plot*2+1].imshow(WordCloud(collocations=False, background_color="white").generate(str(otherMotifs)))
                    axs[i,2+number_motifs_to_plot*2+1].set_axis_off()
                    axs[i,2+number_motifs_to_plot*2+1].set_title("other data base matches:") 
                break 

        
        
        # experimental activities
        data = []
        xTicks= []
        for j in range (0, len(exp_set_ups), 1):   
            df_fimo_temp = df_fimo.loc[(df_fimo['motif_id'].astype(str).str.match(str(motivsToCheck[i])))==True]
            df_fimo_temp= df_fimo_temp.loc[df_fimo_temp["q-value"] < fimo_qual_level]
            #print("Number of seqs with motif", df_fimo_temp.shape)
            df_fimo_temp=df_fimo_temp.drop_duplicates(subset="ID") #drop duplicates damit sequencen in dem motif mehrfach vorkommt nur einmal gezeahlt werden
            #print("Number of seqs with motif without duplicates", df_fimo_temp.shape)
            temp_df_IDs_seqs_reg_labels_hasMotiv=df_IDs_seqs_reg_labels[df_IDs_seqs_reg_labels["ID"].isin(df_fimo_temp["ID"].to_list()) ] #to_list() damit reihenfolge egal ist
            
            #print("Overlap fimo labels", temp_df_IDs_seqs_reg_labels_hasMotiv.shape)
            np_IDs_seqs_reg_labels_hasMotiv=temp_df_IDs_seqs_reg_labels_hasMotiv["mean_"+str(exp_set_ups[j])].to_numpy()
            data.append(np_IDs_seqs_reg_labels_hasMotiv)
            xTicks.append(str(exp_set_ups[j])+" +")


            df_fimo_temp = df_fimo.loc[(df_fimo['motif_id'].astype(str).str.match(str(motivsToCheck[i])))==True]
            df_fimo_temp= df_fimo_temp.loc[df_fimo_temp["q-value"] < fimo_qual_level]
            df_fimo_temp=df_fimo_temp.drop_duplicates(subset="ID")
            temp_df_IDs_seqs_reg_labels_hasNotMotiv=df_IDs_seqs_reg_labels[~df_IDs_seqs_reg_labels["ID"].isin(df_fimo_temp["ID"].to_list()) ]
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
        box=axs[i,0].boxplot(data, showfliers=False, patch_artist=True)
        colors = ['cyan', 'r', "cyan", "r", "cyan", "r"]
        for patch, color in zip(box['boxes'], colors):
            patch.set_facecolor(color)
        #plt.set_ylim=(-2, 100)




        level=1
        # Get the y-axis limits
        bottom, top = axs[i,0].get_ylim()
        y_range = top - bottom
        bar_height = top #(y_range * 0.0002 * level) + top
        bar_tips = bar_height - (y_range * 0.01)
        text_height = bar_height + (y_range * 0.005)

        for j in range (0, (len(exp_set_ups)+1)*2, 2): #+1 wegen diff plots
            #plt siginificance bars
            #level=1
            # Get the y-axis limits
            #bottom, top = axs[i,0].get_ylim()
            #y_range = top - bottom

            # Plot the bar
            #bar_height = top #(y_range * 0.0002 * level) + top
            #bar_tips = bar_height - (y_range * 0.01)
            axs[i,0].plot(
                [j+1, j+1, j+2, j+2],
                [bar_tips, bar_height, bar_height, bar_tips], lw=1, c='k'
            )

            #text_height = bar_height + (y_range * 0.005)
            if stats.kruskal(data[j], data[j+1])[1]<0.01:
                axs[i,0].text((j+1 + j+2) * 0.5, text_height, "*", ha='center', va='bottom', c='k')  

                ### write significant PWMs for selected set up in outpu meme file
                if j==4:
                    for p in range (0, len(data_base_entries), 1):
                        #print(data_base_entries[p].split("\n")[0].split()[0])
                        if data_base_entries[p].split("\n")[0].split()[0] == str(motivsToCheck[i]):
                            append_to_meme_file(data_base_entries[p], outPWM)


            else:
                axs[i,0].text((j+1 + j+2) * 0.5, text_height, "n.s.", ha='center', va='bottom', c='k')  
        
        
         
        axs[i,0].set_xticks([1,2,3,4,5,6], xTicks, rotation=90)
        axs[i,0].set_ylabel("experimental activity")
        axs[i,0].set_title(" n(MOTIF"+str(motivsToCheck[i])+")="+str(data[0].shape[0])+"; "+'n(other)='+str(data[1].shape[0])) #ha

    #plot logo   





    fig.tight_layout()         
    plt.savefig(str(out)+str(exp_set_ups)+".svg")
    outPWM.close()
        #plt.close()

def plot_logo(data_base_entry, axs, axsrow, axscolumn, title):
    temp_PWM=[[],[],[],[]]
    for k in range (2, len(data_base_entry.split("\n"))-2, 1): #2 bis -2 damit header und fuss weg fallen
        if data_base_entry.split("\n")[k].split()[0]!="URL":
            temp_PWM[0].append(float(data_base_entry.split("\n")[k].split()[0]))
            temp_PWM[1].append(float(data_base_entry.split("\n")[k].split()[1]))
            temp_PWM[2].append(float(data_base_entry.split("\n")[k].split()[2]))
            temp_PWM[3].append(float(data_base_entry.split("\n")[k].split()[3]))
                                
    background = 0.25
    cwm_fwd = np.array(temp_PWM) - background
    cwm_rev = cwm_fwd[::-1, ::-1]

    cwm_fwd=cwm_fwd.tolist()
    cwm_rev=cwm_rev.tolist()
                                    

    df_PWM_f=pd.DataFrame()
    df_PWM_f["A"]=cwm_fwd[0]
    df_PWM_f["C"]=cwm_fwd[1]
    df_PWM_f["G"]=cwm_fwd[2]
    df_PWM_f["T"]=cwm_fwd[3]

    df_PWM_r=pd.DataFrame()
    df_PWM_r["A"]=cwm_rev[0]
    df_PWM_r["C"]=cwm_rev[1]
    df_PWM_r["G"]=cwm_rev[2]
    df_PWM_r["T"]=cwm_rev[3]

                                    #print (df_PWM)
                                

                                # create forward Logo object
    PWM_logo = logomaker.Logo(df_PWM_f,
                            shade_below=.5,
                            fade_below=.5,
                            ax=axs[axsrow,axscolumn]
                            )
                                
    # style using Logo methods
    PWM_logo.style_spines(visible=False)
    PWM_logo.style_spines(spines=['left', 'bottom'], visible=True)
    PWM_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

    # style using Axes methods
    PWM_logo.ax.set_ylabel("position frequency", labelpad=-1)
    PWM_logo.ax.xaxis.set_ticks_position('none')
    PWM_logo.ax.xaxis.set_tick_params(pad=-1)   
    PWM_logo.ax.set_title(title + " fwd") 

    # create rev Logo object
    PWM_logo = logomaker.Logo(df_PWM_r,
                            shade_below=.5,
                            fade_below=.5,
                            ax=axs[axsrow,axscolumn+1]
                            )

    # style using Logo methods
    PWM_logo.style_spines(visible=False)
    PWM_logo.style_spines(spines=['left', 'bottom'], visible=True)
    PWM_logo.style_xticks(rotation=90, fmt='%d', anchor=0) 
    
    # style using Axes methods
    PWM_logo.ax.set_ylabel("position frequency", labelpad=-1)
    PWM_logo.ax.xaxis.set_ticks_position('none')
    PWM_logo.ax.xaxis.set_tick_params(pad=-1)    
    PWM_logo.ax.set_title(title + " rev")

def append_to_meme_file(data_base_entry, existingFileName):
    existingFileName.write("MOTIF " + data_base_entry)

        
if __name__ == "__main__":
    cli()
