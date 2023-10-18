"""
Description:        

Inputs:                  see below (click otions)


Outputs:                heatmap 

Example commad:         python processReports.py --report_folder reportMinna_v2

"""

#set seed and load dependencies
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import click

@click.command()
@click.option(
    "--report_folder",
    "report_folder",
    required=True,
    multiple=False,
    type=str,
    help="folder with all tf modisco reports",
)
@click.option(
    "--report_files",
    "report_files",
    required=True,
    multiple=True,
    type=str,
    default=["report_TeloHAEC_CTRL",
            "report_TeloHAEC_IL1b_6h",
            "report_TeloHAEC_IL1b_24h",
            "report_RAW_TGFB", 
            "report_RAW_IL1B", 
            "report_RAW_CTRL", 
            "report_HepG2_untreatedPilot", 
            "report_HASMC_untreatedPilot", 
            "report_HASMC_chol", 
            "report_cell_3T3_undiff_TGFB",
            "report_cell_3T3_diff_CTRL",
            "report_ccell_3T3_undiff_CTRL"],
    help="set ups to consider here. defauilt is all",
)
def cli(report_folder, report_files):
    #print(report_files)
    pwd=os.getcwd()
    df_all=pd.DataFrame()

    for i in range (0, len(report_files), 1):
        # set dir
        os.chdir(pwd+str("/")+report_folder+str("/")+report_files[i])

        #read motifs.html
        df_temp0=pd.DataFrame()
        df_temp0=pd.read_html("motifs.html", flavor='html5lib')[0]

        #remove unnecessary columns
        df_temp=df_temp0[["match0", "num_seqlets"]]

        #rename num_seqlet column for mergeing and to keep info about exp set-up
        df_temp.columns=df_temp.columns+"_"+str(report_files[i])

        #name match0 column back
        df_temp=df_temp.rename(columns={"match0"+"_"+str(report_files[i]) : "match0"})

        #drop TF motif that occur more than once, but keep nan values (https://stackoverflow.com/questions/75624770/pandas-drop-duplicates-without-empty-rows)
        df_temp=df_temp[(~df_temp.duplicated(subset=['match0'], keep="first")) | (df_temp['match0'].isnull())]

        # rename TF motifs that do not match with JASPAR data base
        df_temp=df_temp.fillna("unknown motifs")

        #writing to combined data frame
        if i==0:
            print ("adding "+str(report_files[i]))
            df_all=df_temp
            print (df_all)
        else:
            print ("adding "+str(report_files[i]))
            df_all=df_all.merge(df_temp, how="outer", on=["match0"])
            print (df_all)
            
    #sorting an processing table
    #sorting and storing tsv
    sort_list=[]
    for i in range (0, len(report_files), 1):
        sort_list.append("num_seqlets_"+str(report_files[i]))

    df_all=df_all.sort_values(by=sort_list, ascending=False)
    os.chdir(pwd+str("/")+report_folder)
    df_all.to_csv("reports_sorted.tsv", sep ="\t")



    ## create heatmap
    # create 2d numpy array
    value_list=[]
    for i in range (0, len(report_files), 1):
        value_list.append(df_all["num_seqlets_"+str(report_files[i])].values)
    value_list=np.array(value_list)
    value_list=np.nan_to_num(value_list, copy=True, nan=0, posinf=None, neginf=None).T
    value_list=np.log(value_list)

    TF_list=df_all["match0"].to_list()

    # plot data
    f = plt.figure(figsize=(29, 25))
    plt.matshow(value_list, fignum=f.number)

    plt.xticks(np.arange(len(sort_list)), labels=sort_list, fontsize=6 , rotation=90)
    plt.yticks(np.arange(len(TF_list)), labels=TF_list, fontsize=6)

    # color bar legend
    cb = plt.colorbar(label="log(seqlets)")
    cb.ax.tick_params(labelsize=6)

    plt.title("log(number of seqlets) per TF motif and experimental set up")
    plt.show()
    plt.savefig("TFM_heatmap.svg")
    plt.close()

if __name__ == "__main__":
    cli()



