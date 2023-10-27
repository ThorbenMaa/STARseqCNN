"""
Description:            see readme

Inputs:                 see below (click otions)


Outputs:                heatmap 

Example commad:         python processReports_withMotifClusteringForSpecAna.py --report_folder reportMinna_v2
                        python processReports_withMotifClusteringForSpecAna.py --report_folder ./ --report_files report_diffHEPG2_untreated_pilot_vs_TeloHEAC_CTRL


"""

#set seed and load dependencies
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import click
from tqdm import tqdm
import math

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
            "report_ccell_3T3_undiff_CTRL",
            ],
    help="set ups to consider here. default is all",
)
@click.option(
    "--q_val",
    "q_val",
    required=True,
    multiple=False,
    type=float,
    default=0.01,
    help="set ups to consider here. default is all",
)
def cli(report_folder, report_files, q_val):
    #print(report_files)
    pwd=os.getcwd()
    #df_all=pd.DataFrame()

    all_reports_list=[]
    for i in range (0, len(report_files), 1):
        # set dir
        os.chdir(pwd+str("/")+report_folder+str("/")+report_files[i])

        #read motifs.html
        df_temp0=pd.DataFrame()
        df_temp0=pd.read_html("motifs.html", flavor='html5lib')[0]
        #print(df_temp0)

        # rename motifs not matched with data base
        unknownMotifNumber = 1
        for index, row in df_temp0.iterrows():
            if str(row["match0"])=="nan":
                print (df_temp0.at[index,"match0"])
                print ("unknown motif identified in "+report_files[i])
                df_temp0.at[index, "match0"] = "unkown motif #"+str(unknownMotifNumber)+" "+report_files[i]
                df_temp0.at[index, "match1"] = "unkown motif #"+str(unknownMotifNumber)+" "+report_files[i]
                df_temp0.at[index, "match2"] = "unkown motif #"+str(unknownMotifNumber)+" "+report_files[i]
                unknownMotifNumber = unknownMotifNumber + 1
                print (df_temp0.at[index, "match0"])

        #create
        df_temp0_list=df_temp0.values.tolist()
        print(len(df_temp0_list))
        all_reports_list.append(df_temp0_list)
    print(len(all_reports_list))

    


    # cluster JASPAR motifs
    motif_cluster_list=[]

    # iterate ofer each row in motifs.html of all set ups
    for reportNumber in range (0, len(all_reports_list), 1):
        for row in range(0, len(all_reports_list[reportNumber]), 1):
            
            # create first cluster
            if len(motif_cluster_list) == 0:
                motif_cluster_list.append([])
                for motifNumberInRow in [4,7,10]:
                    #print(all_reports_list[reportNumber][row][motifNumberInRow])
                    if str(all_reports_list[reportNumber][row][motifNumberInRow])!="nan" and all_reports_list[reportNumber][row][motifNumberInRow +1] < q_val: #qval <0.1
                        motif_cluster_list[0].append(all_reports_list[reportNumber][row][motifNumberInRow]) 
                if len(motif_cluster_list[0])==0: # falls alle schlechte qval haben
                    del motif_cluster_list[0]
            
            # create further clusters
            else:
                # check whether one of the motifs of the row is already part of an existing cluster                
                
                added_to_any_cluster=False
                for existing_ClusterNumber in range (0, len(motif_cluster_list), 1):
                    added_to_this_existing_cluster=False
                    for existing_motifNumber in range (0, len(motif_cluster_list[existing_ClusterNumber]), 1):
                        for row_motifNumber  in [4,7,10]: # motif names are in columns 4,7, and 10 in motifs.html
                            if (motif_cluster_list[existing_ClusterNumber][existing_motifNumber]==all_reports_list[reportNumber][row][row_motifNumber]):
                                if str(all_reports_list[reportNumber][row][4]) != "nan" and all_reports_list[reportNumber][row][4+1] < q_val: #quality factor <0.1
                                    motif_cluster_list[existing_ClusterNumber].append(all_reports_list[reportNumber][row][4])
                                if str(all_reports_list[reportNumber][row][7]) != "nan" and all_reports_list[reportNumber][row][7+1] < q_val:
                                    motif_cluster_list[existing_ClusterNumber].append(all_reports_list[reportNumber][row][7])
                                if str(all_reports_list[reportNumber][row][10]) != "nan" and all_reports_list[reportNumber][row][10+1] < q_val:
                                    motif_cluster_list[existing_ClusterNumber].append(all_reports_list[reportNumber][row][10])    
                                added_to_this_existing_cluster = True
                                added_to_any_cluster=True 
                                break

                        # actually, if tfmodisco does its job as I hope, each motif should occur exactly one time in all the clusters.
                        # However, I allow here to look whether it would fit into further clusters and group clusters with overlaps later
                        # this is why there is no break for the next for loop 
                        if added_to_this_existing_cluster == True:
                            break # if they are alreday added to this cluster (they can still be added to another cluster)

                # if motifs have no overlap with any existing cluster, a new cluser is added               
                if added_to_any_cluster == False:
                    motif_cluster_list.append([])
                    for motifNumberInRow in [4,7,10]:
                        #print(all_reports_list[reportNumber][row][motifNumberInRow])
                        if str(all_reports_list[reportNumber][row][motifNumberInRow]) != "nan" and all_reports_list[reportNumber][row][motifNumberInRow +1] < q_val: #qval <0.1:
                            motif_cluster_list[-1].append(all_reports_list[reportNumber][row][motifNumberInRow])      
                            #print(all_reports_list[reportNumber][row][motifNumberInRow])
                            #print(type(all_reports_list[reportNumber][row][motifNumberInRow]))  
                    if len(motif_cluster_list[-1])==0: # falls alle schlechte qval haben
                        del motif_cluster_list[-1]
                            




    print ("number of clusters incl. duplicates", len(motif_cluster_list))
    total_motifs=0
    for i in range (0, len(motif_cluster_list), 1):
        #print (len(motif_cluster_list[i]))
        total_motifs = total_motifs + len(motif_cluster_list[i])
    print ("total number of motifs", total_motifs)


    # check whther there are overlaps between clusters and group them togteher if so and delete one of the overlapping clusters
      
        
    clusterNumber = 0
    len_morif_cluster_list = len(motif_cluster_list)
    """
    for i in range (0, len(motif_cluster_list), 1):
        if i==3 or i == 26:    
            print ("cluster list", motif_cluster_list[i])

    """
    print ("check for overlaps in clusters")
    while clusterNumber < len_morif_cluster_list:
        print ("clusterNumber", clusterNumber)
        index_list_clusters_to_delete=[]
        for motifNumber in range (0, len(motif_cluster_list[clusterNumber]), 1):
            for clusterNumber2 in range (clusterNumber, len(motif_cluster_list), 1):
                if clusterNumber2 not in index_list_clusters_to_delete: # to avoid checking this cluster again
                    doubleCluster=False
                    for motifNumber2 in range (0, len(motif_cluster_list[clusterNumber2]), 1):
                            
                        # exclude comparison of identical clusters
                        if clusterNumber != clusterNumber2:

                            # add all motifs from cluster2 to cluster1 if there is overlap and delete cluster 2
                            if motif_cluster_list[clusterNumber][motifNumber] == motif_cluster_list[clusterNumber2][motifNumber2]:
                                for motifNumber3 in range (0, len(motif_cluster_list[clusterNumber2]), 1):
                                    if motif_cluster_list[clusterNumber2][motifNumber3] != np.nan:
                                        motif_cluster_list[clusterNumber].append(motif_cluster_list[clusterNumber2][motifNumber3])
                                    
                                doubleCluster=True
                                break
    
                    if doubleCluster==True:
                        index_list_clusters_to_delete.append(clusterNumber2) #cannot be deleted directly as this would mess up the for loops
                
        # delete overlapping cluster
        index_list_clusters_to_delete = list(set(index_list_clusters_to_delete)) #remove duplicates 
        #print ("duoble entries", index_list_clusters_to_delete)
        if len(index_list_clusters_to_delete)!=0:
            for index in range(len(index_list_clusters_to_delete)-1, -1, -1): #has to be backwards
                del motif_cluster_list[index_list_clusters_to_delete[index]]
        
        print ("number of motifs in cluster incl. duplicates", len(motif_cluster_list[clusterNumber]))

        # update variables for whiile loop
        len_morif_cluster_list = len(motif_cluster_list)
        clusterNumber = clusterNumber + 1 
        
        
        
    # remove duplicates within cluster
    print ("remove duplicates within cluster")
    for clusterNumber in tqdm(range (0, len(motif_cluster_list), 1)):
        #print ("clusterNumber", clusterNumber)
                 
        motifNumber = 0
        len_cluster = len (motif_cluster_list[clusterNumber])
        while motifNumber < len_cluster:
            index_list_motifs_to_delete=[]
            for motifNumber2 in range (motifNumber, len(motif_cluster_list[clusterNumber]), 1):
                    
                # exclude comparison of identical clusters
                if motifNumber != motifNumber2:

                    # check if it occurs at least twice
                    if motif_cluster_list[clusterNumber][motifNumber] == motif_cluster_list[clusterNumber][motifNumber2]:
                        index_list_motifs_to_delete.append(motifNumber2)
            
            # drop duplicates
            #print ("number of duoble entries", len(index_list_motifs_to_delete))
            if len(index_list_motifs_to_delete)!=0:
                for index in range(len(index_list_motifs_to_delete)-1, -1, -1): #has to be backwards
                    #print(index_list_motifs_to_delete[index])
                    del motif_cluster_list[clusterNumber][index_list_motifs_to_delete[index]]

            # update variables for while loop
            motifNumber = motifNumber + 1
            len_cluster = len (motif_cluster_list[clusterNumber])
    
    for i in range (0, len(motif_cluster_list), 1):
        print ("cluster", i)
        print ("number of motifs in cluster excl. duplicates", len(motif_cluster_list[i]))

    # add up number of seqlets per cluster and set up
    seqlets_per_report_list=[]

    print ("add up number of seqlets per cluster and set up")
    # iterate over motifs html table of articulat exp set up
    for reportNumber in tqdm(range (0, len(all_reports_list), 1)):
        seqlets_per_cluster_in_report_all_clusters=[]
        
        # iterate over clusters and motifs within clusters 
        for clusterNumber in range (0, len(motif_cluster_list), 1):  
            seqlets_per_cluster_in_report = 0
            for motifNumber in range (0, len(motif_cluster_list[clusterNumber]), 1):

                # iterate over rows in html table
                for row in range(0, len(all_reports_list[reportNumber]), 1):
                    
                    # iterate over motifs per row
                    for row_motifNumber  in [4,7,10]:

                        # add up seqlets occuring in this cluster
                        if motif_cluster_list[clusterNumber][motifNumber]==all_reports_list[reportNumber][row][row_motifNumber]:
                            seqlets_per_cluster_in_report = seqlets_per_cluster_in_report + int(all_reports_list[reportNumber][row][1]) # seqlets column
                
            seqlets_per_cluster_in_report_all_clusters.append(seqlets_per_cluster_in_report)

        seqlets_per_report_list.append(seqlets_per_cluster_in_report_all_clusters)
        print (seqlets_per_cluster_in_report_all_clusters)

    print ("plotting...")
    ## create heatmap
    # y ticks
    y_list=[]
    for i in range (0, len(report_files), 1):
        y_list.append(str(report_files[i])[7:])


    # create 2d numpy array
    value_list=seqlets_per_report_list
    value_list=np.array(value_list)
    value_list=np.nan_to_num(value_list, copy=True, nan=0, posinf=None, neginf=None).T
    value_list=np.log(value_list)

    TF_list=motif_cluster_list

    # plot data
    f = plt.figure(figsize=(24, 8))
    plt.matshow(value_list, fignum=f.number)

    plt.xticks(np.arange(len(y_list)), labels=y_list, fontsize=12 , rotation=90)
    
    for i in range (0, len(TF_list), 1):
        TF_list[i]=str(TF_list[i])
    
    #TF_list = [ label.replace(' ', '\n') for label in TF_list ]
    
    plt.yticks(np.arange(len(TF_list)), labels=TF_list, fontsize=12)

    # color bar legend
    cb = plt.colorbar(label="log(seqlets)")
    cb.ax.tick_params(labelsize=12)

    plt.title("log(number of seqlets) per TF motif and experimental set up")
    plt.xlabel("experimental set ups")
    plt.ylabel("TF mtif group")
    plt.tight_layout()
    #plt.show()
    os.chdir(pwd+str("/")+report_folder)
    plt.savefig("TFM_heatmap_withClustering_TeloHEAC_IL1B_6h_vs_24h.svg")
    plt.close()

            
if __name__ == "__main__":
    cli()



