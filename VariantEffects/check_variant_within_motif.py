"""
Description:       
                    

Example commands:   python VariantEffects/check_variant_within_motif.py --output VariantEffects/test
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
    "--fimoFile",
    "fimo_file",
    required=True,
    multiple=False,
    type=str,
    default="fimo_out/fimo.tsv",
    help="fimo file",
)
@click.option(
    "--variantsFile",
    "variants_file",
    required=True,
    multiple=False,
    type=str,
    default="VariantEffects/df_highestVariantEffects_mean_TeloHAEC_IL1b_6h.csv",
    help="e.g. VariantEffects/df_highestVariantEffects_mean_TeloHAEC_IL1b_6h.csv",
)
@click.option(
    "--q_thres",
    "q_thres",
    required=True,
    multiple=False,
    type=float,
    default=0.01,
    help="fimo",
)
@click.option(
    "--p_thres",
    "p_thres",
    required=True,
    multiple=False,
    type=float,
    default=1,
    help="fimo already only gives results with p_val<0.0001",
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
def cli(fimo_file, variants_file, out, q_thres, p_thres):

    #import labels; abs diff activities of seq with ID1 - ID2
    df_activities=pd.read_csv(variants_file, sep=",", low_memory=False)
    print("number of variant pairs :", df_activities.shape[0])
    df_fimo_hits=pd.read_csv(fimo_file, sep="\t",  low_memory=False)

    # select all fimo hits with ID1
    df_fimo_hits1=df_fimo_hits.rename(columns={'sequence_name': 'ID1'})
    #print(df_fimo_hits1)
    df_fimo_hits1=pd.merge(df_fimo_hits1, df_activities, on=["ID1"])
    #print(df_fimo_hits1)
    df_fimo_hits1=df_fimo_hits1[df_fimo_hits1["variantPos"] <= df_fimo_hits1["stop"]]
    df_fimo_hits1=df_fimo_hits1[df_fimo_hits1["start"] <= df_fimo_hits1["variantPos"]]
    df_fimo_hits1=df_fimo_hits1[df_fimo_hits1["q-value"] < q_thres]
    df_fimo_hits1=df_fimo_hits1[df_fimo_hits1["p-value"] < p_thres]
    #print(df_fimo_hits1)

    # select all fimo hits with ID2 (diff sind immer 2 seqs. motif kann theretisch in nur einem der beiden haplotyppes erkannt werden von FIMO)
    df_fimo_hits2=df_fimo_hits.rename(columns={'sequence_name': 'ID2'})
    df_fimo_hits2=pd.merge(df_fimo_hits2, df_activities, on=["ID2"])
    #print(df_fimo_hits2)
    df_fimo_hits2=df_fimo_hits2[df_fimo_hits2["variantPos"] <= df_fimo_hits2["stop"]]
    df_fimo_hits2=df_fimo_hits2[df_fimo_hits2["start"] <= df_fimo_hits2["variantPos"]]
    df_fimo_hits2=df_fimo_hits2[df_fimo_hits2["q-value"] < q_thres]
    df_fimo_hits2=df_fimo_hits2[df_fimo_hits2["p-value"] < p_thres]
    #print(df_fimo_hits2)

    # add fimo hits with ID1 to output
    df_out1=pd.DataFrame()
    df_out1["ID1"] = df_fimo_hits1["ID1"]
    df_out1["ID2"] = df_fimo_hits1["ID2"]
    df_out1["variantPos"] = df_fimo_hits1["variantPos"]
    df_out1["diff_activity mean_TeloHAEC_IL1b_6h"] = df_fimo_hits1["diff_activity mean_TeloHAEC_IL1b_6h"] #Telo HEac hard evcuded but doesn matter, could e every other too
    df_out1["Seq1"] = df_fimo_hits1["Seq1"]
    df_out1["Seq2"] = df_fimo_hits1["Seq2"]
    df_out1["MotifID"] = df_fimo_hits1["motif_id"]
    df_out1["start"] = df_fimo_hits1["start"]
    df_out1["stop"] = df_fimo_hits1["stop"]
    df_out1["strand"] = df_fimo_hits1["strand"]
    df_out1["p-value haplotype effect"] = df_fimo_hits1["P.Value"]
    df_out1["adj. p-value haplotype effect"] = df_fimo_hits1["adj.P.Val"]
    
    #df_out1["TF motif"]=df_fimo_hits1["motif_alt_id"]

    # add fimo hits with ID2 to output. Here, ID2 becomes ID1
    df_out2=pd.DataFrame()
    df_out2["ID1"] = df_fimo_hits2["ID1"]
    df_out2["ID2"] = df_fimo_hits2["ID2"]
    df_out2["variantPos"] = df_fimo_hits2["variantPos"]
    df_out2["diff_activity mean_TeloHAEC_IL1b_6h"] = df_fimo_hits2["diff_activity mean_TeloHAEC_IL1b_6h"]
    df_out2["Seq1"] = df_fimo_hits2["Seq1"]
    df_out2["Seq2"] = df_fimo_hits2["Seq2"]
    df_out2["MotifID"] = df_fimo_hits2["motif_id"]
    df_out2["start"] = df_fimo_hits2["start"]
    df_out2["stop"] = df_fimo_hits2["stop"]
    df_out2["strand"] = df_fimo_hits2["strand"]
    df_out2["p-value haplotype effect"] = df_fimo_hits2["P.Value"]
    df_out2["adj. p-value haplotype effect"] = df_fimo_hits2["adj.P.Val"]
    #df_out2["TF motif"]=df_fimo_hits2["motif_alt_id"]

    # final output
    df_final=pd.concat([df_out1, df_out2])
    #print(df_final)
    print("Number of haplotype effects within tested motifs", df_final.drop_duplicates(subset=["ID1","ID2"]).shape[0])
    df_final=df_final.drop_duplicates(subset=["ID1", "ID2", "variantPos", "diff_activity mean_TeloHAEC_IL1b_6h", "start", "stop", "MotifID"])
    #print("Number of variants within tested motifs", df_final.shape[0])


    df_final.to_csv(str(out)+".csv")
    


   


            
if __name__ == "__main__":
    cli()












