

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os
os.environ["TF_CPP_MIN_LOG_LEVEL"]="2" #to supress warnings with loading tf
import tensorflow as tf
from tensorflow import keras
import random
import click
random.seed(42)#to make things reproducable, 42 is abitrary
np.random.seed(42)#to make things reproducable, 42 is abitrary
tf.random.set_seed(42)#to make things reproducable, 42 is abitrary





@click.command()
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
    "--seqFile",
    "seq_file",
    required=True,
    multiple=False,
    type=str,
    default="starrseq-all-final-toorder_oligocomposition.csv",
    help="e.g. starrseq-all-final-toorder_oligocomposition.csv",
)
@click.option(
    "--sequence_length",
    "sequence_length",
    required=True,
    multiple=False,
    type=int,
    default=198,
    help="length of input seqs",
)
@click.option(
    "--batch_size",
    "batch_size",
    required=True,
    multiple=False,
    type=int,
    default=128,
    help="batch size",
)
@click.option(
    "--model",
    "model_name",
    required=True,
    multiple=False,
    type=str,
    default="TM_s_model_early_stop",
    help="name of model",
)
@click.option(
    "--modelFolder",
    "model_folder",
    required=True,
    multiple=False,
    type=str,
    default="models/",
    help="path to model",
)
def cli(activity_file, seq_file, sequence_length, batch_size, model_name, model_folder):
    
    #import labels
    df_list=[]
    for i in range (0, len(activity_file), 1):
        df_list.append(pd.read_csv(activity_file[i], sep="\t", decimal=',', low_memory=False))
    df_IDs_reg_labels = pd.concat(df_list, axis=0)
    df_IDs_reg_labels=df_IDs_reg_labels.drop_duplicates()
    #parameters


    #import sequences
    df_IDs_Sequences=pd.read_csv(seq_file, sep=",",low_memory=False)

    #remove ">" as first character from ID column
    df_IDs_Sequences["name"]=df_IDs_Sequences["name"].str[1:]

    #convert seqs to list and transform to one hot ecnoceded seqs; add to imported sequences data frame
    sequence_list=df_IDs_Sequences['enhancer'].to_list()
    sequences_tok=[]

    #one-hot-encode sequences and remove sequences that are not of correct length for some reason (only a few)
    for i in range (0, len(sequence_list), 1): 
        if  len(sequence_list[i])==sequence_length:
            sequences_tok.append(one_hot_encode(sequence_list[i]))
        else:
            sequences_tok.append(np.nan)
    df_IDs_Sequences['Seq one hot encoded'] = sequences_tok

    #drop duplicates
    df_IDs_Sequences=df_IDs_Sequences.dropna()

    #merge data frames on name/oligo columns
    df_IDs_Sequences=df_IDs_Sequences.rename(columns={'name': 'ID'})
    df_IDs_reg_labels=df_IDs_reg_labels.rename(columns={'Oligo': 'ID'})
    df_IDs_seqs_reg_labels=pd.merge(df_IDs_Sequences, df_IDs_reg_labels, on=["ID"])

    #select sequences from data set (here, no filtering is done)
    df_IDs_seqs_reg_labels_all=df_IDs_seqs_reg_labels

    #create sequence tensor 
    input_seq_all=tf.convert_to_tensor(df_IDs_seqs_reg_labels_all["Seq one hot encoded"].to_list())

    #load model
    model=keras.models.load_model(str(model_folder)+str(model_name))

    #predictions_ref=model.predict(input_seq_all, batch_size=batch_size, verbose=2)

    #define parameters
    numberOfSeqs=np.asarray(input_seq_all).shape[0]
    seq_length=np.asarray(input_seq_all).shape[1]

    #initialize arrays with hot-hot-encoded sequences as booleans for tfmodisco
    seqs_tfmodisco_format=np.zeros([numberOfSeqs, seq_length, 4])
    for i in range (0, numberOfSeqs, 1):
        for j in range (0, seq_length, 1):
            for k in range (0, 4, 1):
                if input_seq_all[i,j,k]==1:
                    seqs_tfmodisco_format[i, j, k]=True
                else:
                    seqs_tfmodisco_format[i, j, k]=False

    #initialize arrays with hypothetical importance scores for tfmodisco
    hypothetical_contribution_scores_mean_cell_3T3_diff_CTRL=np.zeros([numberOfSeqs, seq_length, 4])
    hypothetical_contribution_scores_mean_ccell_3T3_undiff_CTRL=np.zeros([numberOfSeqs, seq_length, 4])
    hypothetical_contribution_scores_mean_cell_3T3_undiff_TGFB=np.zeros([numberOfSeqs, seq_length, 4])
    hypothetical_contribution_scores_mean_RAW_CTRL=np.zeros([numberOfSeqs, seq_length, 4])
    hypothetical_contribution_scores_mean_RAW_IL1B=np.zeros([numberOfSeqs, seq_length, 4])
    hypothetical_contribution_scores_mean_RAW_TGFB=np.zeros([numberOfSeqs, seq_length, 4])
    hypothetical_contribution_scores_mean_TeloHAEC_CTRL=np.zeros([numberOfSeqs, seq_length, 4])
    hypothetical_contribution_scores_mean_TeloHAEC_IL1b_24h=np.zeros([numberOfSeqs, seq_length, 4])
    hypothetical_contribution_scores_mean_TeloHAEC_IL1b_6h=np.zeros([numberOfSeqs, seq_length, 4])
    hypothetical_contribution_scores_mean_HASMC_untreatedPilot=np.zeros([numberOfSeqs, seq_length, 4])
    hypothetical_contribution_scores_mean_HASMC_Chol=np.zeros([numberOfSeqs, seq_length, 4])
    hypothetical_contribution_scores_mean_HepG2_untreatedPilot=np.zeros([numberOfSeqs, seq_length, 4])

    #fill arrays with hypothetical importance scores for tfmodisco
    for i in range (0, numberOfSeqs, 1):#iterate over numberOfSeqs
        print(str(i)+" of " +str(numberOfSeqs))
        listOfMutSeqs=[]
        for j in range (0, seq_length, 1):#iterate over seq_length
            temp_seq=input_seq_all[i].numpy() #weil tf tensors nicht veraenderlich sind
            for k in range (0, 4, 1):#set all positions to "0"
                if temp_seq[j,k]==1:
                    temp_seq[j,k]=0
                
            for k in range (0, 4, 1):#systematically fill positions with "1" to generate all possible SNV sequences for each particular sequence
                temp_seq[j,k]=1
                listOfMutSeqs.append(tf.convert_to_tensor(temp_seq))
                temp_seq[j,k]=0
            
                #parse all possible SNV sequences for each particular sequence to preloaded model and calculate predictions
        tensorOfMutSeqs=tf.convert_to_tensor(listOfMutSeqs)
        predictions_alt=model.predict(tensorOfMutSeqs, batch_size=batch_size, verbose=0)
            
        #fill arrays with hypothetical importance scores for tfmodisco  
        for j in range (0, seq_length, 1):
            for k in range (0, 4, 1):
                hypothetical_contribution_scores_mean_cell_3T3_diff_CTRL[ i, j ,k] = predictions_alt[j*4+k,0]
                hypothetical_contribution_scores_mean_ccell_3T3_undiff_CTRL[ i, j ,k] = predictions_alt[j*4+k,1]
                hypothetical_contribution_scores_mean_cell_3T3_undiff_TGFB[ i, j ,k] = predictions_alt[j*4+k,2]
                hypothetical_contribution_scores_mean_RAW_CTRL[ i, j ,k] = predictions_alt[j*4+k,3]
                hypothetical_contribution_scores_mean_RAW_IL1B[ i, j ,k] = predictions_alt[j*4+k,4]
                hypothetical_contribution_scores_mean_RAW_TGFB[ i, j ,k] = predictions_alt[j*4+k,5]
                hypothetical_contribution_scores_mean_TeloHAEC_CTRL[ i, j ,k] = predictions_alt[j*4+k,6]
                hypothetical_contribution_scores_mean_TeloHAEC_IL1b_24h[ i, j ,k] = predictions_alt[j*4+k,7]
                hypothetical_contribution_scores_mean_TeloHAEC_IL1b_6h[ i, j ,k] = predictions_alt[j*4+k,8]
                hypothetical_contribution_scores_mean_HASMC_untreatedPilot[ i, j ,k] = predictions_alt[j*4+k,9]
                hypothetical_contribution_scores_mean_HASMC_Chol[ i, j ,k] = predictions_alt[j*4+k,10]
                hypothetical_contribution_scores_mean_HepG2_untreatedPilot[ i, j ,k] = predictions_alt[j*4+k,11]



    #print shapes (number of Seqs x seq length x 4)
    #print(hypothetical_contribution_scores_mean_cell_3T3_diff_CTRL.shape)
    #print(seqs_tfmodisco_format.shape)
    print("data saved to "+str(os.getcwd()))
    #save contribution scores as npz files
    np.savez(str(os.getcwd())+"/"+str(model_name)+'hypothetical_contribution_scores_mean_cell_3T3_diff_CTRL.npz', hypothetical_contribution_scores_mean_cell_3T3_diff_CTRL)
    np.savez(str(os.getcwd())+"/"+str(model_name)+'hypothetical_contribution_scores_mean_ccell_3T3_undiff_CTRL.npz', hypothetical_contribution_scores_mean_ccell_3T3_undiff_CTRL)
    np.savez(str(os.getcwd())+"/"+str(model_name)+'hypothetical_contribution_scores_mean_cell_3T3_undiff_TGFB.npz', hypothetical_contribution_scores_mean_cell_3T3_undiff_TGFB)
    np.savez(str(os.getcwd())+"/"+str(model_name)+'hypothetical_contribution_scores_mean_RAW_CTRL.npz', hypothetical_contribution_scores_mean_RAW_CTRL)
    np.savez(str(os.getcwd())+"/"+str(model_name)+'hypothetical_contribution_scores_mean_RAW_IL1B.npz', hypothetical_contribution_scores_mean_RAW_IL1B)
    np.savez(str(os.getcwd())+"/"+str(model_name)+'hypothetical_contribution_scores_mean_RAW_TGFB.npz', hypothetical_contribution_scores_mean_RAW_TGFB)
    np.savez(str(os.getcwd())+"/"+str(model_name)+'hypothetical_contribution_scores_mean_TeloHAEC_CTRL.npz', hypothetical_contribution_scores_mean_TeloHAEC_CTRL)
    np.savez(str(os.getcwd())+"/"+str(model_name)+'hypothetical_contribution_scores_mean_TeloHAEC_IL1b_24h.npz', hypothetical_contribution_scores_mean_TeloHAEC_IL1b_24h)
    np.savez(str(os.getcwd())+"/"+str(model_name)+'hypothetical_contribution_scores_mean_TeloHAEC_IL1b_6h.npz', hypothetical_contribution_scores_mean_TeloHAEC_IL1b_6h)
    np.savez(str(os.getcwd())+"/"+str(model_name)+'hypothetical_contribution_scores_mean_HASMC_untreatedPilot.npz', hypothetical_contribution_scores_mean_HASMC_untreatedPilot)
    np.savez(str(os.getcwd())+"/"+str(model_name)+'hypothetical_contribution_scores_mean_HASMC_Chol.npz', hypothetical_contribution_scores_mean_HASMC_Chol)
    np.savez(str(os.getcwd())+"/"+str(model_name)+'hypothetical_contribution_scores_mean_HepG2_untreatedPilot.npz', hypothetical_contribution_scores_mean_HepG2_untreatedPilot)

    #save sequences as npz file
    np.savez(str(os.getcwd())+"/"+str(model_name)+'Sequences.npz', seqs_tfmodisco_format)




def one_hot_encode(seq): #taken from https://stackoverflow.com/questions/34263772/how-to-generate-one-hot-encoding-for-dna-sequences?rq=3
    mapping = dict(zip("ACGT", range(4)))    
    seq2 = [mapping[i] for i in seq]
    return np.eye(4)[seq2]

            
if __name__ == "__main__":
    cli()