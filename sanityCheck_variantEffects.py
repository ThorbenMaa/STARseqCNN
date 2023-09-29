"""
Description:        This script has two modes. The "train" mode trains different CNNs on STARseq data provided by the Kaikkonen Lab using two different architectures and several different learning rates, with and without augmentation. 
                    Learning rates can be specified using the `learning_rate` array. The model is evaluated on a hold-out test data set. The "load" mode loads a pre-trained model and evaluates it on a
                    hold-out test data set. Pearson correlations between predicted and experimentally determined activities are calculated and scatter plots are generated.

Inputs:             Input 1 are labels for model training and evaluation. Input 2 are sequences for  model training and evaluation. Input 3 defines the mode (see above) and is either "train" or "load".
                    Input 4 is a pre-trained model that should be evaluated. Only important in the "load" mode of this script. Input 5 is the hold-out chromosme for testing. Input 6 has to be "use_aug"
                    to use augmentation for model training or something else if no augmentation should be used. Input7 like input 1
                    

further parameters: can be specified in the parameters section of this script. 

Example commands:   `python sanityCheck_variantEffects.py 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo1.txt starrseq-all-final-toorder_oligocomposition.csv 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo2.txt allseq-CNN_StarSeq_model_Minna_deepSTAR_lr0.01no_aug chr8 all`
                    
Outputs:            Trained CNNs or evalutaion of a pre-trained CNN.

"""

#set seed and load dependencies
from tqdm import tqdm #progress bar
seed_value=1234
import os
os.environ['PYTHONHASHSEED']=str(seed_value)
os.environ["TF_CPP_MIN_LOG_LEVEL"]="2" #to supress warnings with loading tf
import random
random.seed(seed_value)#to make training reproducable, 1234 is abitrary
import numpy as np
np.random.seed(seed_value)#to make training reproducable, 1234 is abitrary
import matplotlib.pyplot as plt
import pandas as pd
import sys
import tensorflow as tf
from tensorflow import keras
tf.random.set_seed(seed_value)
from scipy import stats


def one_hot_encode(seq): #taken from https://stackoverflow.com/questions/34263772/how-to-generate-one-hot-encoding-for-dna-sequences?rq=3
    mapping = dict(zip("ACGT", range(4)))    
    seq2 = [mapping[i] for i in seq]
    return np.eye(4)[seq2]

def complementary(strand): #adapted from https://codereview.stackexchange.com/questions/193766/generating-complementary-dna-sequence
    complementary_strand = ''
    for dna in strand:
        if dna == 'A':
            complementary_strand += 'T'
        elif dna == 'T':
            complementary_strand += 'A'
        elif dna == 'G':
            complementary_strand += 'C'
        elif dna == 'C':
            complementary_strand += 'G'
    return complementary_strand

#parameters
sequence_length=198
learning_rate=[0.001]#, 0.001, 0.0001]
batch_size=128
epochs=150
number_cell_types=12 #also hot encoded below, also needs to be change below when processing tarining and testing data

#import labels
df_IDs_reg_labels=pd.read_csv(sys.argv[1], sep="\t", decimal=',', low_memory=False)
df_IDs_reg_labels2=pd.read_csv(sys.argv[3], sep="\t", decimal=',', low_memory=False)

df_IDs_reg_labels = pd.concat([df_IDs_reg_labels, df_IDs_reg_labels2], axis=0)

df_IDs_reg_labels=df_IDs_reg_labels.drop_duplicates()

#replace 0 with 1 for log transformation later
#df_IDs_reg_labels=df_IDs_reg_labels.replace(0, 1)

#import sequences
df_IDs_Sequences=pd.read_csv(sys.argv[2], sep=",",low_memory=False)

#remove ">" as first character from ID column
df_IDs_Sequences["name"]=df_IDs_Sequences["name"].str[1:]

#convert seqs to list and transform to one hot ecnoceded seqs; add to imported sequences data frame; drop sequences with wrong sequence length
sequence_list=df_IDs_Sequences['enhancer'].to_list()
sequences_tok=[]

for i in range (0, len(sequence_list), 1): 
    if  len(sequence_list[i])==sequence_length:
        sequences_tok.append(one_hot_encode(sequence_list[i]))
    else:
        sequences_tok.append(np.nan)

#add one-hot-encoded to data frame
df_IDs_Sequences['Seq one hot encoded'] = sequences_tok

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
df_diff=pd.DataFrame(columns=["Seq-hot-enc1",
                              "Seq1",
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
                              "Seq-hot-enc2",
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
                              "mean_HepG2_untreatedPilot2"
                              ])

#idetify haplotype pairs
df_IDs_seqs_reg_labels_list=df_IDs_seqs_reg_labels.values.tolist()
#print(df_IDs_seqs_reg_labels.iloc[[1]])
#print (str(df_IDs_seqs_reg_labels_list[1][0]))#)[:-1])

for i in tqdm(range (len(df_IDs_seqs_reg_labels_list))):#tqdm([0, len(df_IDs_seqs_reg_labels_list)]):
    for j in range (0, len(df_IDs_seqs_reg_labels_list), 1):
        if str(df_IDs_seqs_reg_labels_list[i][0])[:-1] == str(df_IDs_seqs_reg_labels_list[j][0])[:-1]:
            if str(df_IDs_seqs_reg_labels_list[i][0]) != str(df_IDs_seqs_reg_labels_list[j][0]): #also "gleiche" seq aber anderer haplotype
                df_diff.loc[len(df_diff)]=[df_IDs_seqs_reg_labels.iloc[[i]]["Seq one hot encoded"].values[0],
                        df_IDs_seqs_reg_labels.iloc[[i]]["enhancer"].values[0],
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
                        df_IDs_seqs_reg_labels.iloc[[j]]["Seq one hot encoded"].values[0],
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
                        df_IDs_seqs_reg_labels.iloc[[j]]["mean_HepG2_untreatedPilot"].values[0]
                        ]#,
                        #columns= df_diff.columns), ignore_index=True)

                  

                            
                
                #df_diff.loc[len(df_diff)]=df_list
    """
    if i==2:
        break
    """
#calc diff
#print (df_diff)
df_diff["diff_activity mean_cell_3T3_diff_CTRL"]=df_diff["mean_cell_3T3_diff_CTRL1"] - df_diff["mean_cell_3T3_diff_CTRL2"]
df_diff["diff_activity mean_ccell_3T3_undiff_CTRL"]=df_diff["mean_ccell_3T3_undiff_CTRL1"] - df_diff["mean_ccell_3T3_undiff_CTRL2"]
df_diff["diff_activity mean_cell_3T3_undiff_TGFB"]=df_diff["mean_cell_3T3_undiff_TGFB1"] - df_diff["mean_cell_3T3_undiff_TGFB2"]
df_diff["diff_activity mean_RAW_CTRL"]=df_diff["mean_RAW_CTRL1"] - df_diff["mean_RAW_CTRL2"]
df_diff["diff_activity mean_RAW_IL1B"]=df_diff["mean_RAW_IL1B1"] - df_diff["mean_RAW_IL1B2"]
df_diff["diff_activity mean_RAW_TGFB"]=df_diff["mean_RAW_TGFB1"] - df_diff["mean_RAW_TGFB2"]
df_diff["diff_activity mean_TeloHAEC_CTRL"]=df_diff["mean_TeloHAEC_CTRL1"] - df_diff["mean_TeloHAEC_CTRL2"]
df_diff["diff_activity mean_TeloHAEC_IL1b_24h"]=df_diff["mean_TeloHAEC_IL1b_24h1"] - df_diff["mean_TeloHAEC_IL1b_24h2"]
df_diff["diff_activity mean_TeloHAEC_IL1b_6h"]=df_diff["mean_TeloHAEC_IL1b_6h1"] - df_diff["mean_TeloHAEC_IL1b_6h2"]
df_diff["diff_activity mean_HASMC_untreatedPilot"]=df_diff["mean_HASMC_untreatedPilot1"] - df_diff["mean_HASMC_untreatedPilot2"]
df_diff["diff_activity mean_HASMC_Chol"]=df_diff["mean_HASMC_Chol1"] - df_diff["mean_HASMC_Chol2"]
df_diff["diff_activity mean_HepG2_untreatedPilot"]=df_diff["mean_HepG2_untreatedPilot1"] - df_diff["mean_HepG2_untreatedPilot2"]

#use abs value
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


#drop duplicates (weil sonst werte dopplet vorkommen, da a-b == b-a bei abs values). Hier koennte man auch alle machen bei subset, aber wahrscheinlichkeit bei 5 zufaellig gleich zu sein ist schon sehr klein
print("before and after drop duplicates")
print (df_diff)
df_diff=df_diff.drop_duplicates(subset=["diff_activity mean_cell_3T3_diff_CTRL", "diff_activity mean_ccell_3T3_undiff_CTRL", "diff_activity mean_cell_3T3_undiff_TGFB", "diff_activity mean_RAW_CTRL", "diff_activity mean_RAW_IL1B"])
print (df_diff)


#test data
if sys.argv[6]=="test_only":
    df_diff_test=df_diff.loc[(df_diff['ID1'].str.contains(str(sys.argv[5])))==True] #ID1 or ID2 doenst matter
elif sys.argv[6]=="all":
    df_diff_test=df_diff
else:
    print("ERRROOROROROROR")

#experimental differences ref and alt tensor
input_label_test=tf.convert_to_tensor([df_diff_test["diff_activity mean_cell_3T3_diff_CTRL"].to_list(),
                                    df_diff_test["diff_activity mean_ccell_3T3_undiff_CTRL"].to_list(),
                                    df_diff_test["diff_activity mean_cell_3T3_undiff_TGFB"].to_list(),
                                    df_diff_test["diff_activity mean_RAW_CTRL"].to_list(),
                                    df_diff_test["diff_activity mean_RAW_IL1B"].to_list(),
                                    df_diff_test["diff_activity mean_RAW_TGFB"].to_list(),
                                    df_diff_test["diff_activity mean_TeloHAEC_CTRL"].to_list(),
                                    df_diff_test["diff_activity mean_TeloHAEC_IL1b_24h"].to_list(),
                                    df_diff_test["diff_activity mean_TeloHAEC_IL1b_6h"].to_list(),
                                    df_diff_test["diff_activity mean_HASMC_untreatedPilot"].to_list(),
                                    df_diff_test["diff_activity mean_HASMC_Chol"].to_list(),
                                    df_diff_test["diff_activity mean_HepG2_untreatedPilot"].to_list()                                    
])

# no log transform but transpose. the model calculates log transformed activities. The solution here is to back transform the output of the model. (lg(a)-lg(b) != lg(a-b)!!)
input_label_test=tf.transpose(input_label_test)
print("diff tensor (label):")
print (input_label_test)
#input_label_test=tf.math.log(tf.transpose(input_label_test))

#sequences test data1
input_seq_test1=tf.cast(tf.convert_to_tensor(df_diff_test["Seq-hot-enc1"].to_list()), tf.int8)
print("seq 1 tensor:")
print(input_seq_test1)

#sequences test data2
input_seq_test2=tf.cast(tf.convert_to_tensor(df_diff_test["Seq-hot-enc2"].to_list()), tf.int8)
print("seq 2 tensor:")
print(input_seq_test2)


model=keras.models.load_model(str(sys.argv[4]))

#calculate predicted labels1
predictions1=model.predict(input_seq_test1, batch_size=batch_size, verbose=2)
print("predictions 1 tensor before and after exp transformation:")
print(predictions1)
predictions1=tf.math.exp(predictions1)
print(predictions1)

#calculate predicted labels2
predictions2=model.predict(input_seq_test2, batch_size=batch_size, verbose=2)
print("predictions 2 tensor before and after exp transformation:")
print(predictions2)
predictions2=tf.math.exp(predictions2)
print(predictions2)

#calculate diff between ref seq and alt seq and take abs values
predictions_diff=tf.abs(predictions1-predictions2)
print("predictions diff tensor:")
print(predictions_diff)

#correlations of predicted and experimental labels for different cell types
print(stats.pearsonr(predictions_diff[:,0], input_label_test[:,0]))
print(stats.pearsonr(predictions_diff[:,1], input_label_test[:,1]))
print(stats.pearsonr(predictions_diff[:,2], input_label_test[:,2]))
print(stats.pearsonr(predictions_diff[:,3], input_label_test[:,3]))
print(stats.pearsonr(predictions_diff[:,4], input_label_test[:,4]))
print(stats.pearsonr(predictions_diff[:,5], input_label_test[:,5]))
print(stats.pearsonr(predictions_diff[:,6], input_label_test[:,6]))
print(stats.pearsonr(predictions_diff[:,7], input_label_test[:,7]))
print(stats.pearsonr(predictions_diff[:,8], input_label_test[:,8]))
print(stats.pearsonr(predictions_diff[:,9], input_label_test[:,9]))
print(stats.pearsonr(predictions_diff[:,10], input_label_test[:,10]))
print(stats.pearsonr(predictions_diff[:,11], input_label_test[:,11]))

#scatterplots of predicted and experimental labels for different cell types
plt.scatter(predictions_diff[:,0], input_label_test[:,0], s=1)
plt.title('diff_mean_cell_3T3_diff_CTRL'+str(stats.pearsonr(predictions_diff[:,0], input_label_test[:,0])))
plt.savefig(str(sys.argv[4])+"_mean_cell_3T3_diff_CTRL.svg")
plt.close()
plt.scatter(predictions_diff[:,1], input_label_test[:,1], s=1)
plt.title('diff_mean_ccell_3T3_undiff_CTRL'+str(stats.pearsonr(predictions_diff[:,0], input_label_test[:,1])))
plt.savefig(str(sys.argv[4])+"_mean_ccell_3T3_undiff_CTRL.svg")
plt.close()
plt.scatter(predictions_diff[:,2], input_label_test[:,2], s=1)
plt.title('diff_mean_cell_3T3_undiff_TGFB'+str(stats.pearsonr(predictions_diff[:,0], input_label_test[:,2])))
plt.savefig(str(sys.argv[4])+"_mean_cell_3T3_undiff_TGFB.svg")
plt.close()
plt.scatter(predictions_diff[:,3], input_label_test[:,3], s=1)
plt.title('diff_mean_RAW_CTRL'+str(stats.pearsonr(predictions_diff[:,0], input_label_test[:,3])))
plt.savefig(str(sys.argv[4])+"_mean_RAW_CTRL.svg")
plt.close()
plt.scatter(predictions_diff[:,4], input_label_test[:,4], s=1)
plt.title('diff_mean_RAW_IL1B'+str(stats.pearsonr(predictions_diff[:,0], input_label_test[:,4])))
plt.savefig(str(sys.argv[4])+"_mean_RAW_IL1B.svg")
plt.close()
plt.scatter(predictions_diff[:,5], input_label_test[:,5], s=1)
plt.title('diff_mean_cell_3T3_diff_CTRL'+str(stats.pearsonr(predictions_diff[:,0], input_label_test[:,5])))
plt.savefig(str(sys.argv[4])+"_mean_RAW_TGFB.svg")
plt.close()
plt.scatter(predictions_diff[:,6], input_label_test[:,6], s=1)
plt.title('diff_mean_TeloHAEC_CTRL'+str(stats.pearsonr(predictions_diff[:,0], input_label_test[:,6])))
plt.savefig(str(sys.argv[4])+"_mean_TeloHAEC_CTRL.svg")
plt.close()
plt.scatter(predictions_diff[:,7], input_label_test[:,7], s=1)
plt.title('diff_mean_TeloHAEC_IL1b_24h'+str(stats.pearsonr(predictions_diff[:,0], input_label_test[:,7])))
plt.savefig(str(sys.argv[4])+"_mean_TeloHAEC_IL1b_24h.svg")
plt.close()
plt.scatter(predictions_diff[:,8], input_label_test[:,8], s=1)
plt.title('diff_mean_TeloHAEC_IL1b_6h'+str(stats.pearsonr(predictions_diff[:,0], input_label_test[:,8])))
plt.savefig(str(sys.argv[4])+"_mean_TeloHAEC_IL1b_6h.svg")
plt.close()
plt.scatter(predictions_diff[:,9], input_label_test[:,9], s=1)
plt.title('diff_mean_HASMC_untreatedPilot'+str(stats.pearsonr(predictions_diff[:,0], input_label_test[:,9])))
plt.savefig(str(sys.argv[4])+"_mean_HASMC_untreatedPilot.svg")
plt.close()
plt.scatter(predictions_diff[:,10], input_label_test[:,10], s=1)
plt.title('diff_mean_HASMC_Chol'+str(stats.pearsonr(predictions_diff[:,0], input_label_test[:,10])))
plt.savefig(str(sys.argv[4])+"_mean_HASMC_Chol.svg")
plt.close()
plt.scatter(predictions_diff[:,11], input_label_test[:,11], s=1)
plt.title('diff_mean_HepG2_untreatedPilot'+str(stats.pearsonr(predictions_diff[:,0], input_label_test[:,11])))
plt.savefig(str(sys.argv[4])+"_mean_HepG2_untreatedPilot.svg")
plt.close()




















