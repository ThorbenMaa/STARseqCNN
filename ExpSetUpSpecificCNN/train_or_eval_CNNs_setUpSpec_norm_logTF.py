"""
Description:        This script has two modes. The "train" mode trains different CNNs on STARseq data provided by the Kaikkonen Lab using two different architectures and several different learning rates, with and without augmentation. 
                    Learning rates can be specified using the `learning_rate` array. The model is evaluated on a hold-out test data set. The "load" mode loads a pre-trained model and evaluates it on a
                    hold-out test data set. Pearson correlations between predicted and experimentally determined activities are calculated and scatter plots are generated.
                    The aim is training a difference CNN that learns to predict LOG TRANSFORMED AND STANDARDIZED differences in STARRseq activites between 2 experimental set ups

Inputs:             see click

further parameters: can be specified in the parameters section of this script. 


Outputs:            Trained diffCNNs or evalutaion of a pre-trained diffCNN.

Command:            python ExpSetUpSpecificCNN/train_or_eval_CNNs_setUpSpec.py --counts1 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo1.txt  --counts2 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo2.txt  --seqs starrseq-all-final-toorder_oligocomposition.csv  --mode load  --holdOut chr8  --useAug use_aug  --model allseq-CNN_StarSeq_model_Minna_deepSTAR_lr0.01use_augsetUpSpec\('mean_HepG2_untreatedPilot',\ 'mean_TeloHAEC_CTRL'\)


"""

#set seed and load dependencies
seed_value=1234
import os
from types import CellType
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
import click

@click.command()
@click.option(
    "--counts1",
    "counts1_file",
    required=True,
    multiple=False,
    type=str,
    help="first counts file",
)
@click.option(
    "--counts2",
    "counts2_file",
    required=True,
    multiple=False,
    type=str,
    help="second counts file",
)
@click.option(
    "--seqs",
    "seqs_file",
    required=True,
    multiple=False,
    type=str,
    help="sequences file",
)
@click.option(
    "--mode",
    "mode",
    required=True,
    multiple=False,
    type=str,
    help="training (= train) vs evaluation (= load) mode",
)
@click.option(
    "--model",
    "pretrained_model",
    required=True,
    multiple=False,
    type=str,
    default="None",
    help="pretrained model",
)
@click.option(
    "--holdOut",
    "hold_out_chrom",
    required=True,
    multiple=False,
    type=str,
    help="hold out chromosome for testing",
)
@click.option(
    "--useAug",
    "use_aug",
    required=True,
    multiple=False,
    type=str,
    help="use augmentation for training. = use_aug if yes and anything if not",
)
@click.option(
    "--compare",
    "setUps_to_compare",
    required=True,
    multiple=True,
    type=str,
    #default=["mean_TeloHAEC_CTRL", "mean_TeloHAEC_IL1b_6h"],# ["mean_HepG2_untreatedPilot", "mean_TeloHAEC_CTRL"]#
    help="two set ups to be copmpared (oprtzions : ...)",
)
def cli(counts1_file, counts2_file, seqs_file, mode, hold_out_chrom, use_aug, pretrained_model, setUps_to_compare):
    #parameters
    sequence_length=198
    learning_rate=[0.01]#, 0.001, 0.0001]
    batch_size=128
    epochs=100
    number_cell_types=2 #also hot encoded below, also needs to be change below when processing tarining and testing data

    #import labels
    df_IDs_reg_labels=pd.read_csv(counts1_file, sep="\t", decimal=',', low_memory=False)
    df_IDs_reg_labels2=pd.read_csv(counts2_file, sep="\t", decimal=',', low_memory=False)

    df_IDs_reg_labels = pd.concat([df_IDs_reg_labels, df_IDs_reg_labels2], axis=0)

    df_IDs_reg_labels=df_IDs_reg_labels.drop_duplicates()

    #replace 0 with 1 for log transformation later
    #df_IDs_reg_labels=df_IDs_reg_labels.replace(0, 1)

    #import sequences
    df_IDs_Sequences=pd.read_csv(seqs_file, sep=",",low_memory=False)

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

    #average data over replica to generate labels, systematic +1 to be able to normalize by dividing by mean_input2022Dec

    #input for normalization
    df_IDs_seqs_reg_labels['mean_input2022Dec'] = df_IDs_seqs_reg_labels.loc[:, ["input2022Dec_50ng_rep1_2022_12_14", "input2022Dec_50ng_rep2_2022_12_14", "input2022Dec_50ng_rep3_2022_12_14", "input2022Dec_50ng_rep4_2022_12_14"]].mean(axis=1) + 1

    #3T3
    df_IDs_seqs_reg_labels['mean_cell_3T3_diff_CTRL'] = df_IDs_seqs_reg_labels.loc[:, ["cell_3T3_diff_CTRL_rep1_2022_12_14", "cell_3T3_diff_CTRL_rep2_2022_12_14", "cell_3T3_diff_CTRL_rep3_2022_12_14", "cell_3T3_diff_CTRL_rep4_2022_12_14"]].mean(axis=1) + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    df_IDs_seqs_reg_labels['mean_ccell_3T3_undiff_CTRL'] = df_IDs_seqs_reg_labels.loc[:, ["cell_3T3_undiff_CTRL_rep1_2022_12_14", "cell_3T3_undiff_CTRL_rep2_2022_12_14", "cell_3T3_undiff_CTRL_rep3_2022_12_14", "cell_3T3_undiff_CTRL_rep4_2022_12_14"]].mean(axis=1)  + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    df_IDs_seqs_reg_labels['mean_cell_3T3_undiff_TGFB'] = df_IDs_seqs_reg_labels.loc[:, ["cell_3T3_undiff_TGFB_rep1_2022_12_14", "cell_3T3_undiff_TGFB_rep2_2022_12_14", "cell_3T3_undiff_TGFB_rep3_2022_12_14", "cell_3T3_undiff_TGFB_rep4_2022_12_14"]].mean(axis=1)  + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    #RAW
    df_IDs_seqs_reg_labels['mean_RAW_CTRL'] = df_IDs_seqs_reg_labels.loc[:, ["RAW_CTRL_rep1_2022_12_14", "RAW_CTRL_rep2_2022_12_14", "RAW_CTRL_rep3_2022_12_14", "RAW_CTRL_rep4_2022_12_14"]].mean(axis=1)  + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    df_IDs_seqs_reg_labels['mean_RAW_IL1B'] = df_IDs_seqs_reg_labels.loc[:, ["RAW_IL1B_rep1_2022_12_14", "RAW_IL1B_rep2_2022_12_14", "RAW_IL1B_rep3_2022_12_14", "RAW_IL1B_rep4_2022_12_14"]].mean(axis=1)  + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    df_IDs_seqs_reg_labels['mean_RAW_TGFB'] = df_IDs_seqs_reg_labels.loc[:, ["RAW_TGFB_rep1_2022_12_14", "RAW_TGFB_rep2_2022_12_14", "RAW_TGFB_rep3_2022_12_14", "RAW_TGFB_rep4_2022_12_14"]].mean(axis=1)  + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    #TeloHAEC
    df_IDs_seqs_reg_labels['mean_TeloHAEC_CTRL'] = df_IDs_seqs_reg_labels.loc[:, ["TeloHAEC_CTRL_rep1", "TeloHAEC_CTRL_rep2", "TeloHAEC_CTRL_rep3"]].mean(axis=1)  + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    df_IDs_seqs_reg_labels['mean_TeloHAEC_IL1b_24h'] = df_IDs_seqs_reg_labels.loc[:, ["TeloHAEC_IL1b_24h_rep1", "TeloHAEC_IL1b_24h_rep2", "TeloHAEC_IL1b_24h_rep3"]].mean(axis=1)  + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    df_IDs_seqs_reg_labels['mean_TeloHAEC_IL1b_6h'] = df_IDs_seqs_reg_labels.loc[:, ["TeloHAEC_IL1b_6h_rep1", "TeloHAEC_IL1b_6h_rep2", "TeloHAEC_IL1b_6h_rep3"]].mean(axis=1)  + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    #HASMC
    df_IDs_seqs_reg_labels['mean_HASMC_untreatedPilot'] = df_IDs_seqs_reg_labels.loc[:, ["HASMC_untreatedPilot_rep1", "HASMC_untreatedPilot_rep2", "HASMC_untreatedPilot_rep3"]].mean(axis=1)  + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    df_IDs_seqs_reg_labels['mean_HASMC_Chol'] = df_IDs_seqs_reg_labels.loc[:, ["HASMC_Chol_rep1", "HASMC_Chol_rep2", "HASMC_Chol_rep3"]].mean(axis=1) + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']
    #HepG2
    df_IDs_seqs_reg_labels['mean_HepG2_untreatedPilot'] = df_IDs_seqs_reg_labels.loc[:, ["HepG2_untreatedPilot_rep1", "HepG2_untreatedPilot_rep2", "HepG2_untreatedPilot_rep3"]].mean(axis=1)  + 1 / df_IDs_seqs_reg_labels['mean_input2022Dec']




    #split data to train and test data
    df_IDs_seqs_reg_labels_test=df_IDs_seqs_reg_labels.loc[(df_IDs_seqs_reg_labels['ID'].str.contains(str(hold_out_chrom)))==True]
    df_IDs_seqs_reg_labels_train=df_IDs_seqs_reg_labels.loc[(df_IDs_seqs_reg_labels['ID'].str.contains(str(hold_out_chrom)))==False]

    #prepare data sets for model training and testing -> convert to tensors + log transform of labels and data type to int 8 of sequences; with or without augmentation
    if use_aug=="use_aug":
        #add augmentation on training data
        #covert seqs to complementary reverse seqs and add as new column
        df_IDs_seqs_reg_labels_train['compSeq']=df_IDs_seqs_reg_labels_train['enhancer'].apply(complementary)

        #convert to list and create one-hot-encoded
        sequence_list=df_IDs_seqs_reg_labels_train['compSeq'].to_list()
        sequences_tok=[]

        for i in range (0, len(sequence_list), 1): 
            sequences_tok.append(one_hot_encode(sequence_list[i]))
        
        #add one-hot-encoded to data frame
        df_IDs_seqs_reg_labels_train['compSeq one hot encoded'] = sequences_tok

        #label train data (all values 2 times because of augmentation)
        #cell_types=["mean_cell_3T3_diff_CTRL", "mean_ccell_3T3_undiff_CTRL", "mean_cell_3T3_undiff_TGFB", "mean_RAW_CTRL", "mean_RAW_IL1B", "mean_RAW_TGFB", "mean_TeloHAEC_CTRL", "mean_TeloHAEC_IL1b_24h", "mean_TeloHAEC_IL1b_6h", "mean_HASMC_untreatedPilot", "mean_HASMC_Chol", "mean_HepG2_untreatedPilot"]
        setUps_to_compare_list=[]
        
        # standard normalization (instead of log transformation)
        df_IDs_seqs_reg_labels_train["Set up Difference"] = (
                                                             ((np.log(df_IDs_seqs_reg_labels_train[setUps_to_compare[0]]) - np.log(df_IDs_seqs_reg_labels_train[setUps_to_compare[0]].mean())) / np.log(df_IDs_seqs_reg_labels_train[setUps_to_compare[0]].std()) ) 
                                                            -((np.log(df_IDs_seqs_reg_labels_train[setUps_to_compare[1]]) - np.log(df_IDs_seqs_reg_labels_train[setUps_to_compare[1]].mean())) / np.log(df_IDs_seqs_reg_labels_train[setUps_to_compare[1]].std()) )
                                                            )
                                                            
                                                            
        #print(df_IDs_seqs_reg_labels_train["Set up Difference"]) 
        setUps_to_compare_list.append(df_IDs_seqs_reg_labels_train["Set up Difference"].to_list() + df_IDs_seqs_reg_labels_train["Set up Difference"].to_list())
        
        input_label_train=tf.convert_to_tensor(setUps_to_compare_list)

        #transpose
        input_label_train=tf.transpose(input_label_train)
        #print(input_label_train)

        #sequence train data (original sequence + augmented by using complemantary reverse strand sequence)
        input_seq_train=tf.cast(tf.convert_to_tensor(df_IDs_seqs_reg_labels_train["Seq one hot encoded"].to_list() + df_IDs_seqs_reg_labels_train['compSeq one hot encoded'].to_list()), tf.int8)

    else: #i.e. no augmentation
        
        #label train data 
        setUps_to_compare_list=[]
        # log fold 
        df_IDs_seqs_reg_labels_train["Set up Difference"] = (
                                                             ((np.log(df_IDs_seqs_reg_labels_train[setUps_to_compare[0]]) - np.log(df_IDs_seqs_reg_labels_train[setUps_to_compare[0]].mean())) / np.log(df_IDs_seqs_reg_labels_train[setUps_to_compare[0]].std()) ) 
                                                            -((np.log(df_IDs_seqs_reg_labels_train[setUps_to_compare[1]]) - np.log(df_IDs_seqs_reg_labels_train[setUps_to_compare[1]].mean())) / np.log(df_IDs_seqs_reg_labels_train[setUps_to_compare[1]].std()) )
                                                            )
        #print(df_IDs_seqs_reg_labels_train["Set up Difference"]) 
        print ("traiing data mean", df_IDs_seqs_reg_labels_train[setUps_to_compare[1]].mean())
        setUps_to_compare_list.append(df_IDs_seqs_reg_labels_train["Set up Difference"].to_list())
        

        input_label_train=tf.convert_to_tensor(setUps_to_compare_list)

        #transpose
        input_label_train=tf.transpose(input_label_train)
        #print(input_label_train)

        #sequence train data (original sequence + augmented by using complemantary strand sequence)
        input_seq_train=tf.cast(tf.convert_to_tensor(df_IDs_seqs_reg_labels_train["Seq one hot encoded"].to_list()), tf.int8)

    #labels test data
    setUps_to_compare_list=[]
    # log fold 
    df_IDs_seqs_reg_labels_test["Set up Difference"] = (
                                                        ((np.log(df_IDs_seqs_reg_labels_test[setUps_to_compare[0]]) - np.log(df_IDs_seqs_reg_labels_test[setUps_to_compare[0]].mean())) / np.log(df_IDs_seqs_reg_labels_test[setUps_to_compare[0]].std()) ) 
                                                        -((np.log(df_IDs_seqs_reg_labels_test[setUps_to_compare[1]]) - np.log(df_IDs_seqs_reg_labels_test[setUps_to_compare[1]].mean())) / np.log(df_IDs_seqs_reg_labels_test[setUps_to_compare[1]].std()) )
                                                        )
    print(df_IDs_seqs_reg_labels_test["Set up Difference"])

    setUps_to_compare_list.append(df_IDs_seqs_reg_labels_test["Set up Difference"].to_list())

    input_label_test=tf.convert_to_tensor(setUps_to_compare_list)

    #transpose
    input_label_test=tf.transpose(input_label_test)
    #print(input_label_test)

    #sequences test data
    input_seq_test=tf.cast(tf.convert_to_tensor(df_IDs_seqs_reg_labels_test["Seq one hot encoded"].to_list()), tf.int8)


    #train or load model
    #train
    if mode=="train":
        
        #iterate over diferent learning rates specified above
        for lr in learning_rate:
            """
            #standard archtecture of max's CNN pipeline (https://github.com/visze/sequence_cnn_models/tree/master)
            random.seed(1234)#to make training reproducable,
            np.random.seed(1234)#to make training reproducable, 
            tf.random.set_seed(1234)#to make training reproducable,
            model = keras.Sequential([
                keras.layers.InputLayer(input_shape=(sequence_length,4), name="input"),
                keras.layers.Conv1D(250, kernel_size=7, strides=1, activation='relu', name="conv1"),  # 250 7 relu
                keras.layers.BatchNormalization(), 
                keras.layers.Conv1D(250, 8, strides=1, activation='softmax', name="conv2"),
                keras.layers.MaxPooling1D(pool_size=2, strides=None, name="maxpool1"),
                keras.layers.Dropout(0.1),
                keras.layers.Conv1D(250, 3, strides=1, activation='softmax', name="conv3"),
                keras.layers.BatchNormalization(),
                keras.layers.Conv1D(100, 2, strides=1, activation='softmax', name="conv4"),
                keras.layers.BatchNormalization(),
                keras.layers.MaxPooling1D(pool_size=1, strides=None, name="maxpool2"),
                keras.layers.Dropout(0.1),
                keras.layers.Flatten(),
                keras.layers.Dense(300, activation='sigmoid'),
                keras.layers.Dropout(0.3),
                keras.layers.Dense(200, activation='sigmoid'),
                keras.layers.Dense(1, activation='linear'),#predictions, regression, deshalb nur 1 output -> actibity of element
            ])

            print(model.summary())

            #training Max's standard model
            loss= keras.losses.MeanSquaredError()
            optim = tf.keras.optimizers.Adam(learning_rate=lr)
            metrics=["mse", "mae", "mape", "acc"]
            print("compile model")
            model.compile(loss=loss, optimizer=optim, metrics=metrics)

            model.fit(input_seq_train, input_label_train, batch_size=batch_size, epochs=epochs, shuffle=True, verbose=2)
            
            #save model
            print("save model"+str(lr))
            model.save("all_seq-CNN_StarSeq_model_Minna_Max"+str(lr)+str(use_aug)+"setUpSpec"+str(setUps_to_compare))
            
            #quick model evaluation on test data to select models for further evaluation; printed in output of this script
            model.evaluate(input_seq_test, input_label_test, batch_size=batch_size, verbose=2)
            """
        
            #architecture more or less taken from deepstar (https://colab.research.google.com/drive/1Xgak40TuxWWLh5P5ARf0-4Xo0BcRn0Gd#scrollTo=wkq9fhH4wfoD)
            random.seed(seed_value)#to make training reproducable,
            np.random.seed(seed_value)#to make training reproducable,
            tf.random.set_seed(seed_value)#to make training reproducable,
            model = keras.Sequential([
                keras.layers.InputLayer(input_shape=(sequence_length,4), name="input"),
                keras.layers.Conv1D(128, kernel_size=7,  activation='relu', padding='same', name="conv1"),  # 250 7 relu
                keras.layers.BatchNormalization(), 
                keras.layers.MaxPooling1D(2, name="maxpool1"),
                keras.layers.Conv1D(60, kernel_size=3,  activation='relu', padding='same', name="conv2"),  # 250 7 relu
                keras.layers.BatchNormalization(), 
                keras.layers.MaxPooling1D(2, name="maxpool2"),
                keras.layers.Dropout(0.4),
                keras.layers.Conv1D(60, kernel_size=5,  activation='relu', padding='same', name="conv3"),  # 250 7 relu
                keras.layers.BatchNormalization(), 
                keras.layers.MaxPooling1D(pool_size=2,  name="maxpool3"),
                keras.layers.Dropout(0.4),
                keras.layers.Flatten(),
                keras.layers.Dense(64, activation='relu'),
                keras.layers.BatchNormalization(), 
                keras.layers.Dropout(0.4),
                keras.layers.Dense(1, activation='linear'),#predictions, regression, deshalb nur 1 output -> actibity of element
            ])

            print(model.summary())

            #training deepSTAR-like model
            loss= keras.losses.MeanSquaredError()
            optim = tf.keras.optimizers.Adam(learning_rate=lr)
            metrics=["mse", "mae", "mape", "acc"]
            print("compile model")
            model.compile(loss=loss, optimizer=optim, metrics=metrics)
            model.fit(input_seq_train, input_label_train, batch_size=batch_size, epochs=epochs, shuffle=True, verbose=2)
            
            #save model
            print("save model"+str(lr))
            model.save("allseq-CNN_StarSeq_model_Minna_deepSTAR_lr"+str(lr)+str(use_aug)+"setUpSpecNorm_lgTF_"+str(setUps_to_compare))

            #evaluate model on test data
            model.evaluate(input_seq_test, input_label_test, batch_size=batch_size, verbose=2)


    #load model if "load" argument is given
    elif mode=="load":
        print("start evalualtion")
        model=keras.models.load_model(str(pretrained_model))

        # loaded model on test data
        print ("Evaluation:")
        model.evaluate(input_seq_test, input_label_test, batch_size=batch_size, verbose=2)

        #calculate predicted labels
        print ("calculate predictions...")
        predictions=model.predict(input_seq_test, batch_size=batch_size, verbose=2)

        #correlations of predicted and experimental labels for different cell types
        print ("calculate correlations...")
        print(stats.pearsonr(predictions[:,0], input_label_test[:,0]))


        #scatterplots of predicted and experimental labels for different cell types
        plt.scatter(predictions[:,0], input_label_test[:,0], s=1)
        plt.title(str(setUps_to_compare[0])+str(stats.pearsonr(predictions[:,0], input_label_test[:,0])))
        plt.savefig(str(pretrained_model)+str(setUps_to_compare[0])+"_vs_"+str(setUps_to_compare[0])+".svg")
        plt.close()


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




if __name__ == "__main__":
    cli()



