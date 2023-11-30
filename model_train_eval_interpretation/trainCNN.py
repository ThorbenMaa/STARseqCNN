#set seed and load dependencies
seed_value=1234
import os
from xmlrpc.client import boolean
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
    "--augmentation",
    "augmentation",
    required=True,
    multiple=False,
    type=boolean,
    default=False,
    help="use augmentation",
)
@click.option(
    "--earlyStop",
    "earlyStop",
    required=True,
    multiple=False,
    type=boolean,
    default=False,
    help="use early stopping",
)
@click.option(
    "--learningRates",
    "learning_rates",
    required=True,
    multiple=True,
    type=float,
    default=[0.01],
    help="learning rates",
)
@click.option(
    "--epochs",
    "epochs",
    required=True,
    multiple=False,
    type=int,
    default=150,
    help="epochs",
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
    "--n_tasks",
    "n_tasks",
    required=True,
    multiple=False,
    type=int,
    default=12,
    help="output channels of model",
)
@click.option(
    "--modelName",
    "model_name",
    required=True,
    multiple=False,
    type=str,
    default="models/TM_s_model",
    help="name uder that the model is stored",
)
@click.option(
    "--holdOutChr",
    "hold_out_chr",
    required=True,
    multiple=False,
    type=str,
    default="chr8",
    help="e.g., chr8 . used as hold out for testing. cross validation is not advisible if you don't know how the seqs were designed",
)
def cli(activity_file, seq_file, augmentation, learning_rates, epochs, sequence_length, batch_size, n_tasks, model_name, hold_out_chr, earlyStop):

    #import labels
    df_list=[]
    for i in range (0, len(activity_file), 1):
        df_list.append(pd.read_csv(activity_file[i], sep="\t", decimal=',', low_memory=False))
    df_IDs_reg_labels = pd.concat(df_list, axis=0)
    df_IDs_reg_labels=df_IDs_reg_labels.drop_duplicates()


    #import sequences
    df_IDs_Sequences=pd.read_csv(seq_file, sep=",",low_memory=False)

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



    #split data to train and test data
    df_IDs_seqs_reg_labels_test=df_IDs_seqs_reg_labels.loc[(df_IDs_seqs_reg_labels['ID'].str.contains(str(hold_out_chr)))==True]
    df_IDs_seqs_reg_labels_train=df_IDs_seqs_reg_labels.loc[(df_IDs_seqs_reg_labels['ID'].str.contains(str(hold_out_chr)))==False]

    #prepare data sets for model training and testing -> convert to tensors + log transform of labels and data type to int 8 of sequences; with or without augmentation
    if augmentation==True:
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
        input_label_train=tf.convert_to_tensor([df_IDs_seqs_reg_labels_train["mean_cell_3T3_diff_CTRL"].to_list() + df_IDs_seqs_reg_labels_train["mean_cell_3T3_diff_CTRL"].to_list(),
                                            df_IDs_seqs_reg_labels_train["mean_ccell_3T3_undiff_CTRL"].to_list() + df_IDs_seqs_reg_labels_train["mean_ccell_3T3_undiff_CTRL"].to_list(),
                                            df_IDs_seqs_reg_labels_train["mean_cell_3T3_undiff_TGFB"].to_list() + df_IDs_seqs_reg_labels_train["mean_cell_3T3_undiff_TGFB"].to_list(),
                                            df_IDs_seqs_reg_labels_train["mean_RAW_CTRL"].to_list() + df_IDs_seqs_reg_labels_train["mean_RAW_CTRL"].to_list(),
                                            df_IDs_seqs_reg_labels_train["mean_RAW_IL1B"].to_list() + df_IDs_seqs_reg_labels_train["mean_RAW_IL1B"].to_list(),
                                            df_IDs_seqs_reg_labels_train["mean_RAW_TGFB"].to_list() + df_IDs_seqs_reg_labels_train["mean_RAW_TGFB"].to_list(),
                                            df_IDs_seqs_reg_labels_train["mean_TeloHAEC_CTRL"].to_list() + df_IDs_seqs_reg_labels_train["mean_TeloHAEC_CTRL"].to_list(),
                                            df_IDs_seqs_reg_labels_train["mean_TeloHAEC_IL1b_24h"].to_list() + df_IDs_seqs_reg_labels_train["mean_TeloHAEC_IL1b_24h"].to_list(),
                                            df_IDs_seqs_reg_labels_train["mean_TeloHAEC_IL1b_6h"].to_list() + df_IDs_seqs_reg_labels_train["mean_TeloHAEC_IL1b_6h"].to_list(),
                                            df_IDs_seqs_reg_labels_train["mean_HASMC_untreatedPilot"].to_list() + df_IDs_seqs_reg_labels_train["mean_HASMC_untreatedPilot"].to_list(),
                                            df_IDs_seqs_reg_labels_train["mean_HASMC_Chol"].to_list() + df_IDs_seqs_reg_labels_train["mean_HASMC_Chol"].to_list(),
                                            df_IDs_seqs_reg_labels_train["mean_HepG2_untreatedPilot"].to_list() + df_IDs_seqs_reg_labels_train["mean_HepG2_untreatedPilot"].to_list()                                   
        ])

        #log transform adn transpose
        input_label_train=tf.math.log(tf.transpose(input_label_train))

        #sequence train data (original sequence + augmented by using complemantary reverse strand sequence)
        input_seq_train=tf.cast(tf.convert_to_tensor(df_IDs_seqs_reg_labels_train["Seq one hot encoded"].to_list() + df_IDs_seqs_reg_labels_train['compSeq one hot encoded'].to_list()), tf.int8)

    else: #i.e. no augmentation
        
        #label train data 
        input_label_train=tf.convert_to_tensor([df_IDs_seqs_reg_labels_train["mean_cell_3T3_diff_CTRL"].to_list() ,
                                            df_IDs_seqs_reg_labels_train["mean_ccell_3T3_undiff_CTRL"].to_list() ,
                                            df_IDs_seqs_reg_labels_train["mean_cell_3T3_undiff_TGFB"].to_list() ,
                                            df_IDs_seqs_reg_labels_train["mean_RAW_CTRL"].to_list(),
                                            df_IDs_seqs_reg_labels_train["mean_RAW_IL1B"].to_list(),
                                            df_IDs_seqs_reg_labels_train["mean_RAW_TGFB"].to_list(),
                                            df_IDs_seqs_reg_labels_train["mean_TeloHAEC_CTRL"].to_list() ,
                                            df_IDs_seqs_reg_labels_train["mean_TeloHAEC_IL1b_24h"].to_list() ,
                                            df_IDs_seqs_reg_labels_train["mean_TeloHAEC_IL1b_6h"].to_list(),
                                            df_IDs_seqs_reg_labels_train["mean_HASMC_untreatedPilot"].to_list() ,
                                            df_IDs_seqs_reg_labels_train["mean_HASMC_Chol"].to_list(),
                                            df_IDs_seqs_reg_labels_train["mean_HepG2_untreatedPilot"].to_list()                                   
        ])

        #log transform adn transpose
        input_label_train=tf.math.log(tf.transpose(input_label_train))

        #sequence train data (original sequence + augmented by using complemantary strand sequence)
        input_seq_train=tf.cast(tf.convert_to_tensor(df_IDs_seqs_reg_labels_train["Seq one hot encoded"].to_list()), tf.int8)

    #labels test data
    input_label_test=tf.convert_to_tensor([df_IDs_seqs_reg_labels_test["mean_cell_3T3_diff_CTRL"].to_list(),
                                        df_IDs_seqs_reg_labels_test["mean_ccell_3T3_undiff_CTRL"].to_list(),
                                        df_IDs_seqs_reg_labels_test["mean_cell_3T3_undiff_TGFB"].to_list(),
                                        df_IDs_seqs_reg_labels_test["mean_RAW_CTRL"].to_list(),
                                        df_IDs_seqs_reg_labels_test["mean_RAW_IL1B"].to_list(),
                                        df_IDs_seqs_reg_labels_test["mean_RAW_TGFB"].to_list(),
                                        df_IDs_seqs_reg_labels_test["mean_TeloHAEC_CTRL"].to_list(),
                                        df_IDs_seqs_reg_labels_test["mean_TeloHAEC_IL1b_24h"].to_list(),
                                        df_IDs_seqs_reg_labels_test["mean_TeloHAEC_IL1b_6h"].to_list(),
                                        df_IDs_seqs_reg_labels_test["mean_HASMC_untreatedPilot"].to_list(),
                                        df_IDs_seqs_reg_labels_test["mean_HASMC_Chol"].to_list(),
                                        df_IDs_seqs_reg_labels_test["mean_HepG2_untreatedPilot"].to_list()                                    
    ])

    #log transform and transpose
    input_label_test=tf.math.log(tf.transpose(input_label_test))
    #sequences test data
    input_seq_test=tf.cast(tf.convert_to_tensor(df_IDs_seqs_reg_labels_test["Seq one hot encoded"].to_list()), tf.int8)


    #train model
    #iterate over diferent learning rates specified above
    for lr in learning_rates:
        
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
            keras.layers.Dense(n_tasks, activation='linear'),#predictions, regression, deshalb nur 1 output -> actibity of element
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
        model.save("all_seq-CNN_StarSeq_model_Minna_Max"+str(lr)+str(sys.argv[6]))
        
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
            keras.layers.Dense(n_tasks, activation='linear'),#predictions, regression, deshalb nur 1 output -> actibity of element
        ])

        print(model.summary())

        #training deepSTAR-like model
        loss= keras.losses.MeanSquaredError()
        optim = tf.keras.optimizers.Adam(learning_rate=lr)
        metrics=["mse", "mae", "mape", "acc"]
        print("compile model")
        model.compile(loss=loss, optimizer=optim, metrics=metrics)
        if earlyStop==False:
            model.fit(input_seq_train, input_label_train, batch_size=batch_size, epochs=epochs, shuffle=True, verbose=2)
        else:
            print("use early stopping...")
            callback = keras.callbacks.EarlyStopping(monitor='loss', patience=3, restore_best_weights=True) #fiiting stopps if after 3 eopchs no improvement is made
            history = model.fit(input_seq_train, input_label_train, validation_split=0.05, batch_size=batch_size, epochs=epochs, shuffle=True, verbose=2, callbacks=[callback])
        
        #save model
        print("save model"+str(lr))
        model.save(str(model_name))

        #evaluate model on test data
        model.evaluate(input_seq_test, input_label_test, batch_size=batch_size, verbose=2)


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