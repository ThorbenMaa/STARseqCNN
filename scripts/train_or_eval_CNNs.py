"""
Description:        This script has two modes. The "train" mode trains different CNNs on STARseq data provided by the Kaikkonen Lab using two different architectures and several different learning rates. 
                    Learning rates can be specified using the `learning_rate` array. The model is evaluated on a hold-out test data set. The "load" mode loads a pre-trained model and evaluates it on a
                    hold-out test data set. Pearson correlations between predicted and experimentally determined activities are calculated and scatter plots are generated.

Inputs:             Input 1 are labels for model training and evaluation. Input 2 are sequences for  model training and evaluation. Input 3 defines the mode (see above) and is either "train" or "load".
                    Input 4 is a pre-trained model that should be evaluated. Only important in the "load" mode of this script. Input 5 is the hold-out chromosme for testing. 

further parameters: can be specified in the parameters section of this script. 

Example commands:   `python train_or_eval_CNNs.py 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo1.txt starrseq-all-final-toorder_oligocomposition.csv train bla chr8`
                    `python train_or_eval_CNNs.py 2023-01-10_22-29-33\ myCounts.minDNAfilt.depthNorm.keepHaps\ -\ starr.haplotypes.oligo1.txt starrseq-all-final-toorder_oligocomposition.csv load CNN_StarSeq_model_Minna_deepSTAR_lr0.001 chr8`

Outputs:            Trained CNNs or evalutaion of a pre-trained CNN.

"""




import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
import os
os.environ["TF_CPP_MIN_LOG_LEVEL"]="2" #to supress warnings with loading tf
import tensorflow as tf
from tensorflow import keras
import numpy as np
import sys
from scipy import stats
import random
random.seed(1234)#to make training reproducable, 1234 is abitrary
np.random.seed(1234)#to make training reproducable, 1234 is abitrary
tf.random.set_seed(1234)#to make training reproducable, 1234 is abitrary

def one_hot_encode(seq): #taken from https://stackoverflow.com/questions/34263772/how-to-generate-one-hot-encoding-for-dna-sequences?rq=3
    mapping = dict(zip("ACGT", range(4)))    
    seq2 = [mapping[i] for i in seq]
    return np.eye(4)[seq2]


#parameters
sequence_length=198
learning_rate=[0.001, 0.0001, 0.00001]#0.001#0.01#0.1#0.001
batch_size=128
epochs=150
number_cell_types=12 #also hot encoded below, also needs to be change below when processing tarining and testing data

#import labels
df_IDs_reg_labels=pd.read_csv(sys.argv[1], sep="\t", decimal=',', low_memory=False)
df_IDs_reg_labels=df_IDs_reg_labels.drop_duplicates()

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

#split data to train and test data
df_IDs_seqs_reg_labels_test=df_IDs_seqs_reg_labels.loc[(df_IDs_seqs_reg_labels['ID'].str.contains(str(sys.argv[5])))==True]
df_IDs_seqs_reg_labels_train=df_IDs_seqs_reg_labels.loc[(df_IDs_seqs_reg_labels['ID'].str.contains(str(sys.argv[5])))==False]

#prepare data sets for model training and testing -> convert to tensors
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
input_label_test=tf.transpose(input_label_test)

#limit label to 400 to exclude outlyiers (arbitrary chosen, just looked at data plots)
input_label_test= tf.clip_by_value(input_label_test, clip_value_min=0, clip_value_max=400)

#sequences test data
input_seq_test=tf.convert_to_tensor(df_IDs_seqs_reg_labels_test["Seq one hot encoded"].to_list())

#label train data
input_label_train=tf.convert_to_tensor([df_IDs_seqs_reg_labels_train["mean_cell_3T3_diff_CTRL"].to_list(),
                                       df_IDs_seqs_reg_labels_train["mean_ccell_3T3_undiff_CTRL"].to_list(),
                                       df_IDs_seqs_reg_labels_train["mean_cell_3T3_undiff_TGFB"].to_list(),
                                       df_IDs_seqs_reg_labels_train["mean_RAW_CTRL"].to_list(),
                                       df_IDs_seqs_reg_labels_train["mean_RAW_IL1B"].to_list(),
                                       df_IDs_seqs_reg_labels_train["mean_RAW_TGFB"].to_list(),
                                       df_IDs_seqs_reg_labels_train["mean_TeloHAEC_CTRL"].to_list(),
                                       df_IDs_seqs_reg_labels_train["mean_TeloHAEC_IL1b_24h"].to_list(),
                                       df_IDs_seqs_reg_labels_train["mean_TeloHAEC_IL1b_6h"].to_list(),
                                       df_IDs_seqs_reg_labels_train["mean_HASMC_untreatedPilot"].to_list(),
                                       df_IDs_seqs_reg_labels_train["mean_HASMC_Chol"].to_list(),
                                       df_IDs_seqs_reg_labels_train["mean_HepG2_untreatedPilot"].to_list()
                                       
                                    
])
input_label_train=tf.transpose(input_label_train)

#limit label to 400 to exclude outlyiers (arbitrary chosen, just looked at data plots)
input_label_train= tf.clip_by_value(input_label_train, clip_value_min=0, clip_value_max=400)

#sequence train data
input_seq_train=tf.convert_to_tensor(df_IDs_seqs_reg_labels_train["Seq one hot encoded"].to_list())

#train or load model
#train
if sys.argv[3]=="train":

    #iterate over diferent learning rates specified above
    for lr in learning_rate:
        #standard archtecture of max's CNN pipeline (https://github.com/visze/sequence_cnn_models/tree/master)
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
            keras.layers.Dense(number_cell_types, activation='linear'),#predictions, regression, deshalb nur 1 output -> actibity of element

        ])

        print(model.summary())


        #training Max's standard model
        loss= keras.losses.MeanSquaredError()
        optim = tf.keras.optimizers.Adam(learning_rate=lr)
        metrics=["mse", "mae", "mape", "acc"]
        print("compile model")
        model.compile(loss=loss, optimizer=optim, metrics=metrics)
        random.seed(1234)#to make training reproducable, 42 is abitrary
        np.random.seed(1234)#to make training reproducable, 42 is abitrary
        tf.random.set_seed(1234)#to make training reproducable, 42 is abitrary
        model.fit(input_seq_train, input_label_train, batch_size=batch_size, epochs=epochs, shuffle=True, verbose=2)
        
        #save model
        print("save model"+str(lr))
        model.save("CNN_StarSeq_model_Minna_deepSTAR_lr0.001"+str(lr))
        
        #quick model evaluation on test data to select models for further evaluation; printed in output of this script
        model.evaluate(input_seq_test, input_label_test, batch_size=batch_size, verbose=2)

        #architecture more or less taken from deepstar (https://colab.research.google.com/drive/1Xgak40TuxWWLh5P5ARf0-4Xo0BcRn0Gd#scrollTo=wkq9fhH4wfoD)
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
            keras.layers.Dense(number_cell_types, activation='linear'),#predictions, regression, deshalb nur 1 output -> actibity of element

        ])

        print(model.summary())

        #training deepSTAR-like model
        loss= keras.losses.MeanSquaredError()
        optim = tf.keras.optimizers.Adam(learning_rate=lr)
        metrics=["mse", "mae", "mape", "acc"]
        print("compile model")
        model.compile(loss=loss, optimizer=optim, metrics=metrics)
        random.seed(1234)#to make training reproducable, 42 is abitrary
        np.random.seed(1234)#to make training reproducable, 42 is abitrary
        tf.random.set_seed(1234)#to make training reproducable, 42 is abitrary
        model.fit(input_seq_train, input_label_train, batch_size=batch_size, epochs=epochs, shuffle=True, verbose=2)
        
        #save model
        print("save model"+str(lr))
        model.save("CNN_StarSeq_model_Minna_deepSTAR_lr"+str(lr))

        #evaluate model on test data
        model.evaluate(input_seq_test, input_label_test, batch_size=batch_size, verbose=2)


#load model if "load" argument is given
elif sys.argv[3]=="load":
    model=keras.models.load_model(str(sys.argv[4]))

    # loaded model on test data
    model.evaluate(input_seq_test, input_label_test, batch_size=batch_size, verbose=2)

    #calculate predicted labels
    predictions=model.predict(input_seq_test, batch_size=batch_size, verbose=2)

    #correlations of predicted and experimental labels for different cell types
    print(stats.pearsonr(predictions[:,0], input_label_test[:,0]))
    print(stats.pearsonr(predictions[:,1], input_label_test[:,1]))
    print(stats.pearsonr(predictions[:,2], input_label_test[:,2]))
    print(stats.pearsonr(predictions[:,3], input_label_test[:,3]))
    print(stats.pearsonr(predictions[:,4], input_label_test[:,4]))
    print(stats.pearsonr(predictions[:,5], input_label_test[:,5]))
    print(stats.pearsonr(predictions[:,6], input_label_test[:,6]))
    print(stats.pearsonr(predictions[:,7], input_label_test[:,7]))
    print(stats.pearsonr(predictions[:,8], input_label_test[:,8]))
    print(stats.pearsonr(predictions[:,9], input_label_test[:,9]))
    print(stats.pearsonr(predictions[:,10], input_label_test[:,10]))
    print(stats.pearsonr(predictions[:,11], input_label_test[:,11]))

    #scatterplots of predicted and experimental labels for different cell types
    plt.scatter(predictions[:,0], input_label_test[:,0], s=1)
    plt.savefig(str(sys.argv[4])+"_mean_cell_3T3_diff_CTRL.svg")
    plt.close()
    plt.scatter(predictions[:,1], input_label_test[:,1], s=1)
    plt.savefig(str(sys.argv[4])+"_mean_ccell_3T3_undiff_CTRL.svg")
    plt.close()
    plt.scatter(predictions[:,2], input_label_test[:,2], s=1)
    plt.savefig(str(sys.argv[4])+"_mean_cell_3T3_undiff_TGFB.svg")
    plt.close()
    plt.scatter(predictions[:,3], input_label_test[:,3], s=1)
    plt.savefig(str(sys.argv[4])+"_mean_RAW_CTRL.svg")
    plt.close()
    plt.scatter(predictions[:,4], input_label_test[:,4], s=1)
    plt.savefig(str(sys.argv[4])+"_mean_RAW_IL1B.svg")
    plt.close()
    plt.scatter(predictions[:,5], input_label_test[:,5], s=1)
    plt.savefig(str(sys.argv[4])+"_mean_RAW_TGFB.svg")
    plt.close()
    plt.scatter(predictions[:,6], input_label_test[:,6], s=1)
    plt.savefig(str(sys.argv[4])+"_mean_TeloHAEC_CTRL.svg")
    plt.close()
    plt.scatter(predictions[:,7], input_label_test[:,7], s=1)
    plt.savefig(str(sys.argv[4])+"_mean_TeloHAEC_IL1b_24h.svg")
    plt.close()
    plt.scatter(predictions[:,8], input_label_test[:,8], s=1)
    plt.savefig(str(sys.argv[4])+"_mean_TeloHAEC_IL1b_6h.svg")
    plt.close()
    plt.scatter(predictions[:,9], input_label_test[:,9], s=1)
    plt.savefig(str(sys.argv[4])+"_mean_HASMC_untreatedPilot.svg")
    plt.close()
    plt.scatter(predictions[:,10], input_label_test[:,10], s=1)
    plt.savefig(str(sys.argv[4])+"_mean_HASMC_Chol.svg")
    plt.close()
    plt.scatter(predictions[:,11], input_label_test[:,11], s=1)
    plt.savefig(str(sys.argv[4])+"_mean_HepG2_untreatedPilot.svg")
    plt.close()