# STARseqCNN
This repository contains and describes code used to train and evaluate multi task STARseq CNNs trained on STARseq data provided by the Kaikkonen lab.

## Worklflow

### Install dependencies
I recommend to use mamba to create environments and install dependencies:

```
mamba env create --name CNN_TM --file=./envs/CNN_TM.yml
mamba env create --name modisco_lite --file=./envs/modisco_lite.yml
```

### load training labels and sequences (will be added upon publication)

### train CNNs
First, activate the `CNN_TM` environment using `mamba activate CNN_TM`.
You can either directly run `train_or_eval_CNNs.py` in the "train" mode (documentation and example bash commands provided in the script) or run the script using slurm within an sbatch script. An example with reasonable recources is given in `sbatch_scripts/sbtachTrain_CNN_TM.sh` (example bash command given in script). 

The script prints a lot of information. Among others, after each training, it will print the perfromance on a hold-out test data set. The model with the best perfromance can than be used for further analysis.
In the manuscript, the DeepSTAR-like architecture with a learning rate of 0.001 was chosen. 

### further evaluate CNNs
First, activate the `CNN_TM` environment using `mamba activate CNN_TM`.
This can be done using the `train_or_eval_CNNs.py` script in the "load" mode (documentation and example bash commands provided in the script). The script will calculate pearson correlations between 
predicted and experimental STARseq activity and plot corresponding scatter plots.

### In-silico-mutagenisis (ISM)
First, activate the `CNN_TM` environment using `mamba activate CNN_TM`.
You can use the `ism_TM.py` script to generate one-hot-encoded sequences and the corresponding raw CNN scores in a format suitable for tfmodisco-lite for all experimental set-ups explained in the manuscript. Note that the CNN raw scores can not be directly used
for tfmodisco-lite but have to be further processed using `modisco_TM.py` (documentation and example bash commands provided in the respective scripts). To make use of slurm, you can also use an sbatch script provided in `./sbatch_scripts/sbatch_ism.sh`. 

### tfmodisco-lizte analysis
First, activate the `modisco_lite` environment using `mamba activate modisco_lite`.
Use the script `mosidco_TM.py` with the .npz files generated in the previous step by the `ism_TM.py` script as inputs (documentation and example bash commands provided in the respective scripts). The script will create a tfmodisco-lite result.h5 file. 
You can create tfmodisco-lite reports from this results file using e.g.:
```
wget https://jaspar.genereg.net/download/data/2022/CORE/JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
modisco report -i modisco_resultshypothetical_contribution_scores_mean_HASMC_Chol.npz.h5 -o report_HASMC_chol/ -s report_HASMC_Chol/ -m JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme.txt
```
You can also use an sbatch script provided in `./sbatch_scripts/sbatch_tfmodisco.sh` to make use of slurm. Here, all commands for generating tfmodisco results and reports used for the analysis of the manuscript are listed. 

More information on how to use tfmodisco-lite are given at the corresponding github repository https://github.com/jmschrei/tfmodisco-lite/tree/main. 

### Sanity check
Has the CNN really learned motifs that enhance/repress activity in the STARseq experiment? First, activate the `CNN_TM` environment using `mamba activate CNN_TM`. You can use the `sanity_check_modisco_results.py` script to plot experimental activity of sequences containing a motif of interest or not (documentation and example bash command provided in the script). The result will look like this:
![alt text for screen readers](boxplot_HASMC_CholTGAGTCA.svg "Boxplots")

