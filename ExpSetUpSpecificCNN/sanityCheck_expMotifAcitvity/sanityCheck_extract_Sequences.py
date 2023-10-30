"""
Description:            This script takes STARRseq sequences as input and exports them in FASTA format required as input for FIMO analysis


Output:                  sequences in FASTA format

example bash command:   python ./ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/sanityCheck_extract_Sequences.py --seqFile starrseq-all-final-toorder_oligocomposition.csv --output test.fa

"""



import pandas as pd
import click
from tqdm import tqdm

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
    "--output",
    "out",
    required=True,
    multiple=False,
    type=str,
    help="e.g. test.fa",
)
def cli(seq_file, out):

    # import sequences
    df_IDs_Sequences=pd.read_csv(str(seq_file), sep=",",low_memory=False)

    # drop duplicates
    df_IDs_Sequences=df_IDs_Sequences.dropna()

    # save tmp
    df_IDs_Sequences.to_csv(str(seq_file)+".tmp")

    #import sequences
    seqs=open(str(seq_file)+".tmp", "r")
    seq_entries=seqs.readlines()
    seqs.close()

    # write to otfile
    outfile=open(str(out), "w")
    for i in tqdm(range (1, len(seq_entries), 1)):
        outfile.write(str(seq_entries[i].split(",")[1]) + "\n" + str(seq_entries[i].split(",")[4]) + "\n")
    outfile.close()

        
if __name__ == "__main__":
    cli()


