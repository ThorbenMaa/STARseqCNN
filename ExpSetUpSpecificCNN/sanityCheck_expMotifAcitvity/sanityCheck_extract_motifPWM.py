"""
Description:            This script extracts PWMs from e.g. JASPAR data base to pass them to MEME FIMO in a subsequent step

Output:                 Motif PWMs of interest from JASPAR

example bash command:   python ./ExpSetUpSpecificCNN/sanityCheck_expMotifAcitvity/sanityCheck_extract_motifPWM.py --motif MA0029.1_Mecom --motif MA0488.1_JUN --output test.txt
"""


import click

@click.command()
@click.option(
    "--dataBase",
    "data_base_file",
    required=True,
    multiple=False,
    type=str,
    default="JASPAR2022_CORE_vertebrates_non-redundant_pfms_meme_nice.txt",
    help="e.g. JASPAR",
)
@click.option(
    "--motif",
    "motifsToCheck",
    required=True,
    multiple=True,
    type=str,
    help="motifs to check",
)
@click.option(
    "--output",
    "out",
    required=True,
    multiple=False,
    type=str,
    help="output file",
)
def cli(data_base_file, motifsToCheck, out):
    
    # read data base entries
    data_base=open(data_base_file, "r")
    data_base_entries=data_base.read().split("MOTIF")
    data_base.close()


    outfile=open(str(out), "w")

    # define alphabet etc.
    outfile.write(data_base_entries[0])

    

    # look for particular motif and write them to output file
    for i in range (0, len(data_base_entries), 1):
        #print (data_base_entries[i])
        if data_base_entries[i].split("\n")[0].split()[0] in motifsToCheck:
            outfile.write("MOTIF"+str(data_base_entries[i]))
    outfile.close()

            
if __name__ == "__main__":
    cli()
