"""
Description:            This script extracts PWMs from e.g. JASPAR data base to pass them to MEME FIMO in a subsequent step

Output:                 Motif PWMs of interest from JASPAR

example bash command:   python ./allMotifsWithSignificantEffects/trimMotifs.py --outputPWM test.txt
"""


import click
import h5py
import numpy as np

@click.command()
@click.option(
    "--reportFile",
    "report_file",
    required=True,
    multiple=True,
    type=str,
    default=["modisco_resultshypothetical_contribution_scores_mean_diffTeloHEAC_CTRL_vs_6h.npz.h5"],
    help="e.g. modico output",
)
@click.option(
    "--outputPWM",
    "out_pwm",
    required=True,
    multiple=False,
    type=str,
    help="output PWM file",
)
def cli(report_file, out_pwm):
    motifCounter=0
    f = open(out_pwm, 'w')
    background=[0.25, 0.25, 0.25, 0.25]
        
    f.write('MEME version 4\n\n')
    f.write('ALPHABET= ACGT\n\n')
    f.write('strands: + -\n\n')
    f.write('Background letter frequencies (from unknown source):\n')
    f.write('A %.3f C %.3f G %.3f T %.3f\n\n' % tuple(list(background)))
    for k in range (0, len(report_file), 1):
        # read data base entries
        report=h5py.File(report_file[k] , 'r')


        # trim pos patterns
        trimmed_ppms=[]
        #trimmed_hcwms=[]
        for key1 in report.keys(): #pos patterns and neg patterns
            for key2 in report[key1].keys():  
                pwm=report[key1][key2]["sequence"][:]
                cwm=(report[key1][key2]["hypothetical_contribs"][:])

                score = np.sum(np.abs(cwm), axis=1)
                trim_thresh = np.max(score) * 0.3  # Cut off anything less than 30% of max score
                pass_inds = np.where(score >= trim_thresh)[0] #returns indices that are non-zero
                trimmed = pwm[np.min(pass_inds): np.max(pass_inds) + 1]

                trimmed_ppms.append(trimmed)
        report.close()
        
        #write MEME file with trimmed PWM motifs


        for i in range (0, len(trimmed_ppms), 1):
            f.write('\nMOTIF '+ str(motifCounter) +'\n')
            motifCounter = motifCounter + 1
            f.write('letter-probability matrix: alength= 4 w= %d nsites= 1 E= 0e+0\n' % trimmed_ppms[i].shape[0])
            for s in trimmed_ppms[i]:
                f.write('%.5f %.5f %.5f %.5f\n' % tuple(s))
    f.close()
    
    
    
    """
    #write MEME file with trimmed hCWM motifs
    background=[0.25, 0.25, 0.25, 0.25]
    f = open(out_hcwm, 'w')
    f.write('MEME version 4\n\n')
    f.write('ALPHABET= ACGT\n\n')
    f.write('strands: + -\n\n')
    f.write('Background letter frequencies (from unknown source):\n')
    f.write('A %.3f C %.3f G %.3f T %.3f\n\n' % tuple(list(background)))

    for i in range (0, len(trimmed_ppms), 1):
        f.write('\nMOTIF '+ str(i) +'\n')
        f.write('letter-probability matrix: alength= 4 w= %d nsites= 1 E= 0e+0\n' % trimmed_ppms[i].shape[0])
        for s in trimmed_ppms[i]:
            f.write('%.5f %.5f %.5f %.5f\n' % tuple(s))
    f.close()
    """
         
    


            
if __name__ == "__main__":
    cli()
