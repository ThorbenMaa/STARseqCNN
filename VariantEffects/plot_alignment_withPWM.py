"""
Description:       
                    

Example commands:   python VariantEffects/plot_alignment_withPWM.py
Outputs:            fasta sequence file with sequences cooresponding to variants with top x % of all variant effects. CSV with haplo seqs and variant effect

"""

#set seed and load dependencies
import click
import logomaker
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



@click.command()
@click.option(
    "--aligMotifFile",
    "alig_motif_file",
    required=True,
    multiple=False,
    type=str,
    default="results/diffTeloHEAC_CTRL_vs_6h/MPRAlm_Alec_diffTeloHEAC_CTRL_vs_6h_withMotifs.csv",
    help="e.g. starrseq-all-final-toorder_oligocomposition.csv",
)
@click.option(
    "--PWMsFile",
    "PWM_file",
    required=True,
    multiple=False,
    type=str,
    default="results/diffTeloHEAC_CTRL_vs_6h/TeloHAEC_CTRL_sigPWMs.txt",
    help="e.g. 2023-01-10_22-29-33 myCounts.minDNAfilt.depthNorm.keepHaps - starr.haplotypes.oligo1.txt",
)
@click.option(
    "--chunk_size",
    "chunk_size",
    required=True,
    multiple=False,
    type=int,
    default=5,
    help="subfigs",
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
def cli(alig_motif_file, PWM_file, out, chunk_size ):#PWM_file

    #import labels; abs diff activities of seq with ID1 - ID2
    df_aligned=pd.read_csv(alig_motif_file, sep=",", low_memory=False)
    print("number of variant pairs :", df_aligned.shape[0])

    data_base=open(PWM_file, "r")
    data_base_entries=data_base.read().split("MOTIF")
    data_base.close()

    # batch variant motif matches
    chuckNumber=1
    for chunk in np.array_split(df_aligned, chunk_size):
        fig, axs= False, False
        chunk=chunk.reset_index(drop=True)
        chuckNumber=chuckNumber+1
        height_ratio_list=[]
        for i in range (0, chunk.shape[0]*3, 3):
            height_ratio_list.append(3)
            height_ratio_list.append(0.5)
            height_ratio_list.append(0.5)
        fig, axs = plt.subplots(nrows= chunk.shape[0]*3, ncols=1, sharex=True, sharey=False, height_ratios=height_ratio_list, figsize=(16, 3*chunk.shape[0] ))
        if chunk.shape[0]==1: # bug repair, da sonst bei nur einem motif ax[i] mit i=0 nicht gefunden wird. Dann ist axs keine liste sondern nur ein ding
            axs=[axs]
        #fig.subplots_adjust(hspace=0.5)
        for index,row in chunk.iterrows():
            print(index)
            index=index*3
            reverse=False
            title="bla"
            strand=row["strand"]
            start=int(row["start"])-1
            stop=int(row["stop"])-1
            total=len(row["Seq1"])
            variantPos=int(row["variantPos"])-1
            print (total)

            for p in range (0, len(data_base_entries), 1):
                #print(data_base_entries[p].split("\n")[0].split()[0])
                if data_base_entries[p].split("\n")[0].split()[0] == str(row["MotifID"]):

                    plot_logo(data_base_entries[p], axs, index+0, "Motif "+data_base_entries[p].split("\n")[0].split()[0], start, stop, total, variantPos, strand)
                    break

            # seq 1 to pwm and plot

            plot_seq_logo(row["Seq1"], axs, index+1,  "ID: "+str(row["ID1"]),  reverse, False, start, stop, variantPos)

            # seq 2 to pwm and plot
            plot_seq_logo(row["Seq2"], axs, index+2,  "ID: "+str(row["ID2"]),  reverse, True, start, stop, variantPos)

            
        #fig.tight_layout()         
        plt.savefig(str(out)+str(chuckNumber)+".svg")
        plt.close()

            
        
        
        
        
        
        
def plot_logo(data_base_entry, axs, axsrow, title, start, stop, total, variantPos, strand):
    

    #prepare fwd seq
    temp_PWM=[[],[],[],[]]
    for i in range (0, start, 1):
        temp_PWM[0].append(0.25)
        temp_PWM[1].append(0.25)
        temp_PWM[2].append(0.25)
        temp_PWM[3].append(0.25)

    for k in range (2, len(data_base_entry.split("\n"))-2, 1): #2 bis -2 damit header und fuss weg fallen
        if data_base_entry.split("\n")[k].split()[0]!="URL":
            temp_PWM[0].append(float(data_base_entry.split("\n")[k].split()[0]))
            temp_PWM[1].append(float(data_base_entry.split("\n")[k].split()[1]))
            temp_PWM[2].append(float(data_base_entry.split("\n")[k].split()[2]))
            temp_PWM[3].append(float(data_base_entry.split("\n")[k].split()[3]))
    
    for i in range (stop, total, 1):
        temp_PWM[0].append(0.25)
        temp_PWM[1].append(0.25)
        temp_PWM[2].append(0.25)
        temp_PWM[3].append(0.25)
                                
    background = 0.25
    cwm_fwd = np.array(temp_PWM) - background

    #prepare rev seq
    temp_PWM=[[],[],[],[]]
    for i in range (0, start, 1):
        temp_PWM[0].append(0.25)
        temp_PWM[1].append(0.25)
        temp_PWM[2].append(0.25)
        temp_PWM[3].append(0.25)

    for k in range (len(data_base_entry.split("\n"))-3, 1, -1): #2 bis -2 damit header und fuss weg fallen
        if data_base_entry.split("\n")[k].split()[0]!="URL":
            temp_PWM[0].append(float(data_base_entry.split("\n")[k].split()[3]))
            temp_PWM[1].append(float(data_base_entry.split("\n")[k].split()[2]))
            temp_PWM[2].append(float(data_base_entry.split("\n")[k].split()[1]))
            temp_PWM[3].append(float(data_base_entry.split("\n")[k].split()[0]))
    
    for i in range (stop, total, 1):
        temp_PWM[0].append(0.25)
        temp_PWM[1].append(0.25)
        temp_PWM[2].append(0.25)
        temp_PWM[3].append(0.25)
    
    
    
    
    
    
    cwm_rev = np.array(temp_PWM) - background

    cwm_fwd=cwm_fwd.tolist()
    cwm_rev=cwm_rev.tolist()
                                    

    df_PWM_f=pd.DataFrame()
    df_PWM_f["A"]=cwm_fwd[0]
    df_PWM_f["C"]=cwm_fwd[1]
    df_PWM_f["G"]=cwm_fwd[2]
    df_PWM_f["T"]=cwm_fwd[3]

    df_PWM_r=pd.DataFrame()
    df_PWM_r["A"]=cwm_rev[0]
    df_PWM_r["C"]=cwm_rev[1]
    df_PWM_r["G"]=cwm_rev[2]
    df_PWM_r["T"]=cwm_rev[3]

                                    #print (df_PWM)
                                
    if strand=="+":
                                    # create forward Logo object
        PWM_logo = logomaker.Logo(df_PWM_f,
                                shade_below=.5,
                                fade_below=.5,
                                ax=axs[axsrow]
                                )
                                    
        # style using Logo methods
        PWM_logo.style_spines(visible=False)
        PWM_logo.style_spines(spines=['left'], visible=True)
        PWM_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

        # style using Axes methods
        PWM_logo.ax.set_ylabel(title + " fwd" , fontsize=3, rotation="horizontal", labelpad=10)
        PWM_logo.ax.yaxis.set_tick_params(labelsize=3) 
        PWM_logo.ax.xaxis.set_ticks_position('none')
        PWM_logo.highlight_position_range(pmin=start, pmax=stop, color='silver')
        PWM_logo.highlight_position_range(pmin=variantPos, pmax=variantPos, color='cyan')
        #PWM_logo.ax.xaxis.set_tick_params(pad=-1)   
        #PWM_logo.ax.set_title(title + " fwd", fontsize=3) 
    elif strand=="-":
        # create rev Logo object
        PWM_logo = logomaker.Logo(df_PWM_r,
                                shade_below=.5,
                                fade_below=.5,
                                ax=axs[axsrow]
                                )

        # style using Logo methods
        PWM_logo.style_spines(visible=False)
        PWM_logo.style_spines(spines=['left'], visible=True)
        PWM_logo.style_xticks(rotation=90, fmt='%d', anchor=0) 
        
        # style using Axes methods
        PWM_logo.ax.set_ylabel(title + " rev" , fontsize=3, rotation="horizontal", labelpad=10)
        PWM_logo.ax.yaxis.set_tick_params(labelsize=3) 
        PWM_logo.ax.xaxis.set_ticks_position('none')
        PWM_logo.highlight_position_range(pmin=start, pmax=stop, color='silver')
        PWM_logo.highlight_position_range(pmin=variantPos, pmax=variantPos, color='cyan')
        #PWM_logo.ax.xaxis.set_tick_params(pad=-1, labelsize=3)    
        #PWM_logo.ax.set_title(title + " rev", fontsize=3)        
        
def plot_seq_logo(seq, axs, axsrow, title, reverse, last, start, stop, variantPos): 
    seq1PWM=[[],[],[],[]] 
    letterSize=0.5
    for base in seq:
        if base == "A":
            seq1PWM[0].append(letterSize)
            seq1PWM[1].append(0)
            seq1PWM[2].append(0)
            seq1PWM[3].append(0)
        if base == "C":
            seq1PWM[0].append(0)
            seq1PWM[1].append(letterSize)
            seq1PWM[2].append(0)
            seq1PWM[3].append(0)
        if base == "G":
            seq1PWM[0].append(0)
            seq1PWM[1].append(0)
            seq1PWM[2].append(letterSize)
            seq1PWM[3].append(0)
        if base == "T":
            seq1PWM[0].append(0)
            seq1PWM[1].append(0)
            seq1PWM[2].append(0)
            seq1PWM[3].append(letterSize)
    

    background = 0 #0.25
    cwm_fwd = np.array(seq1PWM) - background
    cwm_rev = cwm_fwd[::-1, ::-1]

    cwm_fwd=cwm_fwd.tolist()
    cwm_rev=cwm_rev.tolist()
                                        

    df_PWM_f=pd.DataFrame()
    df_PWM_f["A"]=cwm_fwd[0]
    df_PWM_f["C"]=cwm_fwd[1]
    df_PWM_f["G"]=cwm_fwd[2]
    df_PWM_f["T"]=cwm_fwd[3]

    df_PWM_r=pd.DataFrame()
    df_PWM_r["A"]=cwm_rev[0]
    df_PWM_r["C"]=cwm_rev[1]
    df_PWM_r["G"]=cwm_rev[2]
    df_PWM_r["T"]=cwm_rev[3]


    PWM_logo = logomaker.Logo(df_PWM_f,
                            shade_below=.5,
                            fade_below=.5,
                            ax=axs[axsrow],
                            #color_scheme='dimgray'
                            )
                                
    # style using Logo methods
    PWM_logo.style_spines(visible=False)
    if last == True:
        PWM_logo.style_spines(spines=['bottom'], visible=True)
        #PWM_logo.ax.set_ylabel(title, labelpad=-1, fontsize=3)
    #else:
        #PWM_logo.ax.set_ylabel(title, labelpad=-1, fontsize=3)
        #PWM_logo.style_spines(spines=['left'], visible=True)
    PWM_logo.style_xticks(rotation=90, fmt='%d', anchor=0)
    PWM_logo.ax.set_ylabel(title,  fontsize=3, rotation="horizontal", labelpad=20)
    # style using Axes methods
    
    PWM_logo.ax.yaxis.set_tick_params(left=False, labelsize=3) 
    PWM_logo.ax.set(yticklabels=[])
    PWM_logo.ax.xaxis.set_ticks_position('none')
    PWM_logo.ax.xaxis.set_tick_params(pad=-1, labelsize=2)   
    #PWM_logo.ax.set_title(title + " fwd", fontsize=3) 
    PWM_logo.highlight_position_range(pmin=start, pmax=stop, color='silver')
    PWM_logo.highlight_position_range(pmin=variantPos, pmax=variantPos, color='cyan')

    # create rev Logo object
    if reverse==True:
        PWM_logo = logomaker.Logo(df_PWM_r,
                                shade_below=.5,
                                fade_below=.5,
                                ax=axs[axsrow+1]
                                )

        # style using Logo methods
        PWM_logo.style_spines(visible=False)
        PWM_logo.style_spines(spines=['left', 'bottom'], visible=True)
        PWM_logo.style_xticks(rotation=90, fmt='%d', anchor=0) 
        
        # style using Axes methods
        #PWM_logo.ax.set_ylabel("position frequency", labelpad=-1, fontsize=3)
        #PWM_logo.ax.yaxis.set_tick_params(labelsize=3) 
        PWM_logo.ax.xaxis.set_ticks_position('none')
        PWM_logo.ax.xaxis.set_tick_params(pad=-1)    
        PWM_logo.ax.set_title(title + " rev", fontsize=3)

            
if __name__ == "__main__":
    #aligned_motif_file, PWM_file, out = "MPRAlm_Alec_diffTeloHEAC_CTRL_vs_6h_withMotifs.csv", "TeloHAEC_CTRL_sigPWMs.txt", "bla"
    cli()