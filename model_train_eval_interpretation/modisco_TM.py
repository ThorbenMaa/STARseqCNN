

import modiscolite
import numpy as np
import h5py
import click

@click.command()
@click.option(
    "--contrib_file",
    "contrib_file",
    required=True,
    multiple=False,
    type=str,
    default="hypothetical_contribution_scores_mean_HASMC_Chol.npz",
    help="e.g. hypothetical_contribution_scores_mean_HASMC_Chol.npz",
)
@click.option(
    "--seq_file",
    "seq_file",
    required=True,
    multiple=False,
    type=str,
    default="Sequences.npz",
    help="e.g. Sequences.npz",
)
def cli(contrib_file, seq_file):
        #raw model scores
        #print(str(sys.argv[1]))
        scores = np.load(str(contrib_file)) #"hypothetical_contribution_scores_mean_HASMC_Chol.npz"
        scores= scores.f.arr_0
        scores=scores.astype(np.float32)
        for i in range (0, scores.shape[0], 1):
                for j in range (0, scores.shape[1], 1):
                        #center around 0
                        scores[i,j, :]=(scores[i,j,:]-np.mean(scores[i,j,:]))

        print(scores)
        print(np.array(scores).shape)
        print(scores.dtype)

        #my seqs final
        seqs = np.load(str(seq_file)) #"Sequences.npz"
        seqs= seqs.f.arr_0
        seqs=seqs.astype('bool')
        print(seqs)
        print (np.array(seqs).shape)
        print(seqs.dtype)

        #modisco
        pos_patterns, neg_patterns = modiscolite.tfmodisco.TFMoDISco(
                sliding_window_size=8,
                flank_size=8,
                min_metacluster_size=20, 
                target_seqlet_fdr=0.1,
                hypothetical_contribs=scores,
                one_hot=seqs,
                max_seqlets_per_metacluster=20000,
                trim_to_window_size=10,#20
                n_leiden_runs=2,
                initial_flank_to_add=3,#5
                final_min_cluster_size=30, 
                verbose=True)

        #save results
        modiscolite.io.save_hdf5("modisco_results_v2"+str(contrib_file)+".h5", pos_patterns, neg_patterns)

        #print pos and neg patterns; None or empty array if no patterns are identified.
        print(pos_patterns)
        print(neg_patterns)

if __name__ == "__main__":
    cli()
