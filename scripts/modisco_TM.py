"""
Description:            Script that laods .npz files generated with "ism_TM.py", calulcates per nucleotide hypothetical contribution scores (see also Agrawal and Kelley, 2022 (PMID: 36419176)) from the raw CNN scores, runs tfmodisco-lite, and generates reports.

Inputs:                 Input 1 is file with raw CNN scores as .npz file with numpy array of shape number of sequences x sequence length x 4.  
                        Input 2 is file with one-hot-encoded sequences using the ISM approach (see "ism_TM.py") as .npz file with numpy array of shape number of sequences x sequence length x 4. Both inputs
                        can be gereated using "ism_TM.py".

Output:                 modisco results in as .h5 file

Example command:        python ./scripts/modisco_TM.py ./scripts/hypothetical_contribution_scores_mean_HASMC_Chol.npz Sequences.npz

"""



import modiscolite
import sys
import numpy as np
import h5py

#raw model scores
print(str(sys.argv[1]))
scores = np.load(str(sys.argv[1])) #"hypothetical_contribution_scores_mean_HASMC_Chol.npz"
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
seqs = np.load(str(sys.argv[1])) #"Sequences.npz"
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
        trim_to_window_size=20,
        n_leiden_runs=2,
        initial_flank_to_add=5,
        final_min_cluster_size=30, 
        verbose=True)

#save results
modiscolite.io.save_hdf5("modisco_results"+str(sys.argv[1])+".h5", pos_patterns, neg_patterns)

#print pos and neg patterns; None or empty array if no patterns are identified.
print(pos_patterns)
print(neg_patterns)


