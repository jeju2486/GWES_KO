import numpy as np
from scipy.stats import fisher_exact, chi2
from scipy.optimize import minimize
from scipy.linalg import expm
from Bio import Phylo, SeqIO
import csv
import argparse
import os

def likelihood(params, site1, site2, root1, root2, edgeLengths):
    if len(params) != 5:
        print(f"Error in likelihood: expected 5 parameters but got {len(params)}: {params}")
        return -1e10  # Return a large negative number to signal an error

    # Convert parameters to float
    a1, a2, b1, b2, eps = map(float, params)

    Q = np.array([
        [-(a1+b1), a1, b1, 0], 
        [a2, -(a2+b1*eps), 0, b1*eps], 
        [b2, 0, -(b2+a1*eps), a1*eps], 
        [0, b2/eps, a2/eps, -(b2/eps+a2/eps)]
    ])

    Q2 = np.dot(Q, Q)
    Q3 = np.dot(Q2, Q)
    
    # Define the mapping matrix
    map_ = np.array([
        [0,8,2,10], 
        [4,12,6,14], 
        [1,9,3,11], 
        [5,13,7,15]
    ])

    Id = np.diag([1, 1, 1, 1])
    
    ret = -1e-10
    
    epsilon = 1e-10  # Small constant to ensure numerical stability

    # Ensure combined is a NumPy array for efficient indexing
    combined = np.array([map_[site1[i], site2[i]] for i in range(len(site1))])
    
#    print(f"combined : {combined} Q.flatten()[combined] : {Q.flatten()[combined]}")
    
    # Calculate the power terms for edgeLengths
    edgeLengths_squared = edgeLengths**2
    edgeLengths_cubed = edgeLengths**3
    
#    print(f"Q shape : {Q.flatten()[combined]}, edge shape : {edgeLengths.shape}, edge : {edgeLengths}")
    
    # Calculate log sum terms using vectorized operations
    log_sum_terms = (
        Id.flatten()[combined]
        + Q.flatten()[combined] * edgeLengths
        + Q2.flatten()[combined] * edgeLengths_squared / 2
        + Q3.flatten()[combined] * edgeLengths_cubed / 6
        + epsilon
    )

#    print(f"log_sum = {np.log(log_sum_terms)}")

    # Accumulate the log values
    ret += np.sum(np.log(log_sum_terms))
        
    return -ret

def likelihood0(params, site1, site2, root1, root2, edgeLengths):
    if len(params) != 4:
        print(f"Error in likelihood0: expected 4 parameters but got {len(params)}: {params}")
        return -1e10  # Return a large negative number to signal an error
    # Convert params to list and append 1, then call likelihood
    augmented_params = list(params) + [1]
    return likelihood(augmented_params, site1, site2, root1, root2, edgeLengths)

#Read the tree file and get the edges name and lengths
def get_edges_lengths(tree, seq_dict):
    edges = []
    edge_lengths = []
    for clade in tree.find_clades(order='level'):
        for child in clade.clades:
            if clade.name in seq_dict and child.name in seq_dict:
                edges.append((seq_dict[clade.name], seq_dict[child.name]))
                
                # Add branch length 1e-7 if it's None or 0
                if child.branch_length is None or child.branch_length == 0:
                    edge_lengths.append(1e-7)
                else:
                    edge_lengths.append(child.branch_length)
                    
    return edges, edge_lengths
    

def main(treefile, sequencefile, output_dir):
    print("Reading tree and sequences...")
    tree = Phylo.read(treefile, 'newick')
    seqs = list(SeqIO.parse(sequencefile, 'fasta'))

    labs = [seq.id for seq in seqs]
    print(f"Number of sequences: {len(seqs)}")
    
    # Create a dictionary to map sequence IDs to their indices
    seq_dict = {seq.id: idx for idx, seq in enumerate(seqs)}
    print(f"Sequence dictionary: {seq_dict}")
    
    edges, edge_lengths = get_edges_lengths(tree, seq_dict)
    edges = np.array(edges)
    edge_lengths = np.array(edge_lengths)
    print(f"Number of edges: {len(edges)}")
    max_length = max(edge_lengths)
    print(f"Max length: {max_length}")

    if len(edges) == 0:
        print("No valid edges found. Exiting.")
        return

    l = len(seqs[0].seq)
    print(f"Sequence length: {l}")
    
    mut_file= os.path.join(output_dir, 'mut.txt')
    mut_types_file = os.path.join(output_dir, 'mut_type.txt')
    
    # Load or compute mutation types
    if os.path.exists(str(mut_file)+".npy") and os.path.exists(str(mut_types_file)+".npy"):
        print(f"Loading mutation types from {mut_file} and {mut_types_file}...")
        muts = np.load(str(mut_file)+".npy")
        mutsTyp = np.load(str(mut_types_file)+".npy")
        print("Loading finished")

    else:
        print("Computing mutation types...")
        muts = np.zeros((l, len(edges)), dtype=bool)
        for i in range(l):
            if i % 1000 == 0:
                print(f"Processing site {i+1} of {l}...")
        
            for b in range(len(edges)):
                parent_base = seqs[edges[b, 0]].seq[i].lower()
                child_base = seqs[edges[b, 1]].seq[i].lower()
                if parent_base != child_base:
                    muts[i, b] = True

        print("Determining types of mutations...")
        mutsTyp = np.zeros((l, len(edges)), dtype=int)
        for i in range(l):
            if i % 100 == 0:
                print(f"Processing site {i+1} of {l}...")
        
            for b in range(len(edges)):
                if edges[b, 0] >= len(seqs) or edges[b, 1] >= len(seqs):
                    continue
                parent_base = seqs[edges[b, 0]].seq[i].lower()
                child_base = seqs[edges[b, 1]].seq[i].lower()
        
                if parent_base == 'a' and child_base == 'a':
                    mutsTyp[i, b] = 0
                elif parent_base == 'c' and child_base == 'a':
                    mutsTyp[i, b] = 1
                elif parent_base == 'a' and child_base == 'c':
                    mutsTyp[i, b] = 2
                elif parent_base == 'c' and child_base == 'c':
                    mutsTyp[i, b] = 3

        print("Mutation type determination completed.")
        np.save(mut_file, muts)
        np.save(mut_types_file, mutsTyp)
        
    # Create sequence matrix
    seq_matrix = np.zeros((len(seqs), l))
    for idx, seq in enumerate(seqs):
        seq_matrix[idx] = [1 if base.lower() == 'c' else -1 for base in seq.seq]

    print("Calculating scores...")
    
    # Debugging prints after calculating scores
    print(f"Edges: {edges[1:10,]}")
    print(f"Edge lengths: {edge_lengths}")
    print(f"Mutation types shape: {mutsTyp.shape}")
    print(f"Mutation types: {mutsTyp[:5]}")
    print(f"Mutations shape: {muts.shape}")
    print(f"Mutations: {muts[:5]}") 

    seq_count = len(edges)
    results = []
    s = 0
    tab = np.zeros((2, 2))
    npat = mutsTyp.shape[0]
    store = np.ones(200 * 200 * 200)
    
    # Update the bounds for each parameter
    param_bound = [(1e-5, None)] * 4 + [(1.0, None)]
    
    for i in range(npat):
        if i % 1000 == 0:
            print(f"Processing site {i+1} of {npat}...")
        # Identify columns where muts[i, :] is True
        col_indices = np.where(muts[i, :])[0]
        
            # Create the subset matrix
        muts_subset = muts[i:npat, col_indices]
        
        # Calculate the row sums
        col_sums = np.sum(muts_subset, axis=1)
        
        # Find columns where row_sums >= 4
    
        w = np.where(col_sums >= 4)[0] + i
        w = w[w != i]
        
        print(f"w: {w} and length: {len(w)}")
        
        for j in w:
            print(f"Comparing sites {i} and {j}")
            seq_mult = np.dot(seq_matrix[:, i], seq_matrix[:, j])
            mi = seq_mult / seq_count
            print(f"mi: {mi}")
            if mi < 0.3:
                print(f"skip as {round(mi,2)} < -0.1")
                continue
            m = muts[i, :] * 2 + muts[j, :]
            tab[0, 0] = np.sum(m == 3)
            tab[0, 1] = np.sum(m == 1)
            tab[1, 0] = np.sum(m == 2)
            tab[1, 1] = np.sum(m == 0)
            if tab[0, 0] > 199 or tab[0, 1] > 199 or tab[1, 0] > 199:
                score = round(np.log10(fisher_exact(tab)[1]),2)
            else:
                ind = int(tab[0, 0] * 40000 + tab[0, 1] * 200 + tab[1, 0])
                if store[ind] == 1:
                    store[ind] = round(np.log10(fisher_exact(tab)[1]),2)
                score = store[ind]

            site1 = mutsTyp[i, :]
            site2 = mutsTyp[j, :]
            
            root_candidates = set(edges[:, 0]) - set(edges[:, 1])

            if not root_candidates:
                print(f"No root found for positions {i} and {j}. Skipping.")
                continue
            
            root = list(root_candidates)[0]
            bra = np.where(edges[:, 0] == root)[0][0]

            if site1[bra] == 1 or site1[bra] == 2:
                root1 = 1
            else:
                root1 = 2
            if site2[bra] == 1 or site2[bra] == 2:
                root2 = 1
            else:
                root2 = 2

            score2 = 0
            
#            print(f"site1: {site1}, site2: {site2}, root1: {root1}, root2:{root2} and edge_lengths: {edge_lengths}")

            try:
                fit0 = minimize(lambda params: likelihood0(params, site1, site2, root1, root2, edge_lengths), x0=[1, 1, 1, 1], method='L-BFGS-B', bounds=[(1e-5, None)] * 4)
                fit = minimize(lambda params: likelihood(params, site1, site2, root1, root2, edge_lengths), x0=[fit0.x[0], fit0.x[1], fit0.x[2], fit0.x[3], 1], method='L-BFGS-B', bounds=param_bound)
                lrt = 2 * (-fit.fun + fit0.fun)
                print(f"first model : {fit0.fun}, second model : {fit.fun} with param : {fit.x[0],fit.x[1],fit.x[2],fit.x[3]} eps : {fit.x[4]}, and lrt : {lrt}")
                if lrt < 0:
                    print(f"Negative LRT value: {lrt}, setting to 0")
                    lrt = 0  # To handle cases where lrt is negative due to numerical issues
                p_value = chi2.sf(lrt, df=1)
                
                if p_value <= 0:
                    print(f"Non-positive p-value: {p_value}, setting to smallest positive float")
                    p_value = np.finfo(float).tiny  # Smallest positive float to avoid log10(0)
                    
                score2 = round(np.log10(p_value),2)
                print(f"score: {score} and score2: {score2}")
                
            except Exception as e:
                print(f"Error optimizing likelihood for positions {i} and {j}: {e}")
                continue

            if score <= -10 or score2 <= -3:
                s += 1
                results.append([i, j, tab[0, 0], tab[0, 1], tab[1, 0], tab[1, 1], score, score2])
                print(f"results appended for {i} and {j}. New results: {[i, j, tab[0, 0], tab[0, 1], tab[1, 0], tab[1, 1], score, score2]}")


    print(f"Sorting and saving {len(results)} significant results...")
    results.sort(key=lambda x: (x[6], x[7]))  # Sorting by score and score2

    output_file = os.path.join(output_dir, 'score.tsv')
    with open(output_file, 'w', newline='') as fout:
        writer = csv.writer(fout, delimiter='\t')
        writer.writerow(['i', 'j', 'tab[1,1]', 'tab[1,2]', 'tab[2,1]', 'tab[2,2]', 'score', 'score2'])
        writer.writerows(results)

    print("Process completed successfully.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run CTMC analysis.')
    parser.add_argument('-t', '--treefile', required=True, help='Path to the Newick tree file.')
    parser.add_argument('-s', '--sequencefile', required=True, help='Path to the sequence file in FASTA format.')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory to save the output files.')
   
    # do whatever you do
    args = parser.parse_args()
    main(args.treefile, args.sequencefile, args.output_dir)

        
