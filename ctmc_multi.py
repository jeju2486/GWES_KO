import numpy as np
from scipy.stats import fisher_exact, chi2
from scipy.optimize import minimize
from scipy.linalg import expm
from Bio import Phylo, SeqIO
import csv
import argparse
import os
from concurrent.futures import ThreadPoolExecutor, as_completed

def likelihood(params, site1, site2, root1, root2, edgeLengths):
    if len(params) != 5:
        print(f"Error in likelihood: expected 5 parameters but got {len(params)}: {params}")
        return -1e10  # Return a large negative number to signal an error

    # Convert parameters to float
    a1, a2, b1, b2, eps = map(float, params)

    Q0 = np.array([
        [0, a1, b1, 0], 
        [a2, 0, 0, b1 * eps], 
        [b2, 0, 0, a1 * eps], 
        [0, b2 / eps, a2 / eps, 0]
    ])
    
    Q = Q0 - np.diag(np.sum(Q0, axis=1))

    combined = expm(Q * 3) 

    map_ = np.array([
        [0, 1, 2, 3], 
        [4, 5, 6, 7], 
        [8, 9, 10, 11], 
        [12, 13, 14, 15]
    ]).flatten()
    
    combined = np.array([map_[site1[i] * 4 + site2[i]] for i in range(len(site1))])

    Id = np.diag([1, 1, 1, 1])

    den = a1 * b1 * eps * eps + a1 * b2 + a2 * b1 + a2 * b2
    if den <= 0:
        print(f"Denominator became non-positive with params: {params}")
        return -1e10  # Return a large negative number to signal an error

    if root1 == 1 and root2 == 1:
        num = a2 * b2
    elif root1 == 1 and root2 == 2:
        num = a2 * b1
    elif root1 == 2 and root2 == 1:
        num = a1 * b2
    elif root1 == 2 and root2 == 2:
        num = a1 * b1 * eps * eps

    if num <= 0:
        print(f"Numerator became non-positive with params: {params}")
        return -1e10  # Return a large negative number to signal an error

    ret = np.log(num / den)
    if np.isnan(ret) or np.isinf(ret):
        print(f"Log result became NaN or infinite with params: {params}, num: {num}, den: {den}")
        return -1e10  # Return a large negative number to signal an error

    for idx, P in enumerate(combined):
        P_sum = np.sum(Id + P + 1)
        if P_sum <= 0:
            print(f"Log argument became non-positive with params: {params}, P: {P}, Id: {Id}")
            return -1e10  # Return a large negative number to signal an error
        log_sum = np.log(Id + P + 1)
        if np.isnan(log_sum).any() or np.isinf(log_sum).any():
        
            print(f"Log result became NaN or infinite with params: {params}, P:{P}, Id: {Id}, log_sum: {log_sum}")
            return -1e10  # Return a large negative number to signal an error
        ret += np.sum(log_sum)

    if np.isnan(ret) or np.isinf(ret):
        print(f"Final result became NaN or infinite with params: {params}")
        ret = -1e10
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
                edge_lengths.append(child.branch_length if child.branch_length is not None else 1e-7)
    return edges, edge_lengths

def process_mutation_types(seqs, edges, i):
    muts = np.zeros((len(edges)), dtype=bool)
    for b in range(len(edges)):
        parent_base = seqs[edges[b, 0]].seq[i].lower()
        child_base = seqs[edges[b, 1]].seq[i].lower()
        if parent_base != child_base:
            muts[b] = True
    return muts

def determine_mutation_type(seqs, edges, i):
    mutsTyp = np.zeros((len(edges)), dtype=int)
    for b in range(len(edges)):
        if edges[b, 0] >= len(seqs) or edges[b, 1] >= len(seqs):
            continue
        parent_base = seqs[edges[b, 0]].seq[i].lower()
        child_base = seqs[edges[b, 1]].seq[i].lower()
        
        if parent_base == 'a' and child_base == 'a':
            mutsTyp[b] = 0
        elif parent_base == 'a' and child_base == 'c':
            mutsTyp[b] = 1
        elif parent_base == 'c' and child_base == 'a':
            mutsTyp[b] = 2
        elif parent_base == 'c' and child_base == 'c':
            mutsTyp[b] = 3
    return mutsTyp

def main(treefile, sequencefile, output_dir, num_threads):
    
    print(f"Process starts with CPU number of {num_threads}")
    
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
    print(f"Edges: {edges}")
    print(f"Edge lengths: {edge_lengths}")

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
        mutsTyp = np.zeros((l, len(edges)), dtype=int)

        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            futures = {executor.submit(process_mutation_types, seqs, edges, i): i for i in range(l)}
            for future in as_completed(futures):
                i = futures[future]
                muts[i] = future.result()
                if i % 1000 == 0:
                    print(f"Processed site {i+1} of {l} for mutations")

            futures = {executor.submit(determine_mutation_type, seqs, edges, i): i for i in range(l)}
            for future in as_completed(futures):
                i = futures[future]
                mutsTyp[i] = future.result()
                if i % 100 == 0:
                    print(f"Processed site {i+1} of {l} for mutation types")

        np.save(mut_file, muts)
        np.save(mut_types_file, mutsTyp)
    
    print("Calculating scores...")

    results = []
    s = 0
    tab = np.zeros((2, 2))
    npat = mutsTyp.shape[0]
    store = np.ones(200 * 200 * 200)

    def process_site_pair(i, j):
        site_results = []
        m = muts[i, :] * 2 + muts[j, :]
        tab[0, 0] = np.sum(m == 3)
        tab[0, 1] = np.sum(m == 1)
        tab[1, 0] = np.sum(m == 2)
        tab[1, 1] = np.sum(m == 0)
        if tab[0, 0] > 199 or tab[0, 1] > 199 or tab[1, 0] > 199:
            score = round(np.log10(fisher_exact(tab)[1]))
        else:
            ind = int(tab[0, 0] * 40000 + tab[0, 1] * 200 + tab[1, 0])
            if store[ind] == 1:
                store[ind] = round(np.log10(fisher_exact(tab)[1]))
            score = store[ind]

        site1 = mutsTyp[i, :]
        site2 = mutsTyp[j, :]

        root_candidates = set(edges[:, 0]) - set(edges[:, 1])
        if not root_candidates:
            return site_results

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

        try:
            fit0 = minimize(lambda params: likelihood0(params, site1, site2, root1, root2, edge_lengths), x0=[1, 1, 1, 1], method='L-BFGS-B', bounds=[(1e-8, 10.0)] * 4)
            fit = minimize(lambda params: likelihood(params, site1, site2, root1, root2, edge_lengths), x0=[fit0.x[0], fit0.x[1], fit0.x[2], fit0.x[3], 1], method='L-BFGS-B', bounds=[(1e-8, 10.0)] * 5)
            lrt = 2 * (-fit.fun + fit0.fun)
            if lrt < 0:
                lrt = 0  # To handle cases where lrt is negative due to numerical issues
            p_value = chi2.sf(lrt, df=1)
            if p_value <= 0:
                p_value = np.finfo(float).tiny  # Smallest positive float to avoid log10(0)
            score2 = round(np.log10(p_value))
        except Exception as e:
            print(f"Error optimizing likelihood for positions {i} and {j}: {e}")
            return site_results

        if score <= -10 or score2 <= -10:
            site_results.append([i, j, tab[0, 0], tab[0, 1], tab[1, 0], tab[1, 1], score, score2])

        return site_results

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        for i in range(npat):
            if i % 1000 == 0:
                print(f"Processing site {i+1} of {npat}...")
            col_indices = np.where(muts[i, :])[0]
            muts_subset = muts[i:npat, col_indices]
            row_sums = np.sum(muts_subset, axis=0)
            w = np.where(row_sums >= 4)[0]
            original_w = col_indices[w]
            print(f"w = {original_w}")
            for j in original_w:
                print(f"comparing sites {i} and {j}")
                futures.append(executor.submit(process_site_pair, i, j))

        for future in as_completed(futures):
            results.extend(future.result())

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
    parser.add_argument('-T', '--threads', type=int, default=1, help='Number of threads for multithreading (default: 1).')

    args = parser.parse_args()
    main(args.treefile, args.sequencefile, args.output_dir, args.threads)

