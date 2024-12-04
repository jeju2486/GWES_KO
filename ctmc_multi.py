#!/usr/bin/env python3

import numpy as np
from scipy.stats import fisher_exact, chi2
from scipy.optimize import minimize
from Bio import Phylo, SeqIO
import csv
import argparse
import os
import logging
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm  # Import tqdm for progress bars
import math

# Global variables for shared data
global_muts_filtered = None
global_mutsTyp_filtered = None
global_edge_lengths = None
global_map_matrix = None
global_param_bounds = None

def initializer(muts_filtered, mutsTyp_filtered, edge_lengths, map_matrix, param_bounds):
    """
    Initializer function to set global variables in each worker process.
    """
    global global_muts_filtered
    global global_mutsTyp_filtered
    global global_edge_lengths
    global global_map_matrix
    global global_param_bounds

    global_muts_filtered = muts_filtered
    global_mutsTyp_filtered = mutsTyp_filtered
    global_edge_lengths = edge_lengths
    global_map_matrix = map_matrix
    global_param_bounds = param_bounds

def compute_likelihood(params, site1, site2):
    """
    Compute the negative log-likelihood for given parameters using a truncated series expansion up to Q^3.
    Utilizes global variables for shared data.
    """
    if len(params) != 5:
        logging.error(f"compute_likelihood: expected 5 parameters but got {len(params)}: {params}")
        return -1e10  # Signal an error

    # Unpack parameters
    a1, a2, b1, b2, eps = params

    # Construct the Q matrix
    Q = np.array([
        [-(a1 + b1), a1, b1, 0],
        [a2, -(a2 + b1 * eps), 0, b1 * eps],
        [b2, 0, -(b2 + a1 * eps), a1 * eps],
        [0, b2 / eps, a2 / eps, -(b2 / eps + a2 / eps)]
    ])

    # Precompute Q squared and Q cubed
    Q2 = Q @ Q
    Q3 = Q2 @ Q

    # Flatten matrices for advanced indexing
    Q_flat = Q.flatten()
    Q2_flat = Q2.flatten()
    Q3_flat = Q3.flatten()

    # Compute combined indices for mapping
    combined = global_map_matrix[site1, site2]  # Shape: (num_edges,)

    # Compute terms for the series expansion
    edge_lengths_squared = global_edge_lengths ** 2
    edge_lengths_cubed = global_edge_lengths ** 3

    # Compute the log-sum terms using vectorized operations
    log_sum_terms = (
        1.0  # Identity matrix diagonal is 1
        + Q_flat[combined] * global_edge_lengths
        + Q2_flat[combined] * edge_lengths_squared / 2
        + Q3_flat[combined] * edge_lengths_cubed / 6
        + 1e-10  # Small epsilon for numerical stability
    )

    # Check for non-positive values in log_sum_terms
    if np.any(log_sum_terms <= 0):
        logging.debug(f"Non-positive log_sum_terms encountered. Params: {params}")
        logging.debug(f"log_sum_terms: {log_sum_terms}")
        return 1e10  # Signal an error

    # Sum the log probabilities
    try:
        log_likelihood = np.sum(np.log(log_sum_terms))
    except FloatingPointError as e:
        logging.debug(f"FloatingPointError during log computation: {e}")
        return 1e10

    return -log_likelihood  # Negative log-likelihood for minimization

def compute_likelihood_null(params, site1, site2):
    """
    Compute the negative log-likelihood for the null model (eps = 1).
    Utilizes global variables for shared data.
    """
    if len(params) != 4:
        logging.error(f"compute_likelihood_null: expected 4 parameters but got {len(params)}: {params}")
        return -1e10  # Signal an error

    # Append eps=1 for the null model
    augmented_params = list(params) + [1.0]
    return compute_likelihood(augmented_params, site1, site2)

def get_edges_and_lengths(tree, seq_dict):
    """
    Extract edges and their lengths from the tree.
    """
    edges = []
    edge_lengths = []

    for clade in tree.find_clades(order='level'):
        for child in clade.clades:
            if clade.name in seq_dict and child.name in seq_dict:
                parent_idx = seq_dict[clade.name]
                child_idx = seq_dict[child.name]
                edges.append((parent_idx, child_idx))

                # Use a small value if branch length is missing or zero
                branch_length = child.branch_length if child.branch_length and child.branch_length > 0 else 1e-7
                edge_lengths.append(branch_length)

    return np.array(edges, dtype=int), np.array(edge_lengths, dtype=float)

def compute_mutation_matrices(seq_array, edges):
    """
    Compute mutation presence and types matrices.
    """
    # Extract parent and child sequences based on edges
    parent_seqs = seq_array[edges[:, 0], :]  # Shape: (num_edges, sequence_length)
    child_seqs = seq_array[edges[:, 1], :]   # Shape: (num_edges, sequence_length)

    # Transpose to align sequences for each position
    parent_seqs = parent_seqs.T  # Shape: (sequence_length, num_edges)
    child_seqs = child_seqs.T    # Shape: (sequence_length, num_edges)

    # Compute mutation presence: True where parent != child
    muts = parent_seqs != child_seqs  # Shape: (sequence_length, num_edges)

    # Initialize mutation types matrix
    mutsTyp = np.zeros_like(parent_seqs, dtype=int)  # Shape: (sequence_length, num_edges)

    # Define mutation types
    # Assuming only 'A' and 'C' are relevant; adjust as necessary for other nucleotides
    mutsTyp[(parent_seqs == 'A') & (child_seqs == 'A')] = 0
    mutsTyp[(parent_seqs == 'C') & (child_seqs == 'A')] = 1
    mutsTyp[(parent_seqs == 'A') & (child_seqs == 'C')] = 2
    mutsTyp[(parent_seqs == 'C') & (child_seqs == 'C')] = 3

    return muts, mutsTyp

def process_position_pair(args):
    """
    Process a single position pair (pos_i, pos_j) and return the result.
    Utilizes global shared data.
    """
    idx_i, pos_i, idx_j, pos_j = args

    # Access global shared data
    muts_i = global_muts_filtered[idx_i, :]          # Shape: (num_edges,)
    muts_j = global_muts_filtered[idx_j, :]          # Shape: (num_edges,)
    site1 = global_mutsTyp_filtered[idx_i, :]        # Shape: (num_edges,)
    site2 = global_mutsTyp_filtered[idx_j, :]        # Shape: (num_edges,)

    # Construct contingency table
    m = muts_i * 2 + muts_j           # Combine mutations: 3,1,2,0
    contingency_table = np.array([
        [np.sum(m == 3), np.sum(m == 1)],
        [np.sum(m == 2), np.sum(m == 0)]
    ])

    # Perform Fisher's exact test
    try:
        _, p_value_fisher = fisher_exact(contingency_table)
        score_fisher = round(np.log10(p_value_fisher), 2)
    except Exception as e:
        logging.error(f"Error in Fisher's exact test for positions {pos_i} and {pos_j}: {e}")
        return None

    # Perform likelihood optimizations
    try:
        # Null model optimization
        fit_null = minimize(
            compute_likelihood_null,
            x0=[1.0, 1.0, 1.0, 1.0],
            args=(site1, site2),
            method='L-BFGS-B',
            bounds=[global_param_bounds[0]] * 4  # Use bounds for a1, a2, b1, b2
        )

        if not fit_null.success:
            logging.debug(f"Null model optimization did not converge for positions {pos_i} and {pos_j}. Message: {fit_null.message}")
            return None

        # Alternative model optimization
        fit_alt = minimize(
            compute_likelihood,
            x0=list(fit_null.x) + [1.0],  # Starting point based on null model
            args=(site1, site2),
            method='L-BFGS-B',
            bounds=global_param_bounds
        )

        if not fit_alt.success:
            logging.debug(f"Alternative model optimization did not converge for positions {pos_i} and {pos_j}. Message: {fit_alt.message}")
            return None

        # Compute likelihood ratio test statistic
        lrt_stat = 2 * (-fit_alt.fun + fit_null.fun)
        lrt_stat = max(lrt_stat, 0)  # Ensure non-negative

        # Compute p-value from chi-squared distribution with 1 degree of freedom
        p_value_lrt = chi2.sf(lrt_stat, df=1)
        p_value_lrt = max(p_value_lrt, np.finfo(float).tiny)  # Avoid log10(0)

        if np.isnan(p_value_lrt):
            logging.error(f"Computed p_value_lrt is nan for positions {pos_i} and {pos_j}.")
            return None

        score_lrt = round(np.log10(p_value_lrt), 2)
        logging.debug(f"Likelihood ratio test (LRT) statistic: {lrt_stat}")
        logging.debug(f"Scores: Fisher's p-value score = {score_fisher}, LRT p-value score = {score_lrt}")

    except Exception as e:
        logging.error(f"Error in likelihood optimization for positions {pos_i} and {pos_j}: {e}")
        return None

    # Append significant results based on thresholds
    if score_fisher <= -10 or score_lrt <= -3:
        result = [
            int(pos_i), int(pos_j),
            int(contingency_table[0, 0]), int(contingency_table[0, 1]),
            int(contingency_table[1, 0]), int(contingency_table[1, 1]),
            float(score_fisher), float(score_lrt)
        ]
        logging.debug(f"Significant result for positions {pos_i} and {pos_j}: {result}")
        return result
    else:
        return None

def get_chunk_file_name(output_dir, chunk_id):
    """
    Generate the file name for a given chunk.
    """
    return os.path.join(output_dir, f'score_chunk_{chunk_id}.tsv')

def process_chunk(position_pairs, chunk_id, executor, output_dir):
    """
    Process a single chunk of position pairs using the provided executor.
    """
    results = []
    num_pairs = len(position_pairs)
    logging.info(f"Processing chunk {chunk_id} with {num_pairs} position pairs...")

    # Use tqdm to display progress bar over the completed tasks
    for res in tqdm(executor.map(process_position_pair, position_pairs),
                    total=num_pairs, desc=f"Processing Chunk {chunk_id}"):
        if res is not None:
            results.append(res)

    # Sort results by Fisher's score and then by LRT score
    results.sort(key=lambda x: (x[6], x[7]))  # Sorting by score_fisher, then score_lrt

    # Save results to a TSV file for this chunk
    output_file = get_chunk_file_name(output_dir, chunk_id)
    try:
        with open(output_file, 'w', newline='') as fout:
            writer = csv.writer(fout, delimiter='\t')
            writer.writerow(['Position_i', 'Position_j', 'Tab_00', 'Tab_01', 'Tab_10', 'Tab_11', 'Score_Fisher', 'Score_LRT'])
            writer.writerows(results)
        logging.info(f"Chunk {chunk_id} results saved to '{output_file}'.")
    except Exception as e:
        logging.error(f"Error saving chunk {chunk_id} results to '{output_file}': {e}")

def combine_chunk_results(output_dir, max_chunk_id):
    """
    Combine all chunk result files into a single file.
    """
    combined_results_file = os.path.join(output_dir, 'score_combined.tsv')
    with open(combined_results_file, 'w', newline='') as fout_combined:
        writer = csv.writer(fout_combined, delimiter='\t')
        writer.writerow(['Position_i', 'Position_j', 'Tab_00', 'Tab_01', 'Tab_10', 'Tab_11', 'Score_Fisher', 'Score_LRT'])
        for chunk_id in range(1, max_chunk_id + 1):
            chunk_file = get_chunk_file_name(output_dir, chunk_id)
            if os.path.exists(chunk_file):
                with open(chunk_file, 'r') as fin_chunk:
                    reader = csv.reader(fin_chunk, delimiter='\t')
                    next(reader)  # Skip header
                    writer.writerows(reader)
                logging.info(f"Chunk {chunk_id} results appended to '{combined_results_file}'.")
            else:
                logging.warning(f"Chunk file '{chunk_file}' not found.")

def main(treefile, sequencefile, output_dir, num_threads):
    logging.info("Reading tree and sequences...")

    # Read the phylogenetic tree
    try:
        tree = Phylo.read(treefile, 'newick')
    except Exception as e:
        logging.error(f"Error reading tree file '{treefile}': {e}")
        return

    # Read sequences from the FASTA file
    try:
        seqs = list(SeqIO.parse(sequencefile, 'fasta'))
    except Exception as e:
        logging.error(f"Error reading sequence file '{sequencefile}': {e}")
        return

    num_seqs = len(seqs)
    logging.info(f"Number of sequences: {num_seqs}")

    # Validate sequences are of equal length
    sequence_lengths = set(len(seq.seq) for seq in seqs)
    if len(sequence_lengths) != 1:
        logging.error("Error: Not all sequences are of the same length.")
        return
    sequence_length = sequence_lengths.pop()
    logging.info(f"Sequence length: {sequence_length}")

    # Map sequence IDs to indices
    seq_dict = {seq.id: idx for idx, seq in enumerate(seqs)}
    logging.info(f"Sequence dictionary created with {len(seq_dict)} entries.")

    # Get edges and edge lengths from the tree
    edges, edge_lengths = get_edges_and_lengths(tree, seq_dict)
    num_edges = len(edges)
    if num_edges == 0:
        logging.error("No valid edges found in the tree. Exiting.")
        return
    logging.info(f"Number of edges: {num_edges}")
    logging.info(f"Maximum branch length: {np.max(edge_lengths)}")

    # Paths for mutation files
    mut_file = os.path.join(output_dir, 'mut.npy')
    mut_types_file = os.path.join(output_dir, 'mut_type.npy')

    # Convert sequences to a NumPy array
    logging.info("Converting sequences to NumPy array...")
    try:
        seq_array = np.array([list(str(seq.seq)) for seq in seqs])  # Shape: (num_sequences, sequence_length)
    except Exception as e:
        logging.error(f"Error converting sequences to NumPy array: {e}")
        return

    # Load or compute mutation matrices
    if os.path.exists(mut_file) and os.path.exists(mut_types_file):
        logging.info(f"Loading mutation matrices from '{mut_file}' and '{mut_types_file}'...")
        try:
            muts = np.load(mut_file)           # Shape: (sequence_length, num_edges)
            mutsTyp = np.load(mut_types_file)  # Shape: (sequence_length, num_edges)
            logging.info("Mutation matrices loaded successfully.")
        except Exception as e:
            logging.error(f"Error loading mutation matrices: {e}")
            return
    else:
        logging.info("Computing mutation matrices...")
        try:
            muts, mutsTyp = compute_mutation_matrices(seq_array, edges)  # Shape: (sequence_length, num_edges)
            # Save mutation matrices
            os.makedirs(output_dir, exist_ok=True)
            np.save(mut_file, muts)
            np.save(mut_types_file, mutsTyp)
            logging.info(f"Mutation matrices computed and saved to '{mut_file}' and '{mut_types_file}'.")
        except Exception as e:
            logging.error(f"Error computing mutation matrices: {e}")
            return

    # Create sequence presence matrix: 1 for 'C', 0 otherwise
    logging.info("Creating sequence presence matrix...")
    try:
        seq_presence = (seq_array == 'C').astype(int)  # Shape: (num_sequences, sequence_length)
    except Exception as e:
        logging.error(f"Error creating sequence presence matrix: {e}")
        return

    # Compute fraction of 'C's at each site
    logging.info("Calculating fraction of 'C's at each site...")
    try:
        frac_present = np.mean(seq_presence, axis=0)  # Shape: (sequence_length,)
    except Exception as e:
        logging.error(f"Error computing fraction of 'C's: {e}")
        return

    # Identify positions with at least 5% presence
    positions_to_consider = np.where(frac_present >= 0.05)[0]
    num_positions = len(positions_to_consider)
    logging.info(f"Number of positions with at least 5% 'C' presence: {num_positions}")

    if num_positions == 0:
        logging.error("No positions meet the 5% presence threshold. Exiting.")
        return

    # Filter mutations to include only positions to consider
    muts_filtered = muts[positions_to_consider, :]        # Shape: (num_positions, num_edges)
    mutsTyp_filtered = mutsTyp[positions_to_consider, :]  # Shape: (num_positions, num_edges)

    # Prepare the mapping matrix for likelihood calculations
    map_matrix = np.array([
        [0, 2, 8, 10],
        [1, 3, 9, 11],
        [4, 6, 12, 14],
        [5, 7, 13, 15]
    ])  # Shape: (4, 4)

    # Define parameter bounds for optimization based on maximum branch length
    max_branch_length = np.max(edge_lengths)
    scaling_factor = 0.5
    upper_boundary = scaling_factor * (1.0 / max_branch_length)
    logging.info(f"Setting upper boundary for substitution rates to {upper_boundary:.6f}")

    global_param_bounds = [
        (1e-5, upper_boundary),  # a1
        (1e-5, upper_boundary),  # a2
        (1e-5, upper_boundary),  # b1
        (1e-5, upper_boundary),  # b2
        (1.0, 100.0)             # eps (allowing flexibility)
    ]

    # Create a mapping from position to index for efficient lookups
    pos_map = {pos: idx for idx, pos in enumerate(positions_to_consider)}

    logging.info("Preparing position pairs for analysis...")

    # Initialize variables for chunk processing
    position_pairs = []
    chunk_counter = 1
    total_pairs_processed = 0

    # Initialize ProcessPoolExecutor once with shared data
    logging.info(f"Starting parallel processing with {num_threads} threads...")
    with ProcessPoolExecutor(max_workers=num_threads, initializer=initializer,
                             initargs=(muts_filtered, mutsTyp_filtered, edge_lengths, map_matrix, global_param_bounds)) as executor:
        # Add progress bar to the position pair preparation loop
        with tqdm(total=num_positions, desc="Collecting Position Pairs") as pbar:
            for idx_i, pos_i in enumerate(positions_to_consider):
                muts_i = muts_filtered[idx_i, :]
                mut_edges_i = np.where(muts_i)[0]

                if mut_edges_i.size == 0:
                    pbar.update(1)
                    continue  # Skip if no mutations at position i

                muts_subset = muts_filtered[idx_i:, mut_edges_i]
                col_sums = np.sum(muts_subset, axis=1)

                w_indices = np.where(col_sums >= 4)[0] + idx_i
                w_indices = w_indices[w_indices != idx_i]

                for idx_j in w_indices:
                    pos_j = positions_to_consider[idx_j]
                    # Only pass minimal data: indices
                    args = (idx_i, pos_i, idx_j, pos_j)
                    position_pairs.append(args)
                    total_pairs_processed += 1

                    # Check if we've reached the chunk size
                    if len(position_pairs) >= 10000:  # Adjust chunk size as needed
                        # Check if the result file for this chunk already exists
                        chunk_file = get_chunk_file_name(output_dir, chunk_counter)
                        if os.path.exists(chunk_file):
                            logging.info(f"Chunk {chunk_counter} already processed. Skipping.")
                        else:
                            logging.info(f"Processing chunk {chunk_counter} with {len(position_pairs)} position pairs...")
                            process_chunk(position_pairs, chunk_counter, executor, output_dir)
                            logging.info(f"Chunk {chunk_counter} processed.")
                        position_pairs = []
                        chunk_counter += 1

                pbar.update(1)

    # After iterating through all positions, process any remaining position_pairs
    if position_pairs:
        chunk_file = get_chunk_file_name(output_dir, chunk_counter)
        if os.path.exists(chunk_file):
            logging.info(f"Chunk {chunk_counter} already processed. Skipping.")
        else:
            logging.info(f"Processing final chunk {chunk_counter} with {len(position_pairs)} position pairs...")
            process_chunk(position_pairs, chunk_counter, executor, output_dir)
            logging.info(f"Chunk {chunk_counter} processed.")

    # Combine all chunk result files into a single file
    combine_chunk_results(output_dir, chunk_counter)

    logging.info("Process completed successfully.")

def process_chunk(position_pairs, chunk_id, executor, output_dir):
    """
    Process a single chunk of position pairs using the provided executor.
    """
    results = []
    num_pairs = len(position_pairs)
    logging.info(f"Processing chunk {chunk_id} with {num_pairs} position pairs...")

    # Use tqdm to display progress bar over the completed tasks
    for res in tqdm(executor.map(process_position_pair, position_pairs),
                    total=num_pairs, desc=f"Processing Chunk {chunk_id}"):
        if res is not None:
            results.append(res)

    # Sort results by Fisher's score and then by LRT score
    results.sort(key=lambda x: (x[6], x[7]))  # Sorting by score_fisher, then score_lrt

    # Save results to a TSV file for this chunk
    output_file = get_chunk_file_name(output_dir, chunk_id)
    try:
        with open(output_file, 'w', newline='') as fout:
            writer = csv.writer(fout, delimiter='\t')
            writer.writerow(['Position_i', 'Position_j', 'Tab_00', 'Tab_01', 'Tab_10', 'Tab_11', 'Score_Fisher', 'Score_LRT'])
            writer.writerows(results)
        logging.info(f"Chunk {chunk_id} results saved to '{output_file}'.")
    except Exception as e:
        logging.error(f"Error saving chunk {chunk_id} results to '{output_file}': {e}")

def combine_chunk_results(output_dir, max_chunk_id):
    """
    Combine all chunk result files into a single file.
    """
    combined_results_file = os.path.join(output_dir, 'score_combined.tsv')
    with open(combined_results_file, 'w', newline='') as fout_combined:
        writer = csv.writer(fout_combined, delimiter='\t')
        writer.writerow(['Position_i', 'Position_j', 'Tab_00', 'Tab_01', 'Tab_10', 'Tab_11', 'Score_Fisher', 'Score_LRT'])
        for chunk_id in range(1, max_chunk_id + 1):
            chunk_file = get_chunk_file_name(output_dir, chunk_id)
            if os.path.exists(chunk_file):
                with open(chunk_file, 'r') as fin_chunk:
                    reader = csv.reader(fin_chunk, delimiter='\t')
                    next(reader)  # Skip header
                    writer.writerows(reader)
                logging.info(f"Chunk {chunk_id} results appended to '{combined_results_file}'.")
            else:
                logging.warning(f"Chunk file '{chunk_file}' not found.")

if __name__ == '__main__':
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,  # Set to DEBUG for more detailed output
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler()
        ]
    )

    parser = argparse.ArgumentParser(description='Run CTMC analysis with chunk processing and checkpointing.')
    parser.add_argument('-t', '--treefile', required=True, help='Path to the Newick tree file.')
    parser.add_argument('-s', '--sequencefile', required=True, help='Path to the sequence file in FASTA format.')
    parser.add_argument('-o', '--output_dir', required=True, help='Directory to save the output files.')
    parser.add_argument('-n', '--threads', type=int, default=1, help='Number of threads (processes) to use.')

    args = parser.parse_args()
    main(args.treefile, args.sequencefile, args.output_dir, args.threads)
