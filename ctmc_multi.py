#!/usr/bin/env python3

import argparse
import csv
import logging
import math
import os

import matplotlib.pyplot as plt  # If you need plotting later
import numpy as np
from Bio import Phylo, SeqIO
from concurrent.futures import ProcessPoolExecutor
from numba import njit
from scipy.optimize import minimize
from scipy.stats import chi2 as chi2_dist
from tqdm import tqdm

# ---------------------------
# Global variables for parallel LRT stage
# ---------------------------
global_muts = None         # Not strictly needed if you only use global_mutsTyp
global_mutsTyp = None
global_edge_lengths = None
global_map_matrix = None
global_param_bounds = None
global_seq_array = None
global_bra = None  # Index of the "root branch" in edges


def initializer(muts_input, mutsTyp_input, edge_len, map_mat, param_bounds, seqarr):
    """
    Store large data structures in global variables so that process_position_pair
    can read them without capturing unpickleable local objects.
    """
    global global_muts, global_mutsTyp, global_edge_lengths
    global global_map_matrix, global_param_bounds, global_seq_array
    global_muts = muts_input
    global_mutsTyp = mutsTyp_input
    global_edge_lengths = edge_len
    global_map_matrix = map_mat
    global_param_bounds = param_bounds
    global_seq_array = seqarr


@njit
def Q_matrix(a1, a2, b1, b2, eps):
    """
    Build the 4x4 Q matrix given parameters.
    (A->A)=0, (C->A)=1, (A->C)=2, (C->C)=3 representation.
    """
    Q = np.empty((4, 4), dtype=np.float64)
    # Row 0
    Q[0, 0] = -(a1 + b1)
    Q[0, 1] = a1
    Q[0, 2] = b1
    Q[0, 3] = 0.0
    # Row 1
    Q[1, 0] = a2
    Q[1, 1] = -(a2 + b1 * eps)
    Q[1, 2] = 0.0
    Q[1, 3] = b1 * eps
    # Row 2
    Q[2, 0] = b2
    Q[2, 1] = 0.0
    Q[2, 2] = -(b2 + a1 * eps)
    Q[2, 3] = a1 * eps
    # Row 3
    Q[3, 0] = 0.0
    Q[3, 1] = b2 / eps
    Q[3, 2] = a2 / eps
    Q[3, 3] = -(b2 / eps + a2 / eps)
    return Q


@njit
def dQ_deps(a1, a2, b1, b2, eps):
    """
    Partial derivative of the Q matrix w.r.t. eps.
    """
    dQ = np.zeros((4, 4), dtype=np.float64)
    # Only entries depending on eps:
    dQ[1, 1] = -b1
    dQ[1, 3] = b1
    dQ[2, 2] = -a1
    dQ[2, 3] = a1
    dQ[3, 1] = -b2 / (eps * eps)
    dQ[3, 2] = -a2 / (eps * eps)
    dQ[3, 3] = (b2 + a2) / (eps * eps)
    return dQ


@njit
def compute_nll_and_grad_core_numba(params, site1, site2, root1, root2, edge_lengths, map_matrix):
    """
    Numba-compiled function that returns (negative log-likelihood, gradient).
    We have placeholders for partial derivatives wrt (a1,a2,b1,b2), but do compute wrt eps.
    """
    a1, a2, b1, b2, eps = params

    # Construct Q and its derivative w.r.t. eps
    Q = Q_matrix(a1, a2, b1, b2, eps)
    dQeps = dQ_deps(a1, a2, b1, b2, eps)

    # Compute Q^2 and Q^3
    Q2 = Q @ Q
    Q3 = Q2 @ Q

    # Compute derivatives for Q^2 and Q^3
    dQ2 = Q @ dQeps + dQeps @ Q
    dQ3 = Q2 @ dQeps + Q @ dQ2 + dQeps @ Q2

    # Flatten
    Qf = Q.ravel()
    Q2f = Q2.ravel()
    Q3f = Q3.ravel()
    dQf_eps = dQeps.ravel()
    dQ2f_eps = dQ2.ravel()
    dQ3f_eps = dQ3.ravel()

    ll_sum = 0.0
    dll_da1 = dll_da2 = dll_db1 = dll_db2 = dll_deps = 0.0
    n_edges = edge_lengths.shape[0]

    # Summation over edges
    for i in range(n_edges):
        idx = map_matrix[site1[i], site2[i]]
        t = edge_lengths[i]
        t2 = t * t
        t3 = t2 * t

        fi = 1.0 \
             + Qf[idx]   * t \
             + Q2f[idx]  * (t2 * 0.5) \
             + Q3f[idx]  * (t3 / 6.0) \
             + 1e-10
        if fi <= 0.0:
            return 1e10, np.zeros(5, dtype=np.float64)
        ll_sum += math.log(fi)

        # derivative wrt eps
        dfi_deps = (dQf_eps[idx] * t
                    + dQ2f_eps[idx] * (t2 * 0.5)
                    + dQ3f_eps[idx] * (t3 / 6.0))
        inv_fi = 1.0 / fi
        dll_deps += inv_fi * dfi_deps

        # placeholders for a1,a2,b1,b2
        dll_da1 += 0.0
        dll_da2 += 0.0
        dll_db1 += 0.0
        dll_db2 += 0.0

    # Root correction term: log(num/den)
    den = a1 * b1 * (eps ** 2) + a1 * b2 + a2 * b1 + a2 * b2
    if den <= 0.0:
        return 1e10, np.zeros(5, dtype=np.float64)

    # figure out the "root states"
    if root1 == 1 and root2 == 1:
        num = a2 * b2
        dnum_da1_ = 0.0
        dnum_da2_ = b2
        dnum_db1_ = 0.0
        dnum_db2_ = a2
        dnum_deps_ = 0.0
    elif root1 == 1 and root2 == 2:
        num = a2 * b1
        dnum_da1_ = 0.0
        dnum_da2_ = b1
        dnum_db1_ = a2
        dnum_db2_ = 0.0
        dnum_deps_ = 0.0
    elif root1 == 2 and root2 == 1:
        num = a1 * b2
        dnum_da1_ = b2
        dnum_da2_ = 0.0
        dnum_db1_ = 0.0
        dnum_db2_ = a1
        dnum_deps_ = 0.0
    elif root1 == 2 and root2 == 2:
        num = a1 * b1 * (eps ** 2)
        dnum_da1_ = b1 * (eps ** 2)
        dnum_da2_ = 0.0
        dnum_db1_ = a1 * (eps ** 2)
        dnum_db2_ = 0.0
        dnum_deps_ = 2.0 * a1 * b1 * eps
    else:
        # If for some reason out of these states, set num near-zero
        num = 1e-10
        dnum_da1_ = dnum_da2_ = dnum_db1_ = dnum_db2_ = dnum_deps_ = 0.0

    if num <= 0.0:
        return 1e10, np.zeros(5, dtype=np.float64)

    ll_sum += math.log(num / den)

    # partials for den
    dden_da1 = b1 * (eps ** 2) + b2
    dden_da2 = b1 + b2
    dden_db1 = a1 * (eps ** 2) + a2
    dden_db2 = a1 + a2
    dden_deps = 2.0 * a1 * b1 * eps

    inv_num = 1.0 / num
    inv_den = 1.0 / den

    dlog_num_da1 = inv_num * dnum_da1_
    dlog_num_da2 = inv_num * dnum_da2_
    dlog_num_db1 = inv_num * dnum_db1_
    dlog_num_db2 = inv_num * dnum_db2_
    dlog_num_deps = inv_num * dnum_deps_

    dlog_den_da1 = inv_den * dden_da1
    dlog_den_da2 = inv_den * dden_da2
    dlog_den_db1 = inv_den * dden_db1
    dlog_den_db2 = inv_den * dden_db2
    dlog_den_deps = inv_den * dden_deps

    dll_da1 += (dlog_num_da1 - dlog_den_da1)
    dll_da2 += (dlog_num_da2 - dlog_den_da2)
    dll_db1 += (dlog_num_db1 - dlog_den_db1)
    dll_db2 += (dlog_num_db2 - dlog_den_db2)
    dll_deps += (dlog_num_deps - dlog_den_deps)

    nll = -ll_sum
    grad = np.array([-dll_da1, -dll_da2, -dll_db1, -dll_db2, -dll_deps], dtype=np.float64)
    return nll, grad


def process_position_pair(args):
    """
    Each 'args' is:
      (v_idx, w_idx, extra_data_tuple)

    Where:
      - v_idx, w_idx are the 0-based positions in the alignment
      - extra_data_tuple holds the 7 columns: (distance, flag, score, count, M2, min_dist, max_dist).

    We do the LRT:
      - "null" model = fix eps=1
      - "alt" model  = eps free
    Then return [v, w, distance, flag, score, count, M2, min_dist, max_dist, LRT_p, log10_LRT_p]
    """
    v_idx, w_idx, extra_cols = args
    dist_val, flag_val, score_val, count_val, m2_val, min_d, max_d = extra_cols

    # site vectors from the global_mutsTyp array
    site1 = global_mutsTyp[v_idx, :]
    site2 = global_mutsTyp[w_idx, :]

    # figure out root states in that "root branch"
    bra = global_bra
    root1 = 1 if site1[bra] in [1,2] else 2
    root2 = 1 if site2[bra] in [1,2] else 2

    def likelihood_null_and_grad(param4):
        # fix eps=1.0
        full_5 = np.array([param4[0], param4[1], param4[2], param4[3], 1.0], dtype=np.float64)
        nll_val, grad_5 = compute_nll_and_grad_core_numba(
            full_5, site1, site2, root1, root2, global_edge_lengths, global_map_matrix
        )
        return nll_val, grad_5[:4]

    def likelihood_alt_and_grad(param5):
        return compute_nll_and_grad_core_numba(
            param5, site1, site2, root1, root2, global_edge_lengths, global_map_matrix
        )

    # Minimization
    try:
        # Null: 4 parameters (a1,a2,b1,b2), fix eps=1
        fit_null = minimize(
            fun=likelihood_null_and_grad,
            x0=[1.0, 1.0, 1.0, 1.0],
            method="L-BFGS-B",
            jac=True,
            bounds=[global_param_bounds[0]]*4,
        )
        if not fit_null.success:
            # fallback
            LRT_p = 1.0
            log10_LRT = 0.0
            return [v_idx, w_idx, dist_val, flag_val, score_val,
                    count_val, m2_val, min_d, max_d, LRT_p, log10_LRT]

        # Alternative: 5 parameters (a1,a2,b1,b2,eps)
        alt_init = list(fit_null.x) + [1.0]
        fit_alt = minimize(
            fun=likelihood_alt_and_grad,
            x0=alt_init,
            method="L-BFGS-B",
            jac=True,
            bounds=global_param_bounds,
        )
        if not fit_alt.success:
            LRT_p = 1.0
            log10_LRT = 0.0
            return [v_idx, w_idx, dist_val, flag_val, score_val,
                    count_val, m2_val, min_d, max_d, LRT_p, log10_LRT]

        # LRT statistic
        lrt_stat = 2.0 * (-fit_alt.fun + fit_null.fun)
        if lrt_stat < 0:
            lrt_stat = 0.0
        LRT_p = chi2_dist.sf(lrt_stat, df=1)
        if LRT_p < np.finfo(float).tiny:
            LRT_p = np.finfo(float).tiny
        log10_LRT = round(math.log10(LRT_p), 6)

    except Exception as e:
        logging.error(f"Error in LRT for positions {v_idx} {w_idx}: {e}")
        LRT_p = 1.0
        log10_LRT = 0.0

    return [v_idx, w_idx, dist_val, flag_val, score_val,
            count_val, m2_val, min_d, max_d, LRT_p, log10_LRT]


def get_edges_and_lengths(tree, seq_dict):
    edges = []
    edge_lengths = []
    for clade in tree.find_clades(order="level"):
        for child in clade.clades:
            if clade.name in seq_dict and child.name in seq_dict:
                parent_idx = seq_dict[clade.name]
                child_idx = seq_dict[child.name]
                edges.append((parent_idx, child_idx))
                branch_length = child.branch_length if child.branch_length and child.branch_length > 0 else 1e-7
                edge_lengths.append(branch_length)
    return np.array(edges, dtype=int), np.array(edge_lengths, dtype=float)


def compute_mutation_matrices(seq_array, edges):
    # parent->child
    parent_seqs = seq_array[edges[:, 0], :]
    child_seqs = seq_array[edges[:, 1], :]
    # transpose so shape is [n_positions, n_edges]
    parent_seqs = parent_seqs.T
    child_seqs = child_seqs.T

    # mismatch matrix (not used here, but computed anyway)
    muts = (parent_seqs != child_seqs)

    # 4-state matrix (A->A=0, C->A=1, A->C=2, C->C=3)
    mutsTyp = np.zeros_like(parent_seqs, dtype=int)
    # fill
    mutsTyp[(parent_seqs == "A") & (child_seqs == "A")] = 0
    mutsTyp[(parent_seqs == "C") & (child_seqs == "A")] = 1
    mutsTyp[(parent_seqs == "A") & (child_seqs == "C")] = 2
    mutsTyp[(parent_seqs == "C") & (child_seqs == "C")] = 3

    return muts, mutsTyp


def get_chunk_file_name(output_dir, chunk_id):
    return os.path.join(output_dir, f"score_chunk_{chunk_id}.tsv")


def process_chunk(pairs_for_chunk, chunk_id, executor, output_dir):
    """
    pairs_for_chunk is a list of (v_idx, w_idx, extra_data_tuple).
    We'll run process_position_pair in parallel.
    """
    results = []
    num_pairs = len(pairs_for_chunk)
    logging.info(f"Processing chunk {chunk_id} with {num_pairs} position pairs...")

    for res in tqdm(
        executor.map(process_position_pair, pairs_for_chunk),
        total=num_pairs,
        desc=f"Chunk {chunk_id}",
    ):
        if res is not None:
            results.append(res)

    # Sort if you like. For example, by p-value
    # results.sort(key=lambda x: x[-2]) # LRT_p is second-last

    out_file = get_chunk_file_name(output_dir, chunk_id)
    try:
        with open(out_file, "w", newline="") as fout:
            writer = csv.writer(fout, delimiter="\t")
            # Write header
            writer.writerow([
                "v", "w", "distance", "flag", "score", "count",
                "M2", "min_distance", "max_distance", "LRT_p", "log10_LRT_p"
            ])
            writer.writerows(results)
        logging.info(f"Chunk {chunk_id} saved to '{out_file}'.")
    except Exception as e:
        logging.error(f"Error saving chunk {chunk_id} results: {e}")


def combine_chunk_results(output_dir, max_chunk_id):
    combined_path = os.path.join(output_dir, "score_combined.tsv")
    with open(combined_path, "w", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow([
            "v", "w", "distance", "flag", "score", "count",
            "M2", "min_distance", "max_distance", "LRT_p", "log10_LRT_p"
        ])
        for chunk_id in range(1, max_chunk_id + 1):
            chunk_file = get_chunk_file_name(output_dir, chunk_id)
            if os.path.exists(chunk_file):
                with open(chunk_file, "r") as fin:
                    reader = csv.reader(fin, delimiter="\t")
                    next(reader)  # skip header
                    writer.writerows(reader)
    logging.info(f"Combined results saved to '{combined_path}'.")


def main(treefile, sequencefile, output_dir, num_threads, pairs_file):
    """
    1) Read the tree, build the edges.
    2) Read the alignment, build the (muts, mutsTyp) arrays.
    3) Read your pairs from the `pairs_file`, which has 9 columns:
       v, w, distance, flag, score, count, M2, min_distance, max_distance
       All are 0-based for v,w.

    4) For each pair, do the LRT and produce:
       v, w, distance, flag, score, count, M2, min_d, max_d, LRT_p, log10_LRT_p
    """

    logging.info("Reading tree & sequences...")
    try:
        tree = Phylo.read(treefile, "newick")
    except Exception as e:
        logging.error(f"Error reading tree '{treefile}': {e}")
        return
    try:
        seqs = list(SeqIO.parse(sequencefile, "fasta"))
    except Exception as e:
        logging.error(f"Error reading seqs '{sequencefile}': {e}")
        return

    num_seqs = len(seqs)
    logging.info(f"# of sequences = {num_seqs}")
    lens = set(len(s.seq) for s in seqs)
    if len(lens) != 1:
        logging.error("Not all sequences have the same length.")
        return
    seq_len = lens.pop()
    logging.info(f"Alignment length = {seq_len}")

    # Build seq_dict
    seq_dict = {seq.id: i for i, seq in enumerate(seqs)}

    edges, edge_lengths = get_edges_and_lengths(tree, seq_dict)
    logging.info(f"# edges = {len(edges)}. Max branch length = {np.max(edge_lengths):.6f}")

    # Identify the root "branch"
    all_parents = np.unique(edges[:, 0])
    all_children = np.unique(edges[:, 1])
    root_nodes = np.setdiff1d(all_parents, all_children)
    if len(root_nodes) == 0:
        logging.warning("No root node found, continuing anyway.")
    else:
        root = root_nodes[0]
        bra_indices = np.where(edges[:, 0] == root)[0]
        if bra_indices.size > 0:
            global global_bra
            global_bra = bra_indices[0]
            logging.info(f"Root branch index = {global_bra}")
        else:
            logging.warning("No child from root node found.")

    # Convert to numpy
    seq_array = np.array([list(str(s.seq)) for s in seqs], dtype="<U1")

    # Precompute Muts, MutsTyp
    mut_file = os.path.join(output_dir, "mut.npy")
    mut_type_file = os.path.join(output_dir, "mut_type.npy")

    if os.path.exists(mut_file) and os.path.exists(mut_type_file):
        logging.info("Loading mutation matrices from disk...")
        try:
            muts = np.load(mut_file)
            mutsTyp = np.load(mut_type_file)
        except Exception as e:
            logging.error(f"Error loading precomputed arrays: {e}")
            return
    else:
        logging.info("Computing mutation matrices for all positions...")
        os.makedirs(output_dir, exist_ok=True)
        try:
            muts, mutsTyp = compute_mutation_matrices(seq_array, edges)
            np.save(mut_file, muts)
            np.save(mut_type_file, mutsTyp)
            logging.info("Saved Muts/MutsTyp to disk.")
        except Exception as e:
            logging.error(f"Error computing mutation matrices: {e}")
            return

    # Read pairs (v,w, distance, flag, score, count, M2, min_d, max_d)
    logging.info(f"Reading pairs from '{pairs_file}'...")
    pairs_raw = []
    with open(pairs_file, "r") as fin:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 9:
                continue  # or raise error
            v_idx = int(fields[0])  # 0-based
            w_idx = int(fields[1])  # 0-based
            dist_val = fields[2]
            flag_val = fields[3]
            score_val = fields[4]
            count_val = fields[5]
            m2_val = fields[6]
            min_d = fields[7]
            max_d = fields[8]

            # convert them as needed if they are integers/floats
            # e.g. dist_val = float(fields[2]) if you want
            # For now we keep them as strings or parse some as needed:
            # dist_val = float(fields[2])  # example
            pairs_raw.append((
                v_idx,
                w_idx,
                (dist_val, flag_val, score_val, count_val, m2_val, min_d, max_d)
            ))

    logging.info(f"Loaded {len(pairs_raw)} pairs.")

    # Build map_matrix
    # This matches the state coding: 4 states => 16 possible transitions
    map_matrix = np.array([
        [0,  2,  8,  10],  # parent=0 => child=0..3
        [1,  3,  9,  11],  # parent=1 => child=0..3
        [4,  6,  12, 14],  # parent=2 => child=0..3
        [5,  7,  13, 15],  # parent=3 => child=0..3
    ])

    # Bounds
    min_branch_length = np.min(edge_lengths)
    scaling_factor = 0.5
    upper_boundary = scaling_factor * (1.0 / min_branch_length)
    logging.info(f"Setting max bound for (a1,a2,b1,b2) = {upper_boundary:.6f}")

    global global_param_bounds
    global_param_bounds = [
        (1e-5, upper_boundary),  # a1
        (1e-5, upper_boundary),  # a2
        (1e-5, upper_boundary),  # b1
        (1e-5, upper_boundary),  # b2
        (1.0, 100.0)             # eps
    ]

    logging.info(f"Performing LRT on {len(pairs_raw)} pairs using {num_threads} threads...")

    # Parallel chunking
    lrt_chunk_size = 20000  # adjust as desired
    total_pairs = len(pairs_raw)
    chunk_counter = 1

    with ProcessPoolExecutor(
        max_workers=num_threads,
        initializer=initializer,
        initargs=(muts, mutsTyp, edge_lengths, map_matrix, global_param_bounds, seq_array),
    ) as executor:

        for i in range(0, total_pairs, lrt_chunk_size):
            chunk_pairs = pairs_raw[i : i + lrt_chunk_size]
            out_file = get_chunk_file_name(output_dir, chunk_counter)
            if os.path.exists(out_file):
                logging.info(f"Chunk {chunk_counter} already exists; skipping.")
                chunk_counter += 1
                continue

            process_chunk(chunk_pairs, chunk_counter, executor, output_dir)
            chunk_counter += 1

    # Combine results
    combine_chunk_results(output_dir, chunk_counter - 1)
    logging.info("All done.")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler()],
    )
    parser = argparse.ArgumentParser(
        description="CTMC-based LRT for pairs given in a 9-column file, appending (LRT_p, log10_LRT_p)."
    )
    parser.add_argument("-t", "--treefile", required=True, help="Path to Newick tree.")
    parser.add_argument("-s", "--sequencefile", required=True, help="Path to alignment in FASTA.")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for results.")
    parser.add_argument("-n", "--threads", type=int, default=1, help="Number of threads.")
    parser.add_argument("-p", "--pairs_file", required=True, help="File with 9 columns: v w distance flag score count M2 min_dist max_dist.")

    args = parser.parse_args()
    main(args.treefile, args.sequencefile, args.output_dir, args.threads, args.pairs_file)
