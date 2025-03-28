#!/bin/bash

#SBATCH --job-name=cuttlefish
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --partition=short
#SBATCH --cluster=all
#SBATCH --time=0-12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kell7366@ox.ac.uk
#SBATCH --error=/data/biol-micro-genomics/kell7366/kmer_gwes/simulation_result/ctmc_test_error.txt
#SBATCH --output=/data/biol-micro-genomics/kell7366/kmer_gwes/simulation_result/ctmc_test_output.txt

# Load necessary modules
module purge
module load Anaconda3/2024.02-1
source activate $DATA/kmer_gwes_env

# Prefix for output files
input_dir="efcls_assemblies"
input_prefix="output_asr"
output_prefix="cdbg"

ls -d "$input_dir"/* > efcls_assemblies.txt

cuttlefish build --list efcls_assemblies.txt \
                --kmer-len 61 \
                --output "${output_prefix}" \
                --threads $SLURM_CPUS_PER_TASK \
                -f 1
                
#This code will parse the gfa1 file generate
# --cdg_paths: path info of de brujin graph
# cdbg.counts: count information
# cdbg.edges: edge info
# cdbg.fasta: kmer presence(C)/absence(A) fasta file 
# cdbg.paths: path file info
# cdbg.unitigs: uniigs sequence info (0-based)
gfa1_parser cdbg.gfa1 cdbg

