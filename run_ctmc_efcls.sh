#!/bin/bash

#SBATCH --job-name=ctmc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=short,medium,long
#SBATCH --cluster=all
#SBATCH --time=0-03:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kell7366@ox.ac.uk
#SBATCH --error=/data/biol-micro-genomics/kell7366/kmer_gwes/simulation_result/ctmc_ud_error.txt
#SBATCH --output=/data/biol-micro-genomics/kell7366/kmer_gwes/simulation_result/ctmc_ud_output.txt

# Load necessary modules
module purge
module load Anaconda3/2024.02-1
source activate $DATA/kmer_gwes_env

sequence_file="cdbg.fasta"
spydrpick_output="./efcls_spydrpick_output"
output_dir="./ctmc_efcls_result"
tree_file="${spydrpick_output}/output_asr.treefile"

mkdir -p "$output_dir" "$spydrpick_output"

awk '/^>/{print} !/^>/{print toupper($0)}' "$spydrpick_output"/cdbg.fasta > "$spydrpick_output"/cdbg_upper.fasta

cat "$spydrpick_output"/cdbg_upper.fasta "$spydrpick_output"/output_asr.fasta > "$output_dir"/output_asr_combined.fasta
#output format: pos1 pos2 genome_distance ARACNE MI
SpydrPick --alignmentfile "$output_dir"/output_asr_combined.fasta --maf-threshold 0.05 --mi-values 50000000 --threads $SLURM_CPUS_PER_TASK --verbose

##v w distance (flag) (score) (count) (M2) (min_distance) (max_distance)
unitig_distance --unitigs-file cdbg.unitigs \
                --edges-file cdbg.edges \
                --k-mer-length 61 \
                --sgg-paths-file cdbg.paths \
                --queries-file output_asr.*spydrpick_couplings*edges \
                --threads $SLURM_CPUS_PER_TASK \
                --queries-one-based \
                --run-sggs-only \
                --output-stem cdbg \
                --verbose

#Run CTMC
conda deactivate
source activate $DATA/python3_12_2
                
python ctmc_multi.py -t "$tree_file" \
                    -s "$output_dir"/output_asr_combined.fasta \
                    -p "$spydrpick_output"/cdbg.ud_sgg_0_based \
                    -o "$output_dir" \
                    -n $SLURM_CPUS_PER_TASK 

conda deactivate
module purge
module load R/4.4.1-gfbf-2023b      
         
Rscript gwes_plot.r cdbg.ud_sgg_0_based