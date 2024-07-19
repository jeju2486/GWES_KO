#!/bin/bash

#SBATCH --job-name=ctmc_run
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=medium
#SBATCH --cluster=htc
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=kell7366@ox.ac.uk
#SBATCH --error=/data/biol-micro-genomics/kell7366/kmer_gwes/ctmc_error.txt
#SBATCH --output=/data/biol-micro-genomics/kell7366/kmer_gwes/ctmc_output.txt

# Load necessary modules
module purge
module load Anaconda3/2024.02-1
source activate $DATA/kmer_gwes_env

tree_file="/data/biol-micro-genomics/kell7366/kmer_gwes/cfml_result/saureus.output.chrom_tree.newick"
sequence_file="/data/biol-micro-genomics/kell7366/kmer_gwes/sccmec.fasta"
output_dir="/data/biol-micro-genomics/kell7366/kmer_gwes/ctmc_result"
temp_file="/data/biol-micro-genomics/kell7366/kmer_gwes/ctmc_result"
thread_num="36"

mkdir -p "$output_dir"

#modeltest-ng -i "$sequence_file" -d nt

log_file="${sequence_file}.log"

#echo "$log_file"
#
## Extract the line containing '> iqtree'
#iqtree_command=$(grep "> iqtree" "$log_file" | head -1 | sed -e "s|> iqtree|iqtree|" -e "s|\$| -te $tree_file -asr -T $thread_num|")
#
#$iqtree_command

#iqtree -s "$sequence_file" -m TVM+G4 -te "$tree_file" -asr -T 36


# Define input and output files
state_file="${sequence_file}.state"
output_fasta="output_multifasta_file.fasta"

## Run awk to process the state file and generate the FASTA file
#awk '
#BEGIN {
#    FS="\t";
#}
#
#/^#/ { next }  # Skip comment lines
#
## Skip the header line
#NR == 1 { next }
#
#{
#    if ($1 != current_node) {
#        if (current_node != "") {
#            print sequence > output_fasta;
#        }
#        current_node = $1;
#        print ">" current_node > output_fasta;
#        sequence = "";
#    }
#    sequence = sequence $3;
#}
#
#END {
#    if (current_node != "") {
#        print sequence > output_fasta;
#    }
#}
#' output_fasta="$output_fasta" "$state_file"
#
#cat output_multifasta_file.fasta | awk '/^>/{if (seqname) print seqname "\t" length(seq); seqname=$1; seq=""; next} {seq = seq $0} END {print seqname "\t" length(seq)}' | sed 's/>//; s/ / /' > ./sccmec_inner_genome_size.txt

#cat "$sequence_file" "$output_fasta" > "${sequence_file}.combined"

#run python code
#python ctmc.py -t "$tree_file" -s "${sequence_file}.combined" -o "$output_dir"
python ctmc_multi.py -t "$tree_file" -s "${sequence_file}.combined" -o "$output_dir" -T 24