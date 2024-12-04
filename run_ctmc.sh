#!/bin/bash

#SBATCH --job-name=ctmc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --partition=medium
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kell7366@ox.ac.uk
#SBATCH --error=/data/biol-micro-genomics/kell7366/kmer_gwes/ctmc_test_error.txt
#SBATCH --output=/data/biol-micro-genomics/kell7366/kmer_gwes/ctmc_test_output.txt

# Load necessary modules
module purge
module load Anaconda3/2024.02-1
source activate $DATA/kmer_gwes_env

tree_file="/data/biol-micro-genomics/kell7366/kmer_gwes/cfml_result/saureus.output.chrom_tree.newick"
new_tree_file="/data/biol-micro-genomics/kell7366/kmer_gwes/cfml_result/saureus.output.chrom_tree.copy.newick"
sequence_file="/data/biol-micro-genomics/kell7366/kmer_gwes/sccmec.fasta"
output_dir="/data/biol-micro-genomics/kell7366/kmer_gwes/ctmc_result_test"
temp_file="/data/biol-micro-genomics/kell7366/kmer_gwes/ctmc_result"
thread_num="12"

mkdir -p "$output_dir"

#modeltest-ng -i "$sequence_file" -d nt
#
## Generate the outgroup sequence
#new_sequence=$(printf 'a%.0s' $(seq 1 519347))
#new_sequence_file="/data/biol-micro-genomics/kell7366/kmer_gwes/sccmec.fasta.new"

#cp "$sequence_file" "$new_sequence_file"
#
## Append the outgroup to the sccmec.fasta file
#echo -e ">outgroup\n$new_sequence" >> "$new_sequence_file"
#
#log_file="${sequence_file}.log"
#
##echo "$log_file"
##
### Extract the line containing '> iqtree'
##iqtree_command=$(grep "> iqtree" "$log_file" | head -1 | sed -e "s|> iqtree|iqtree|" -e "s|\$| -te $tree_file -asr -T $thread_num|")
##
##$iqtree_command
#
#iqtree -s "$new_sequence_file" -m TVM+G4 -te "$new_tree_file" -asr 
#
#
## Define input and output files
#state_file="${new_sequence_file}.state"
#output_fasta="output_multifasta_file.fasta"
#
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
#
#cat "$sequence_file" "$output_fasta" > "${sequence_file}.combined"

# Convert sequences to uppercase
#awk '/^>/ {print $0; next} {print toupper($0)}' "${sequence_file}.combined" > "${sequence_file}_upper.combined"

#run python code
echo "Number of threads: $thread_num"

#python ctmc.py -t "$new_tree_file" -s "${sequence_file}_upper.combined" -o "$output_dir"
stdbuf -oL -eL python ctmc_multi.py -t "$new_tree_file" -s "${sequence_file}_upper.combined" -o "$output_dir" -n 36