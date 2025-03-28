#!/bin/bash

#SBATCH --job-name=iqtree
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=short,medium,long
#SBATCH --cluster=all
#SBATCH --time=0-12:00:00
#SBATCH --error=/data/biol-micro-genomics/kell7366/kmer_gwes/simulation_result/iqtree_error.txt
#SBATCH --output=/data/biol-micro-genomics/kell7366/kmer_gwes/simulation_result/iqtree_output.txt

# Load necessary modules
module purge
module load Anaconda3/2024.02-1
source activate $DATA/kmer_gwes_env

# Input files
sequence_file="cdbg.fasta"

# Choose an appropriate substitution model for your data, for example GTR+I+G for DNA.
MODEL="GTR+I+G"

# Prefix for output files
output_prefix="output_asr"
tree_file="renamed_tree.nwk"

# Run IQ-tree with the ancestral state reconstruction option (-asr)
#iqtree -s /data/biol-micro-genomics/kell7366/kmer_gwes/simulation_result/pirate_output/core_alignment.fasta -st DNA -m ${MODEL} -pre ${output_prefix} -bb 1000 -nt $SLURM_CPUS_PER_TASK 
#iqtree -s /data/biol-micro-genomics/kell7366/kmer_gwes/simulation_result/pirate_output/core_alignment.fasta -st DNA -m ${MODEL} -pre ${output_prefix} -bb 1000 -asr -nt $SLURM_CPUS_PER_TASK 
iqtree -s ${sequence_file} -te ${tree_file} -asr -st DNA -m ${MODEL} -pre ${output_prefix} -nt $SLURM_CPUS_PER_TASK 


# Define input and output files
state_file="${output_prefix}.state"

echo "Converting IQ tree outputfile to multifasta format"

# Run awk to process the state file and generate the FASTA file
awk -v output_fasta="${output_prefix}.fasta" '
BEGIN {
    FS="\t";
}

/^#/ { next }  # Skip comment lines

# Skip header line with "Node" in the first field
$1 == "Node" { next }

{
    if ($1 != current_node) {
        if (current_node != "") {
            print sequence > output_fasta;
        }
        current_node = $1;
        print ">" current_node > output_fasta;
        sequence = "";
    }
    sequence = sequence $3;
}

END {
    if (current_node != "") {
        print sequence > output_fasta;
    }
}
' "$state_file"

echo "Recording the genome sizes"

cat "${output_prefix}.fasta" | awk '/^>/{if (seqname) print seqname "\t" length(seq); seqname=$1; seq=""; next} {seq = seq $0} END {print seqname "\t" length(seq)}' | sed 's/>//; s/ / /' > "$output_prefix"_genome_size.txt

echo "Sorting the output files and removing the temp files"

cat "$sequence_file" "${output_prefix}.fasta" > "${output_prefix}_combined.fasta"

#python split_fasta.py -i "${output_prefix}.fasta" -f efcls_fasta

gzip "$state_file"

##rm ${output_prefix}.ckp.gz ${output_prefix}.state ${output_prefix}.fasta ${output_prefix}.log 
#mv ${output_prefix}.treefile ${output_prefix}_treefile.nwk