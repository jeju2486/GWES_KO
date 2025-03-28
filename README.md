# GWES with Evolutionary DATA

This is add-on tool to PAN-GWES by using Continuous Time Markov Chain (CTMC) model.
This is the alignment-free pangenome-wise Epistasis Anaylsis which includes the phylogenetic data.

Please fine the individual repo for more details [PAN-GWES](https://github.com/Sudaraka88/PAN-GWES?tab=readme-ov-file) and original [CTMC-model](https://github.com/xavierdidelot/campy/tree/main)

## Requirements

Before using `maskGWAS`, ensure the following software are installed:

1. **[PAN-GWES](https://github.com/Sudaraka88/PAN-GWES?tab=readme-ov-file)** meat-package which includes
    * **[SpydrPick](https://github.com/santeripuranen/SpydrPick)** (version 1.2.0 or higher)
    * **[cuttlefish](https://www.cog-genomics.org/plink/1.9/)** (version 2.2.0 or higher)
    * **[gfa1_parser](https://github.com/jurikuronen/PANGWES/tree/main/gfa1_parser)** 
2. **[IQ-TREE](https://github.com/iqtree/iqtree2)** (version 2.3.5 or higher)

#Can you please make similar as above from below? skip the unimportant parts

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

## Installation

Clone the repository to your local directory:

```ruby
git clone https://github.com/jeju2486/GWES_KO
cd GWES_KO
```
## Input

To run this analysis pipeline you will need:
- tree file of your isolate
- folder that contains all your assemblies fasta files

## Usage

### 1. Build the pangenome graph (coloured de Bruijn graph) 

```ruby
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
```

### 2. Parse the gfa1-formatted output files

```ruby
#This code will parse the gfa1 file generate
# --cdg_paths: path info of de brujin graph
# cdbg.counts: count information
# cdbg.edges: edge info
# cdbg.fasta: kmer presence(C)/absence(A) fasta file 
# cdbg.paths: path file info
# cdbg.unitigs: uniigs sequence info (0-based)
gfa1_parser cdbg.gfa1 cdbg
```

### 3. build the inter node sequences from tree file

```ruby
# Input files
sequence_file="cdbg.fasta"

# Choose an appropriate substitution model for your data, for example GTR+I+G for DNA.
MODEL="GTR+I+G"

# Prefix for output files
output_prefix="output_asr"
tree_file="renamed_tree.nwk"

# Run IQ-tree with the ancestral state reconstruction option (-asr)
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
```

### 4. Run the Spydrpick and Unitig distance calculation

```ruby
sequence_file="cdbg.fasta"
spydrpick_output="./efcls_spydrpick_output"
output_dir="./ctmc_efcls_result"
tree_file="${spydrpick_output}/output_asr.treefile"

mkdir -p "$output_dir" "$spydrpick_output"

awk '/^>/{print} !/^>/{print toupper($0)}' "$spydrpick_output"/cdbg.fasta > "$spydrpick_output"/cdbg_upper.fasta

cat "$spydrpick_output"/cdbg_upper.fasta "$spydrpick_output"/output_asr.fasta > "$output_dir"/output_asr_combined.fasta
#output format: pos1 pos2 genome_distance ARACNE MI
SpydrPick --alignmentfile "$output_dir"/output_asr_combined.fasta --maf-threshold 0.05 --mi-values 50000000 --threads $SLURM_CPUS_PER_TASK --verbose

##output format v w distance (flag) (score) (count) (M2) (min_distance) (max_distance)
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
```

### 5. Run the CTMC

```ruby
python ctmc_multi.py -t "$tree_file" \
                    -s "$output_dir"/output_asr_combined.fasta \
                    -p "$spydrpick_output"/cdbg.ud_sgg_0_based \
                    -o "$output_dir" \
                    -n $SLURM_CPUS_PER_TASK 
```

### 6. Plot the data

```ruby
Rscript gwes_plot.r cdbg.ud_sgg_0_based
```