# GWES with Evolutionary DATA

This is the alignment-free pangenome-wise Epistasis Anaylsis which includes the phylogenetic data.

# Input file
*`unitig.fasta` : output file from the cuttlefish/gfa1_parser. It is the mock fasta file that shows the presence of unitig as 'c' and absence as 'a'.
*`/unitig_path` : folder that contains the path data of pangenome graph.
*`name.unitigs` : output ifle from the cuttlefish/gfa1_parser. shows the name/order of unitig and their sequence.
*`treefile.nwk` : phylogenetic tree file. newick format

# Introduction
1. It reads the input files and apply the treefile to generate the middle node tree sequence.
2. Read the tree file and allocate unitig presence/absence for each node and leaf of three.
3. It calcaulte the model performance of null and epistasis CTMC model and perform the Likelihood-ratio test (LRT)
4. ##Todo-list Build the accurate way to calculate the distance and plot it with LRT result

# How to run
`python ctmc_multi.py -t "tree_file" -s "sequence_file" -o "$output_dir" -n 36`
(**warning** it is not runable generally now, I will change the code in the earliest convenience)
