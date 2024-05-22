#!/bin/bash

#SBATCH --job-name=test_run
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=short
#SBATCH --time=00-01:00:00
#SBATCH --error=/data/biol-micro-genomics/kell7366/kmer_gwes/test_error.txt
#SBATCH --output=/data/biol-micro-genomics/kell7366/kmer_gwes/test_output.txt


##module loading
module purge
module load Anaconda3
source activate $DATA/newenv

##make the input.txt
#find /data/biol-micro-genomics/kell7366/minimap/ref_fas/ref_aureus -type f -name '*.fas' > input.txt
#shuf -n 100 input.txt > input_100.txt

###run cuttlefish ## 15mins
#cuttlefish build --list      input_100.txt \
#                 --kmer-len  61        \
#                 --output    sccmec     \
#                 --threads   8    \
#                 --format 1
#
## Parse the resulting Graphical Fragment Assembly (GFA) formatted file. ## 13mins
#./PANGWES/gfa1_parser/bin/gfa1_parser sccmec.gfa1 sccmec 

## Calculate mutual information scores and obtain the list of top candidate unitig pairs. ##2560.33s
#SpydrPick --alignmentfile  sccmec.fasta \
#            --maf-threshold  0.05           \
#            --mi-values      50000000       \
#            --threads        8              \
#            --verbose

### Calculate the graph distances between unitigs in the single genome subgraphs.
#./PANGWES/unitig_distance/bin/unitig_distance --unitigs-file    sccmec.unitigs                \
#                                              --edges-file      sccmec.edges                  \
#                                              --k-mer-length    61                              \
#                                              --sgg-paths-file  sccmec.paths                  \
#                                              --queries-file    sccmec.*spydrpick_couplings*edges \
#                                              --threads         8                               \
#                                              --queries-one-based                               \
#                                              --run-sggs-only                                   \
#                                              --output-stem sccmec                                \
#                                              --verbose

##module load
module purge
module load R/4.3.2-gfbf-2023a

## Draw a GWES Manhattan plot to visualize the results.

Rscript ./PANGWES/scripts/gwes_plot.r /data/biol-micro-genomics/kell7366/kmer_gwes/sccmec.ud_sgg_0_based /data/biol-micro-genomics/kell7366/kmer_gwes/sccmec_sgg_mean.png 0 "Mutual information" \
                    "" 0 "sccmec_sgg_mean (no sgg counts filtering)" \
