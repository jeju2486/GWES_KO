# GWES with Evolutionary DATA

This repository provides an add-on to **[PAN-GWES](https://github.com/Sudaraka88/PAN-GWES?tab=readme-ov-file)**, incorporating a Continuous Time Markov Chain (CTMC) model for alignment-free, pangenome-wide epistasis analysis that integrates phylogenetic data. It references the original [CTMC-model](https://github.com/xavierdidelot/campy/tree/main) by Xavier Didelot.

## Table of Contents

1. [Overview](#overview)  
2. [Requirements](#requirements)  
3. [Installation](#installation)  
4. [Usage](#usage)  
   1. [Build the pangenome graph](#1-build-the-pangenome-graph)  
   2. [Parse the GFA1 output files](#2-parse-the-gfa1-output-files)  
   3. [ASR with IQ-TREE](#3-ancestral-state-reconstruction-asr-with-iq-tree)  
   4. [Run SpydrPick and Unitig distance](#4-run-spydrpick-and-unitig-distance-calculations)  
   5. [Run the CTMC pipeline](#5-run-the-ctmc-pipeline)  
   6. [Plot the data](#6-plot-the-data)  
5. [Acknowledgments](#acknowledgments)  
6. [Contact](#contact)


## Overview

**GWES with Evolutionary DATA** is designed to:

- Perform alignment-free pangenome-wide epistasis analysis.  
- Integrate phylogenetic data using a CTMC model.  
- Leverage the [PAN-GWES](https://github.com/Sudaraka88/PAN-GWES?tab=readme-ov-file) framework with additional steps for reconstructing ancestral states and applying evolutionary models.


## Requirements

Before running this workflow, ensure you have installed or can access:

1. **[PAN-GWES](https://github.com/Sudaraka88/PAN-GWES?tab=readme-ov-file)**
   - **[SpydrPick](https://github.com/santeripuranen/SpydrPick)** (≥ v1.2.0)  
   - **[cuttlefish](https://github.com/COMBINE-lab/cuttlefish)** (≥ v2.2.0)  
   - **[gfa1_parser](https://github.com/jurikuronen/PANGWES/tree/main/gfa1_parser)**  
2. **[IQ-TREE](https://github.com/iqtree/iqtree2)** (≥ v2.3.5)  
3. **Python** (≥ 3.7), including packages:
   - `numpy`, `scipy`, `biopython`, `matplotlib`, `tqdm`, `numba`

> **Note:**  
> - The instructions below assume a Linux environment.  
> - Future releases may include a conda environment for simpler installation.

## Installation

Clone this repository into your local environment:

```bash
git clone https://github.com/jeju2486/GWES_KO
cd GWES_KO
```

## Usage

Below is a step-by-step example of how to run the pipeline. Adapt file paths, filenames, and parameters to your own data.  
In these examples, we use **`cdbg`** as the output prefix, but you can choose another as needed.

### 1. Build the pangenome graph

First, list all assemblies in a file:

```bash
ls -d "$input_dir"/* > assemblies_list.txt
```

Then build the coloured de Bruijn graph (CDBG) with a chosen k-mer length (e.g., `61`):

```bash
cuttlefish build --list assemblies_list.txt \
                 --kmer-len 61 \
                 --output cdbg \
                 --threads $SLURM_CPUS_PER_TASK \
                 -f 1
```

This produces a GFA1-formatted pangenome graph: `cdbg.gfa1`.

### 2. Parse the GFA1 output files

Use **`gfa1_parser`** to generate additional outputs for readability:

```bash
gfa1_parser cdbg.gfa1 cdbg
```

This command produces:

- **`cdbg.counts`** – Count information  
- **`cdbg.edges`** – Edge data of the pangenome graph  
- **`cdbg.fasta`** – Kmer presence(C)/absence(A) binary FASTA file  
- **`cdbg.paths`** – Path file info  
- **`cdbg.unitigs`** – Unitig sequence info (0-based)

### 3. Ancestral State Reconstruction (ASR) with IQ-TREE

Run **IQ-TREE** to reconstruct ancestral states. For more models and options, see the [IQ-TREE documentation](https://github.com/iqtree/iqtree2).

```bash
iqtree -s cdbg.fasta \
       -te ${tree_file} \
       -asr \
       -st DNA \
       -m GTR+I+G \
       -pre output_asr \
       -nt $SLURM_CPUS_PER_TASK
```

The file `output_asr.state` shows state probabilities for each tree node but is not in FASTA format. Convert it to a multi-FASTA:

```bash
awk -v output_fasta="output_asr.fasta" '
BEGIN { FS="\t"; }
/^#/ || $1 == "Node" { next }
{
   if ($1 != current_node) {
      if (current_node != "") print sequence > output_fasta;
      current_node = $1;
      print ">" current_node > output_fasta;
      sequence = "";
   }
   sequence = sequence $3;
}
END { if (current_node != "") print sequence > output_fasta; }
' "output_asr.state"
```

Finally, normalize/uppercase and combine ancestral sequences:

```bash
awk '/^>/{print} !/^>/{print toupper($0)}' cdbg.fasta > cdbg_upper.fasta
cat cdbg_upper.fasta output_asr.fasta > output_asr_combined.fasta
```

### 4. Run SpydrPick and Unitig distance calculations

Identify couplings (e.g., correlated variations) with **SpydrPick** and compute distances with **unitig_distance**. The `k-mer-length` should match what you used in **cuttlefish**.

```bash
SpydrPick --alignmentfile output_asr_combined.fasta \
          --maf-threshold 0.05 \
          --mi-values 50000000 \
          --threads $SLURM_CPUS_PER_TASK \
          --verbose

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

The key output from `unitig_distance` is `cdbg.ud_sgg_0_based`, which typically has columns:

```
unitig_i  unitig_j  distance  ARACNE  MI  count  M2  min_distance  max_distance
```

### 5. Run the CTMC pipeline

Apply the CTMC model to the results:

```bash
python ctmc_multi.py -t "$tree_file" \
                     -s output_asr_combined.fasta \
                     -p cdbg.ud_sgg_0_based \
                     -o "$output_dir" \
                     -n $SLURM_CPUS_PER_TASK
```

The resultant output typically includes:

```
unitig_i  unitig_j  distance  ARACNE  MI  count  M2  ...
...       LRT_p-value  -log_10_LRT_p-value
```

### 6. Plot the data (Optional)

Use your preferred plotting environment. For example, in R:

```bash
Rscript gwes_plot.r cdbg.ud_sgg_0_based
```

This script interprets `ARACNE == 1` as a direct correlation and `ARACNE == 0` as indirect.

## Acknowledgments

- **[Original CTMC code](https://github.com/xavierdidelot/campy/tree/main)** by Xavier Didelot  
- **[PAN-GWES](https://github.com/Sudaraka88/PAN-GWES?tab=readme-ov-file)** for pangenome-wide epistasis analysis


## Contact

For questions, bug reports, or troubleshooting, please open an issue in this repository or contact the repository owner.