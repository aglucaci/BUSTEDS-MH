# This repository contains software, data, and results for the BUSTED[S]-MH manuscript.

## About

A link for our preprint version will be provided here.

Conventional models may fail to accurately depict the evolutionary pressures and constraints at individual sites and along branches for a protein-coding sequence. [We have previously demonstrated widespread empirical support for models which account for instantaneous multiple nucleotide changes within a codon (MH)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0248337). Here, we present an extension of the commonly used branch-site framework of [BUSTED[S]](https://academic.oup.com/mbe/article/37/8/2430/5739973?login=false), a method for the inference of gene-wide episodic diversifying selection, developed in order to account for the impact of MH events (BUSTED[S]-MH). Our new model contributes towards an increased interest in the advances our understanding of biological realism within codon models and offers an attractive alternative to traditional models used for genome-wide scans of adaptive evolution. Comprehensive analysis of the effects of MH mutations which may act on a subset of branches in a phylogeny and at a subset of sites within the gene may shed light on the additional pathways available for genes to embark upon during the evolutionary process. 

## Requirements

1. `HyPhy`
2. `Snakemake`

## Installation (Conda-based)

1. `git clone https://github.com/aglucaci/BUSTEDS-MH.git`
2. `cd BUSTEDS-MH`
3. `conda env create -f environment.yml`.  
   This will create a conda environment called (BUSTEDS-MH) with the necessary dependencies.
4. At this point, run `conda activate BUSTEDS-MH` and your environment will be ready to go.
5. We also require the [HyPhy-analyses](https://github.com/veg/hyphy-analyses) repository for the BUSTED[S]-MH model batch file, you can clone this repository within the working direcory with `git clone https://github.com/veg/hyphy-analyses

## Datasets

We use a number of small empirical datasets, those available are available here, [Empirical 14 datasets](https://github.com/aglucaci/BUSTEDS-MH/tree/main/data/Empirical_14_datasets), and described in detail in the manuscript.

We use a number of large empirical datasets, those available are below, and described in detail in the manuscript.
<ol>
  <li>[Selectome, Euteleostomi, Version 6](https://github.com/aglucaci/BUSTEDS-MH/tree/develop/data/Empirical_Selectome)</li>
</ol>

We use a number of simulated datasets, those available are below, and described in detail in the manuscript.
<ol>
  <li>[Simulated 16/31 sequences, Power simulations](https://github.com/aglucaci/BUSTEDS-MH/tree/develop/data/SIMULATED_16_31)</li>
  <li>[Null Simulations](https://github.com/aglucaci/BUSTEDS-MH/tree/develop/data/SIMULATED_null)</li>
</ol>

## Minimum Configuration

Everything is controlled through the `cluster.json` and `config.yml` files.

The cluster JSON file contains information about where to submit your jobs on a HPC system.

The `config.yml` specifies the following variables:
`DataDirectory`, the full path to your data directory
`SubDirectory`, the subdirectory of that folder to analyze
`FileEnding`, the file ending to search for, for example `.fasta` or `.nex`
`OutputDirectory`, the ouput directory for your HyPhy JSON results
`BUSTEDSMHBF`, the full path to your BUSTED[S]-MH HyPhy Batch file, typically this is located within `./hyphy-analyses/BUSTED-MH/BUSTED-MH.bf`


