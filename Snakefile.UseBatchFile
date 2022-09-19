# Snakefile for BUSTEDS-MH Analysis
# @Author: Alexander G Lucaci

# Imports
import os
import sys
import json
import csv
from pathlib import Path
import glob

# Declares ------------------------------------------------------------
with open("cluster.json", "r") as in_c:
  cluster = json.load(in_c)

# Set the data directory 
#FASTA_DIR = "/home/aglucaci/BUSTEDS-MH/data/Enard"
#TREE_DIR = "/home/aglucaci/BUSTEDS-MH/data/Enard/Trees/BioNJ/ForHyPhy"


FASTA_DIR = "/home/aglucaci/BUSTEDS-MH/data/Shultz/compgen_alignments"
TREE_DIR = "/home/aglucaci/BUSTEDS-MH/data/Shultz/compgen_alignments/Trees/BioNJ/ForHyPhy"

# glob all of the files
FASTAS = glob.glob(FASTA_DIR + '/*.phy')
TREES = glob.glob(TREE_DIR + '/*.nwk')

FASTA_filenames = [os.path.basename(x) for x in glob.glob(FASTA_DIR + '/*.phy')]
TREE_filenames = [os.path.basename(x) for x in glob.glob(TREE_DIR + '/*.nwk')]

# Report to user
print("# Number of fasta files to process:", len(FASTAS))
print("# Number of accompanying tree files to process:", len(TREES))

#OUTDIR = "/home/aglucaci/BUSTEDS-MH/analysis/Enard"
OUTDIR = "/home/aglucaci/BUSTEDS-MH/analysis/Shultz"
OUTDIR_BUSTEDSMH = os.path.join(OUTDIR, "BUSTEDS-MH")
OUTDIR_BUSTEDS    = os.path.join(OUTDIR, "BUSTEDS")

# Report to user
print("# Files for BUSTEDS-MH and BUSTEDS will be saved in:", OUTDIR)

# Create output dir.
Path(OUTDIR).mkdir(parents=True, exist_ok=True)
Path(OUTDIR_BUSTEDSMH).mkdir(parents=True, exist_ok=True)
Path(OUTDIR_BUSTEDS).mkdir(parents=True, exist_ok=True)

# Settings, these can be passed in or set in a config.json type file
PPN = cluster["__default__"]["ppn"] 

BUSTEDSMH_bf = "/home/aglucaci/hyphy-analyses/BUSTED-MH/BUSTED-MH.bf"

hyphy = "HYPHYMPI"

# Examples:
# ENSG00000120158.fa-BioNJ_tree_hyphy.nwk
# ENSG00000120451.fa


# Rule all
rule all:
    input:
        expand(os.path.join(OUTDIR_BUSTEDSMH, "{fasta}.BUSTEDS-MH.json"), fasta=FASTA_filenames),
        expand(os.path.join(OUTDIR_BUSTEDS, "{fasta}.BUSTEDS.json"), fasta=FASTA_filenames)

# Individual rules

rule BUSTEDSMH:
    input:
        fasta = os.path.join(FASTA_DIR, "{sample}"),
        tree  = os.path.join(TREE_DIR, "{sample}-BioNJ_tree_hyphy.nwk") 
    output:
        output = os.path.join(OUTDIR_BUSTEDSMH, "{sample}.BUSTEDS-MH.json")
    conda: 'environment.yml'        
    shell:
        "mpirun --use-hwthread-cpus -np {PPN} {hyphy} {BUSTEDSMH_bf} --alignment {input.fasta} --tree {input.tree} --output {output.output} --starting-points 10"
    #end shell
#end fule BUSTEDSMH

rule BUSTEDS:
    input:
        fasta = os.path.join(FASTA_DIR, "{sample}"),
        tree  = os.path.join(TREE_DIR, "{sample}-BioNJ_tree_hyphy.nwk") 
    output:
        output = os.path.join(OUTDIR_BUSTEDS, "{sample}.BUSTEDS.json")
    conda: 'environment.yml'
    shell:
        "mpirun --use-hwthread-cpus -np {PPN} {hyphy} BUSTED --alignment {input.fasta} --tree {input.tree} --output {output.output} --starting-points 10"
    #end shell
#end rule BUSTEDS


# End of file
