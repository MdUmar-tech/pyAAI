# pyAAI
Python implementation of the AAI algorithm following the methodology of Rodriguez-R & Konstantinidis (2016), using reciprocal best hits for genome similarity estimation.



#Install Dependencies
conda install -c bioconda diamond prodigal
pip install biopython pandas numpy
## ⚙️ Requirements

- Python ≥ 3.8
- DIAMOND or BLAST+
- Prodigal (for protein prediction)

Python packages:
- biopython
- pandas
- numpy

---

## 🔬 Overview

AAI (Average Amino Acid Identity) is a genome-wide similarity metric widely used for:

- Prokaryotic taxonomy
- Species delineation
- Phylogenomic comparisons
- Genome similarity estimation

pyAAI:

- Uses reciprocal best hits (RBH)
- Supports DIAMOND or BLASTP
- Performs all-vs-all comparisons
- Generates detailed AAI statistics


Step 1: Predict Proteins
If starting from genome FASTA files:

python multiprodigal.py -i genomes_folder -o proteins_folder
This will generate .faa protein files for each genome.

Step 2: Run pyAAI
python pyAAI.py \
    -i proteins_folder \
    -o all_vs_all_aai.csv \
    -p diamond \
    -id 20 \
    -l 0 \
    -s 50 \
    -L 0.0 \
    -n 50 \
    -t 8



About

Original Implementation

This Python implementation is inspired by the original aai.rb script from the enveomics collection:

Original Ruby script:
https://github.com/lmrodriguezr/enveomics/blob/master/Scripts/aai.rb

Rodriguez-R, L.M.; Konstantinidis, K.T. (2016).
The enveomics collection: A toolbox for specialized analyses of microbial genomes and metagenomes. PeerJ Preprints 4:e1900v1.
