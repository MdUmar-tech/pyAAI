#!/usr/bin/env python3

import argparse
import subprocess
from Bio import SeqIO
import pandas as pd
import numpy as np
import sys

# =====================================================
# ARGPARSE
# =====================================================
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Calculate Average Amino Acid Identity (AAI) between two genomes using BLASTP"
    )

    parser.add_argument("-1", "--seq1", required=True,
                        help="Genome 1 protein FASTA")
    parser.add_argument("-2", "--seq2", required=True,
                        help="Genome 2 protein FASTA")

    parser.add_argument("-l", "--len", type=int, default=0,
                        help="Minimum alignment length")
    parser.add_argument("-L", "--len_fraction", type=float, default=0.0,
                        help="Minimum alignment length fraction")
    parser.add_argument("-i", "--id", type=float, default=20,
                        help="Minimum percent identity")
    parser.add_argument("-s", "--bitscore", type=float, default=0,
                        help="Minimum bitscore")
    parser.add_argument("-n", "--hits", type=int, default=50,
                        help="Minimum RBH hits required")

    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Threads")
    parser.add_argument("-d", "--dec", type=int, default=2,
                        help="Decimal places")
    parser.add_argument("-R", "--rbm",
                        help="Output RBH table")
    parser.add_argument("-a", "--auto", action="store_true",
                        help="Only print final AAI")
    parser.add_argument("-q", "--quiet", action="store_true",
                        help="Quiet mode")

    return parser.parse_args()


args = parse_arguments()

FASTA1 = args.seq1
FASTA2 = args.seq2

MIN_LEN = args.len
MIN_LEN_FRAC = args.len_fraction
MIN_ID = args.id
MIN_BITSCORE = args.bitscore
MIN_HITS = args.hits
THREADS = args.threads
DECIMALS = args.dec

TMP1 = "tmp1.faa"
TMP2 = "tmp2.faa"
DB1 = "db1"
DB2 = "db2"
BLAST1 = "1_vs_2.tsv"
BLAST2 = "2_vs_1.tsv"

# =====================================================
# FASTA RENUMBERING (Ruby behavior)
# =====================================================
def renumber_fasta(input_fasta, output_fasta):
    id_map = {}
    lengths = {}
    counter = 1

    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            out.write(f">{counter}\n{str(record.seq)}\n")
            id_map[str(counter)] = record.id
            lengths[str(counter)] = len(record.seq)
            counter += 1

    return id_map, lengths


if not args.quiet:
    print("Renumbering FASTA files...")

id_map1, len1 = renumber_fasta(FASTA1, TMP1)
id_map2, len2 = renumber_fasta(FASTA2, TMP2)

# =====================================================
# BUILD BLAST DATABASES
# =====================================================
if not args.quiet:
    print("Building BLAST databases...")

subprocess.run([
    "makeblastdb",
    "-in", TMP1,
    "-dbtype", "prot",
    "-out", DB1
], check=True)

subprocess.run([
    "makeblastdb",
    "-in", TMP2,
    "-dbtype", "prot",
    "-out", DB2
], check=True)

# =====================================================
# RUN BLASTP
# =====================================================
if not args.quiet:
    print("Running BLASTP...")

subprocess.run([
    "blastp",
    "-query", TMP1,
    "-db", DB2,
    "-num_threads", str(THREADS),
    "-outfmt", "6",
    "-max_target_seqs", "1",
    "-out", BLAST1
], check=True)

subprocess.run([
    "blastp",
    "-query", TMP2,
    "-db", DB1,
    "-num_threads", str(THREADS),
    "-outfmt", "6",
    "-max_target_seqs", "1",
    "-out", BLAST2
], check=True)

# =====================================================
# LOAD RESULTS
# =====================================================
cols = [
    "qseqid","sseqid","pident","length","mismatch","gapopen",
    "qstart","qend","sstart","send","evalue","bitscore"
]

df1 = pd.read_csv(BLAST1, sep="\t", names=cols)
df2 = pd.read_csv(BLAST2, sep="\t", names=cols)

# Ensure numeric columns
numeric_cols = cols[2:]
df1[numeric_cols] = df1[numeric_cols].apply(pd.to_numeric)
df2[numeric_cols] = df2[numeric_cols].apply(pd.to_numeric)

df1["qseqid"] = df1["qseqid"].astype(str)
df1["sseqid"] = df1["sseqid"].astype(str)
df2["qseqid"] = df2["qseqid"].astype(str)
df2["sseqid"] = df2["sseqid"].astype(str)

# =====================================================
# FIRST VALID HIT (Ruby behavior)
# =====================================================
def get_first_valid_hits(df, lenA, lenB):
    seen = set()
    valid_rows = []

    for _, row in df.iterrows():
        q = row["qseqid"]
        s = row["sseqid"]

        if q in seen:
            continue

        min_len = min(lenA.get(q, 0), lenB.get(s, 0))

        if (row["length"] >= MIN_LEN and
            row["length"] >= MIN_LEN_FRAC * min_len and
            row["pident"] >= MIN_ID and
            row["bitscore"] >= MIN_BITSCORE):

            seen.add(q)
            valid_rows.append(row)

    return pd.DataFrame(valid_rows)


df1_best = get_first_valid_hits(df1, len1, len2)
df2_best = get_first_valid_hits(df2, len2, len1)

# =====================================================
# ONE-WAY AAI
# =====================================================
def compute_aai(df):
    n = len(df)
    if n == 0:
        return 0.0, 0.0, 0
    id_sum = df["pident"].sum()
    sq_sum = (df["pident"] ** 2).sum()
    mean = id_sum / n
    sd = np.sqrt(sq_sum/n - mean**2)
    return mean, sd, n


mean1, sd1, n1 = compute_aai(df1_best)
mean2, sd2, n2 = compute_aai(df2_best)

# =====================================================
# RECIPROCAL BEST HITS
# =====================================================
forward = {row["qseqid"]: row["sseqid"]
           for _, row in df1_best.iterrows()}

id_sum = 0
sq_sum = 0
n_hits = 0
rbh_pairs = []

for _, row in df2_best.iterrows():
    q2 = row["qseqid"]
    s2 = row["sseqid"]

    if s2 in forward and forward[s2] == q2:
        pid = row["pident"]
        id_sum += pid
        sq_sum += pid**2
        n_hits += 1
        rbh_pairs.append([id_map1[s2], id_map2[q2], pid])

if n_hits < MIN_HITS:
    if not args.auto:
        print(f"Insufficient hits ({n_hits} < {MIN_HITS})")
    sys.exit(1)

mean_aai = id_sum / n_hits
sd_aai = np.sqrt(sq_sum/n_hits - mean_aai**2)

# =====================================================
# OUTPUT
# =====================================================
if args.auto:
    print(f"{mean_aai:.{DECIMALS}f}")
else:
    print("\n==============================")
    print(f"! One-way AAI 1: {mean1:.{DECIMALS}f}% (SD: {sd1:.{DECIMALS}f}%), from {n1} proteins.")
    print(f"! One-way AAI 2: {mean2:.{DECIMALS}f}% (SD: {sd2:.{DECIMALS}f}%), from {n2} proteins.")
    print(f"! Two-way AAI  : {mean_aai:.{DECIMALS}f}% (SD: {sd_aai:.{DECIMALS}f}%), from {n_hits} proteins.")
    print("==============================")

# =====================================================
# SAVE RBH
# =====================================================
if args.rbm:
    rbh_df = pd.DataFrame(
        rbh_pairs,
        columns=["genome1_protein","genome2_protein","pident"]
    )
    rbh_df.to_csv(args.rbm, index=False)

#python aai_blast.py \
#-1 proteins/Streptomyces_sp._MMCC_100.proteins.faa \
#-2 proteins/Streptomyces_rubrogriseus_NBRC_15455.proteins.faa \
#-i 20 \
#-l 0 \
#-s 50 \
#-L 0.0 \
#-n 50 \
#-t 8
