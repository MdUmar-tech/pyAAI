#!/usr/bin/env python3

import argparse
import subprocess
from Bio import SeqIO
import pandas as pd
import numpy as np
import os
import sys

# =====================================================
# ARGPARSE
# =====================================================
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Calculates the Average Amino Acid Identity (AAI) between two genomes"
    )

    parser.add_argument("-1", "--seq1", required=True, help="FASTA file of genome 1 (proteins)")
    parser.add_argument("-2", "--seq2", required=True, help="FASTA file of genome 2 (proteins)")
    parser.add_argument("-l", "--len", type=int, default=0, help="Minimum alignment length (default: 0)")
    parser.add_argument("-L", "--len-fraction", type=float, default=0.0, help="Minimum alignment length fraction (0-1)")
    parser.add_argument("-i", "--id", type=float, default=20, help="Minimum identity percentage (default: 20)")
    parser.add_argument("-s", "--bitscore", type=float, default=0, help="Minimum bitscore (default: 0)")
    parser.add_argument("-n", "--hits", type=int, default=50, help="Minimum number of RBH hits (default: 50)")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads (default: 1)")
    parser.add_argument("-d", "--dec", type=int, default=2, help="Decimal positions (default: 2)")
    parser.add_argument("-R", "--rbm", help="Save reciprocal best matches to file")
    parser.add_argument("-a", "--auto", action="store_true",help="Only print AAI value")
    parser.add_argument("-q", "--quiet", action="store_true", help="Run quietly")

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
# BUILD DIAMOND DATABASES
# =====================================================
if not args.quiet:
    print("Building DIAMOND databases...")

subprocess.run(["diamond", "makedb", "--in", TMP1, "-d", DB1], check=True)
subprocess.run(["diamond", "makedb", "--in", TMP2, "-d", DB2], check=True)

# =====================================================
# RUN DIAMOND
# =====================================================
if not args.quiet:
    print("Running DIAMOND...")

subprocess.run([
    "diamond", "blastp",
    "--query", TMP1,
    "--db", DB2,
    "--threads", str(THREADS),
    "--sensitive",
    "--outfmt", "6",
    "--out", BLAST1
], check=True)

subprocess.run([
    "diamond", "blastp",
    "--query", TMP2,
    "--db", DB1,
    "--threads", str(THREADS),
    "--sensitive",
    "--outfmt", "6",
    "--out", BLAST2
], check=True)

# =====================================================
# LOAD RESULTS
# =====================================================
cols = [
    "qseqid","sseqid","pident","length","mismatch","gapopen",
    "qstart","qend","sstart","send","evalue","bitscore"
]

df1 = pd.read_csv(BLAST1, sep="\t", names=cols, dtype=str)
df2 = pd.read_csv(BLAST2, sep="\t", names=cols, dtype=str)

# Convert numeric columns explicitly
numeric_cols = cols[2:]
df1[numeric_cols] = df1[numeric_cols].apply(pd.to_numeric)
df2[numeric_cols] = df2[numeric_cols].apply(pd.to_numeric)

# Ensure IDs stay as clean strings (no .0)
df1["qseqid"] = df1["qseqid"].astype(str)
df1["sseqid"] = df1["sseqid"].astype(str)
df2["qseqid"] = df2["qseqid"].astype(str)
df2["sseqid"] = df2["sseqid"].astype(str)

# =====================================================
# FIRST VALID HIT (Ruby-style behavior)
# =====================================================
def get_first_valid_hits(df, lenA, lenB):
    seen = set()
    valid_rows = []

    for _, row in df.iterrows():
        q = str(row["qseqid"])
        s = str(row["sseqid"])

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
# ONE-WAY AAI 1
# =====================================================
n1 = len(df1_best)
if n1 > 0:
    id_sum1 = df1_best["pident"].sum()
    sq_sum1 = (df1_best["pident"] ** 2).sum()
    mean1 = id_sum1 / n1
    sd1 = np.sqrt(sq_sum1/n1 - mean1**2)
else:
    mean1 = 0.0
    sd1 = 0.0

# =====================================================
# ONE-WAY AAI 2
# =====================================================
n2 = len(df2_best)
if n2 > 0:
    id_sum2 = df2_best["pident"].sum()
    sq_sum2 = (df2_best["pident"] ** 2).sum()
    mean2 = id_sum2 / n2
    sd2 = np.sqrt(sq_sum2/n2 - mean2**2)
else:
    mean2 = 0.0
    sd2 = 0.0

# =====================================================
# RECIPROCAL BEST HITS (Two-way AAI)
# =====================================================
forward = {str(row["qseqid"]): str(row["sseqid"])
           for _, row in df1_best.iterrows()}

id_sum = 0.0
sq_sum = 0.0
n_hits = 0
rbh_pairs = []

for _, row in df2_best.iterrows():
    q2 = str(row["qseqid"]).split(".")[0]
    s2 = str(row["sseqid"]).split(".")[0]

    if s2 in forward and forward[s2] == q2:
        pid = row["pident"]
        id_sum += pid
        sq_sum += pid**2
        n_hits += 1

        rbh_pairs.append([
            id_map1[s2],
            id_map2[q2],
            pid
        ])

# =====================================================
# FINAL OUTPUT
# =====================================================
if n_hits < MIN_HITS:
    if not args.auto:
        print(f"Insufficient hits ({n_hits} < {MIN_HITS})")
    sys.exit(1)

mean_aai = id_sum / n_hits
sd_aai = np.sqrt(sq_sum/n_hits - mean_aai**2)

if args.auto:
    print(f"{mean_aai:.{DECIMALS}f}")
else:
    print("\n==============================")
    print(f"! One-way AAI 1: {mean1:.{DECIMALS}f}% "
          f"(SD: {sd1:.{DECIMALS}f}%), from {n1} proteins.")
    print(f"! One-way AAI 2: {mean2:.{DECIMALS}f}% "
          f"(SD: {sd2:.{DECIMALS}f}%), from {n2} proteins.")
    print(f"! Two-way AAI  : {mean_aai:.{DECIMALS}f}% "
          f"(SD: {sd_aai:.{DECIMALS}f}%), from {n_hits} proteins.")
    print("==============================")

# =====================================================
# SAVE RBM
# =====================================================
if args.rbm:
    rbh_df = pd.DataFrame(
        rbh_pairs,
        columns=["genome1_protein","genome2_protein","pident"]
    )
    rbh_df.to_csv(args.rbm, index=False)
#################################################################
#################################################################
# Example usage:
# python aai_diamond.py \
#   -1 proteins_1.faa \
#   -2 proteins_2.proteins.faa \
#   -i 20 -l 0 -s 50 -L 0.0 -n 50 -t 8
#################################################################
