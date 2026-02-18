#!/usr/bin/env python3

import argparse
import subprocess
import pandas as pd
import numpy as np
import os
import sys
import glob
from Bio import SeqIO
from itertools import combinations

# =====================================================
# ARGPARSE
# =====================================================
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="All-vs-All AAI calculator (DIAMOND/BLAST)"
    )

    parser.add_argument("-i", "--input_folder", required=True,
                        help="Folder containing .faa files")

    parser.add_argument("-o", "--output_csv", required=True,
                        help="Output CSV file")

    parser.add_argument("-p", "--program", choices=["diamond", "blast"],
                        default="diamond")

    parser.add_argument("-l", "--len", type=int, default=0)
    parser.add_argument("-L", "--len_fraction", type=float, default=0.0)
    parser.add_argument("-id", "--identity", type=float, default=20)
    parser.add_argument("-s", "--bitscore", type=float, default=0)
    parser.add_argument("-n", "--hits", type=int, default=50)
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("-d", "--decimals", type=int, default=3)

    return parser.parse_args()


# =====================================================
# FASTA RENUMBERING
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

    return id_map, lengths, counter - 1


# =====================================================
# BUILD DATABASE
# =====================================================
def build_database(program, fasta, dbname):
    if program == "diamond":
        subprocess.run(["diamond", "makedb", "--in", fasta, "-d", dbname], check=True)
    else:
        subprocess.run(["makeblastdb", "-in", fasta,
                        "-dbtype", "prot", "-out", dbname], check=True)


# =====================================================
# RUN SEARCH
# =====================================================
def run_search(program, query, db, out, threads):
    if program == "diamond":
        subprocess.run([
            "diamond", "blastp",
            "--query", query,
            "--db", db,
            "--threads", str(threads),
            "--sensitive",
            "--outfmt", "6",
            "--out", out
        ], check=True)
    else:
        subprocess.run([
            "blastp",
            "-query", query,
            "-db", db,
            "-num_threads", str(threads),
            "-outfmt", "6",
            "-max_target_seqs", "1",
            "-out", out
        ], check=True)


# =====================================================
# FILTER FIRST VALID HIT
# =====================================================
def get_first_valid_hits(df, lenA, lenB, MIN_LEN, MIN_LEN_FRAC,
                         MIN_ID, MIN_BITSCORE):

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


# =====================================================
# COMPUTE AAI
# =====================================================
def compute_aai(df):
    n = len(df)
    if n == 0:
        return 0.0, 0.0, 0

    mean = df["pident"].mean()
    sd = df["pident"].std(ddof=0)
    return mean, sd, n


# =====================================================
# CORE AAI BETWEEN TWO GENOMES
# =====================================================
def calculate_pairwise_aai(file1, file2, args):

    TMP1 = "tmp1.faa"
    TMP2 = "tmp2.faa"
    DB1 = "db1"
    DB2 = "db2"
    BLAST1 = "1_vs_2.tsv"
    BLAST2 = "2_vs_1.tsv"

    id_map1, len1, total1 = renumber_fasta(file1, TMP1)
    id_map2, len2, total2 = renumber_fasta(file2, TMP2)

    build_database(args.program, TMP1, DB1)
    build_database(args.program, TMP2, DB2)

    run_search(args.program, TMP1, DB2, BLAST1, args.threads)
    run_search(args.program, TMP2, DB1, BLAST2, args.threads)

    cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
            "qstart","qend","sstart","send","evalue","bitscore"]

    df1 = pd.read_csv(BLAST1, sep="\t", names=cols)
    df2 = pd.read_csv(BLAST2, sep="\t", names=cols)

    df1[cols[2:]] = df1[cols[2:]].apply(pd.to_numeric)
    df2[cols[2:]] = df2[cols[2:]].apply(pd.to_numeric)

    df1_best = get_first_valid_hits(
        df1, len1, len2,
        args.len, args.len_fraction,
        args.identity, args.bitscore
    )

    df2_best = get_first_valid_hits(
        df2, len2, len1,
        args.len, args.len_fraction,
        args.identity, args.bitscore
    )

    mean1, sd1, n1 = compute_aai(df1_best)
    mean2, sd2, n2 = compute_aai(df2_best)

    forward = dict(zip(df1_best["qseqid"], df1_best["sseqid"]))

    rbh = df2_best[
        df2_best.apply(
            lambda x: forward.get(x["sseqid"]) == x["qseqid"],
            axis=1
        )
    ]

    mean_rbh, sd_rbh, n_rbh = compute_aai(rbh)

    if n_rbh < args.hits:
        return None

    return {
        "Genome1": os.path.basename(file1).replace(".faa", ""),
        "Genome2": os.path.basename(file2).replace(".faa", ""),
        "Total_seq1": total1,
        "Total_seq2": total2,
        "One-way_AAI_1": mean1,
        "Selected_seq1": n1,
        "SD1": sd1,
        "One-way_AAI_2": mean2,
        "Selected_seq2": n2,
        "SD2": sd2,
        "Two-way_AAI": mean_rbh,
        "Orthologous_seq": n_rbh,
        "SD": sd_rbh
    }


# =====================================================
# MAIN ALL-VS-ALL
# =====================================================
def main():
    args = parse_arguments()

    # -------------------------------------------------
    # Ensure .csv extension
    # -------------------------------------------------
    if not args.output_csv.endswith(".csv"):
        args.output_csv += ".csv"

    # -------------------------------------------------
    # Collect FASTA files
    # -------------------------------------------------
    faa_files = glob.glob(os.path.join(args.input_folder, "*.faa"))

    if len(faa_files) < 2:
        print("Error: Need at least 2 .faa files in input folder.")
        sys.exit(1)

    print(f"Total genomes found: {len(faa_files)}")

    # -------------------------------------------------
    # Header handling (append-safe)
    # -------------------------------------------------
    header_written = os.path.exists(args.output_csv)

    # -------------------------------------------------
    # All-vs-All comparisons
    # -------------------------------------------------
    for file1, file2 in combinations(faa_files, 2):
        print(f"Processing {os.path.basename(file1)} vs {os.path.basename(file2)}")

        res = calculate_pairwise_aai(file1, file2, args)

        if res:
            df_row = pd.DataFrame([res])

            df_row.to_csv(
                args.output_csv,
                mode="a",
                header=not header_written,
                index=False
            )

            header_written = True
        else:
            print("Skipped (insufficient RBH hits)")

    print("\nAll-vs-All AAI completed.")
    print(f"Results saved to: {args.output_csv}")


if __name__ == "__main__":
    main()

