import os
import subprocess

# Define input directory
input_dir = "genomes"

# Define separate output folders
gene_dir = "genes"
protein_dir = "proteins"

# Create folders if they don't exist
os.makedirs(gene_dir, exist_ok=True)
os.makedirs(protein_dir, exist_ok=True)

# List all files in input directory
input_files = os.listdir(input_dir)

# Process each .fna file
for file in input_files:
    if file.endswith(".fna"):
        
        base_name = os.path.splitext(file)[0]

        gene_output = os.path.join(gene_dir, base_name + ".genes")
        protein_output = os.path.join(protein_dir, base_name + ".proteins.faa")

        command = [
            "prodigal",
            "-i", os.path.join(input_dir, file),
            "-o", gene_output,
            "-a", protein_output
        ]

        subprocess.run(command)

print("All genomes processed successfully.")
