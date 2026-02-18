library(tidyverse)
library(pheatmap)

# Load file
df <- read.delim("AAI_OUT.csv", sep = ",", header = TRUE)

# Rename problematic column
df <- df %>%
  rename(Two_way_AAI = `Two.way_AAI`)

# Trim whitespace
df$Genome1 <- trimws(df$Genome1)
df$Genome2 <- trimws(df$Genome2)

# Keep required columns
df2 <- df %>%
  select(Genome1, Genome2, Two_way_AAI)

# Add reverse combinations
df_rev <- df2 %>%
  rename(Genome1 = Genome2,
         Genome2 = Genome1)

df_full <- bind_rows(df2, df_rev)

# Get all genomes
genomes <- unique(c(df_full$Genome1, df_full$Genome2))

# Add diagonal
diag_df <- data.frame(
  Genome1 = genomes,
  Genome2 = genomes,
  Two_way_AAI = 100
)

df_full <- bind_rows(df_full, diag_df)

# Remove duplicates
df_full <- df_full %>%
  distinct(Genome1, Genome2, .keep_all = TRUE)

# Create matrix
aai_matrix <- df_full %>%
  pivot_wider(names_from = Genome2,
              values_from = Two_way_AAI)

aai_matrix <- as.data.frame(aai_matrix)
rownames(aai_matrix) <- aai_matrix$Genome1
aai_matrix$Genome1 <- NULL
aai_matrix <- as.matrix(aai_matrix)

# Replace NA
aai_matrix[is.na(aai_matrix)] <- 100

# Sort
aai_matrix <- aai_matrix[order(rownames(aai_matrix)),
                         order(colnames(aai_matrix))]

# Save matrix
write.table(aai_matrix,
            file = "AAI_matrix.tsv",
            sep = "\t",
            quote = FALSE)

cat("Matrix created successfully\n")

# Heatmap
pheatmap(
  aai_matrix,
  display_numbers = TRUE,
  number_format = "%.1f",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  filename = "AAI_heatmap.png",
  width = 13,
  height = 10
)

cat("Heatmap saved successfully\n")
