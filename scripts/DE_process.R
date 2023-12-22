# SETUP ------------------------------------------------------------------------
# Load packages
library(readxl)
library(tidyverse)

# Define the input files
infile <- "results/tracy/XopAL1_RNAseq(allDEGs).xlsx"

# Define the output files
outfile <- "results/tracy/DE.tsv"

# Read all DE results into a single dataframe
sheet_names <- excel_sheets(infile)
names(sheet_names) <- sub("\\(all\\)", "", sheet_names)
DE_res <- map_dfr(.x = sheet_names, .f = read_excel, .id = "contrast",
                  path = infile) |>
  select(contrast, gene = `Gene ID`, lfc = log2FoldChange, padj) |>
  mutate(gene = sub("\\..*", "", gene))

# Write the output files
write_tsv(DE_res, outfile)
