# SETUP ------------------------------------------------------------------------
# Load packages
library(readxl)
library(tidyverse)

# Define the input files
infile <- "metadata/Xtt&Xtu_RNAseq(wheat_mapped).xlsx"
col_types <- c("text", "numeric", "numeric", rep("text", 4))

# Define the output files
DE_file <- "metadata/DE.tsv"
annot_file <- "metadata/annot_DE.tsv"

# Read all DE results into a single dataframe
sheet_names <- excel_sheets(infile)[2:length(excel_sheets(infile))]
names(sheet_names) <- sub("∆", "D", sheet_names)
DE_res_raw <- map_dfr(
  .x = sheet_names, .f = read_excel, .id = "contrast",
  path = infile, col_types = col_types, range = "A2:G50000"
  ) |>
  janitor::clean_names()

# Separately
DE_res <- DE_res_raw |>
  select(contrast, gene = gene_id, lfc = log2fold_change, padj)

annot <- DE_res_raw |>
  distinct(gene_id, .keep_all = TRUE) |> 
  select(gene = gene_id,
         arabidopsis_ortholog, arabidopsis_description = annotation,
         rice_ortholog, rice_description = rice_annotation)

# Write the output files
write_tsv(DE_res, DE_file)
write_tsv(annot, annot_file)
