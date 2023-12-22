# SETUP ------------------------------------------------------------------------
# Load packages
pkgs <- c("clusterProfiler", "enrichplot", "org.At.tair.db",
          "ggpubr", "here", "tidyverse", "readxl")
pacman::p_load(char = pkgs)

# Source script with helper functions
source(here("mcic-scripts/rnaseq/rfuns/enrich_funs.R"))
source(here("mcic-scripts/rnaseq/rfuns/DE_funs.R"))

# Settings
set.seed(1)               # Needed for reproducible GSEA results due to permutations
OrgDb <- "org.At.tair.db" # BioConductor organismal database
keyType <- "TAIR"         # Gene ID type in our DE results to match with GO db gene IDs

# Define input files
annot_file <- here("data/ref/Chinese Spring Genome.xlsx")
DE_file <- here("results/tracy/DE.tsv")

# Define output files
GO_ORA_file <- here("results/GO/GO_ORA_res.tsv")
GO_GSEA_file <- here("results/GO/GO_GSEA_res.tsv")


# READ AND PREP INPUT FILES ----------------------------------------------------
# Read the input - annotation
annot <- read_xlsx(annot_file) |>
  dplyr::select(gene = locusName,
                arabi_id = `Best-hit-arabi-name`,
                arabi_name = `Best-hit-arabi-defline`) |>
    dplyr::filter(!is.na(arabi_id)) |>
    distinct(gene, .keep_all = TRUE)

# Read the input - DE results
DE_res <- read_tsv(DE_file, show_col_types = FALSE) |>
  mutate(isDE = ifelse(padj < 0.05, TRUE, FALSE)) |>
  left_join(annot, by = "gene") |>
  dplyr::rename(gene_wheat = gene) |>
  mutate(gene = ifelse(!is.na(arabi_id), arabi_id, gene_wheat))

# Define a vector of contrasts to iterate over
contrasts <- unique(DE_res$contrast)


# RUN THE GO ORA ANALYSIS ------------------------------------------------------
# Prep df with all combinations of contrasts and DE direction
combs <- expand_grid(contrasts, DE_dirs = c("down", "up"))

# Run
enrich_res <- map2_dfr(
  .x = combs$contrasts, .y = combs$DE_dirs, .f = run_enrich,
  OrgDb = OrgDb, keyType = keyType, ontology_type = "GO",
  DE_res = DE_res, allow_dups = TRUE, return_df = TRUE,
  )

# Write output to file
write_tsv(enrich_res, GO_ORA_file)


# RUN THE GO GSEA ANALYSIS -----------------------------------------------------
# NOTE: Because Arabidopsis genes often correspond to multiple wheat genes,
#       in order to use Log-fold change values, we'll have to compute means across such 'duplicated' genes.
#       The function below does this automatically when 'allow_dups=TRUE', but this is something to keep in mind!
gsea_res <- map_dfr(
  .x = contrasts, .f = run_gsea,
  OrgDb = OrgDb, keyType = keyType,
  DE_res = DE_res, allow_dups = TRUE, return_df = TRUE
  )

# Write output to file
write_tsv(gsea_res, GO_GSEA_file)
