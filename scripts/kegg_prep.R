# This script will create a wheat gene to KEGG pathway lookup file
# The input is the gene-to-KO (Kegg Ortholog) info in the wheat annotation file
# We can download the wheat KEGG pathway info using clusterProfiler,
# and download the KO-to-pathway links for the wheat pathways using KEGGREST.

# KEGG does have an entry for wheat:
# https://www.kegg.jp/kegg-bin/show_organism?menu_type=pathway_maps&org=taes
# But I can't find a way to directly work with the wheat KEGG annotations via clusterProfiler,
# since we have Ensembl gene IDs instead of NCBI gene IDs.
# No conversions seem to exist - tried biomaRt (no wheat data) and 'gget info TraesCS1A03G0007900'.
# This is likely because it is a completely separate annotation.

# SETUP ------------------------------------------------------------------------
pkgs <- c("readxl", "tidyverse", "clusterProfiler", "KEGGREST")
pacman::p_load(char = pkgs)

# Source script with functions
source("mcic-scripts/rnaseq/rfuns/kegg-db_funs.R")

# Define input files
annot_file <- "data/ref/Chinese Spring Genome.xlsx" # Wheat annotation file with K-annotations

# Define output files
outdir <- "results/kegg"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
kegg_outfile <- file.path(outdir, "kegg_map.tsv")


# READ INPUT FILES -------------------------------------------------------------
# Wheat annotation
annot <- read_xlsx(annot_file) |>
  select(gene = locusName, KO) |>
  drop_na() |>
  distinct() # duplicates due to transcripts
nrow(annot) # 24,775

# Get wheat-specific pathways with clusterProfiler
pathway_info_raw <- clusterProfiler::download_KEGG(species = "taes")
pathway_info <- pathway_info_raw$KEGGPATHID2NAME |>
  rename(pathway = from, description = to)

# For each wheat-specific pathway, get the KOs linked to it
pw2ko_df_raw <- map_dfr(unique(pathway_info$pathway), pw2kos)
pw2ko_df <- pw2ko_df_raw |> select(pathway, KO) |> drop_na()
# Pathways without KOs:
#>taes01100 ... Metabolic pathways - Triticum aestivum (bread wheat) - No KOs found!
#>taes01110 ... Biosynthesis of secondary metabolites - Triticum aestivum (bread wheat) - No KOs found!
#>taes01200 ... Carbon metabolism - Triticum aestivum (bread wheat) - No KOs found!
#>taes01210 ... 2-Oxocarboxylic acid metabolism - Triticum aestivum (bread wheat) - No KOs found!
#>taes01212 ... Fatty acid metabolism - Triticum aestivum (bread wheat) - No KOs found!
#>taes01230 ... Biosynthesis of amino acids - Triticum aestivum (bread wheat) - No KOs found!
#>taes01232 ... Nucleotide metabolism - Triticum aestivum (bread wheat) - No KOs found!
#>taes01250 ... Biosynthesis of nucleotide sugars - Triticum aestivum (bread wheat) - No KOs found!
#>taes01240 ... Biosynthesis of cofactors - Triticum aestivum (bread wheat) - No KOs found!

# Replace the pathway IDs with the tomato-specific pathways
kegg_map <- inner_join(annot, pw2ko_df, by = "KO", relationship = "many-to-many") |>
  left_join(pathway_info, by = "pathway") |>
  arrange(pathway) |>
  select(pathway, gene, description, KO)

# Report a few stats
message("Nr of genes with at least one pathway: ", length(unique(kegg_map$gene))) #> 14,277
message("Nr of pathways: ", length(unique(kegg_map$pathway))) #> 134

# Write the output file
write_tsv(kegg_map, kegg_outfile)
