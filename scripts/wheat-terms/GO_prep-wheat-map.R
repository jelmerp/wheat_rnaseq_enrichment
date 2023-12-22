# SETUP ------------------------------------------------------------------------
# Load packages
pkgs <- c("readxl", "tidyverse", "GO.db")
pacman::p_load(char = pkgs)

# Define the input files
infile <- "data/ref/Chinese Spring Genome.xlsx"

# Define the output files
outdir <- "results/GO"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(outdir, "GO_map.tsv")

# Read the input files
annot <- read_xlsx(infile) 


# CREATE THE GO MAP ------------------------------------------------------------
# Create a df with term-to-description info
go_info <- AnnotationDbi::select(GO.db,
                                 columns = c("GOID", "TERM", "ONTOLOGY"),
                                 keys = keys(GO.db, keytype = "GOID"),
                                 keytype = "GOID") |>
  dplyr::rename(term = GOID, description = TERM, ontology = ONTOLOGY)

# Create a df with gene-to-term relationships
go_map <- annot |>
  dplyr::select(term = GO, gene = locusName) |>
  filter(!is.na(term)) |>
  separate_longer_delim(term, delim = " ") |>
  distinct() |>
  left_join(go_info, by = "term") |>
  arrange(term) |>
  dplyr::select(term, gene, description, ontology)
nterms_before <- length(unique(go_map$term))

# Remove deprecated GO terms with no description
go_map <- go_map |> filter(!is.na(description))
n_deprecated <- nterms_before - length(unique(go_map$term))

# Report a few stats
# NOTE: Keeping deprecated GO terms here, but may want to remove them
message("Nr of genes with at least one GO term: ", length(unique(go_map$gene)))
message("Nr of gene-to-GOterm entries: ", nrow(go_map))
message("Nr of deprecated GO terms removed: ", n_deprecated)

# Write the final df to file
write_tsv(go_map, outfile)
