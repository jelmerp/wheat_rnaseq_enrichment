# Load packages
library(DESeq2)
library(tidyverse)
library(readxl)

# Define the input files
counts_file <- "results/tracy/XopAL1_WheatRNAseqCounts.txt"
annot_infile <- "data/ref/Chinese Spring Genome.xlsx"

# Define the output files
meta_file <- "metadata/meta.tsv"
count_mat_file <- "results/DE/count_mat.tsv"
norm_mat_file <- "results/DE/norm_mat.tsv"
dds_file <- "results/DE/dds.rds"
annot_outfile <- "data/ref/arabi_annot.tsv"
  
# Settings
count_colnames <- c("gene",
                    "Mock_R1", "Mock_R2", "Mock_R3",
                    "Xtt_R1",  "Xtt_R2", "Xtt_R3",
                    "XttDxop_R1", "XttDxop_R2", "XttDxop_R3",
                    "Xtu_R1", "Xtu_R2", "Xtu_R3", 
                    "XtuXopAL1_R1", "XtuXopAL1_R2", "XtuXopAL1_R3")

# Read the input files
count_mat <- read_tsv(counts_file, skip = 1, col_names = count_colnames,
                      show_col_types = FALSE) |>
  column_to_rownames("gene") |>
  as.matrix()
annot_raw <- read_xlsx(annot_infile)

# Create a simple annotation df
annot <- annot_raw |>
  dplyr::select(gene = locusName, description = `Best-hit-arabi-defline`) |>
  distinct(gene, .keep_all = TRUE) |>
  column_to_rownames("gene")

# Create a metadata df
meta <- data.frame(sample = count_colnames[2:length(count_colnames)]) |>
  mutate(treatment = sub("_R\\d", "", sample),
         species = sub("(Xt[tu]).*", "\\1", treatment),
         mutation = case_when(
           treatment == "XttDxop" ~ "DeltaXop",
           treatment == "XtuXopAL1" ~ "PlusXop",
           treatment %in% c("Xtt", "Xtu") ~ "wt",
           .default = "Mock"
         ))
rownames(meta) <- meta$sample

# Create a DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = meta, design = ~ 1)

# Create a normalized count matrix
norm_mat <- assay(rlog(dds))

# Create the output dirs
dir.create("metadata", showWarnings = FALSE, recursive = TRUE)
dir.create("results/DE", showWarnings = FALSE, recursive = TRUE)

# Save the output files
write.table(meta, meta_file, quote = FALSE, sep = "\t")
write.table(count_mat, count_mat_file, quote = FALSE, sep = "\t")
write.table(norm_mat, norm_mat_file, quote = FALSE, sep = "\t")
write.table(annot, annot_outfile, quote = FALSE, sep = "\t")
saveRDS(dds, dds_file)
