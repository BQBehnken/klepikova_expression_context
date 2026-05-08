#!/usr/bin/env Rscript
# aggregate_klepikova.R
# =====================
# Import featureCounts raw counts, normalize with DESeq2,
# average biological replicates per tissue, log2-transform,
# and export the expression matrix for BAT formatting.
#
# Dependencies: DESeq2, dplyr, readr
#   install.packages("BiocManager")
#   BiocManager::install("DESeq2")
#   install.packages(c("dplyr", "readr"))
#
# Usage:
#   Rscript scripts/aggregate_klepikova.R

# I like using pacman, but you can do this via install.packages noted above 
pacman::p_load(BiocManager, DESeq2, dplyr, readr)

# -------------
# Configuration
# -------------

DATADIR     <- "data/klepikova"
COUNTS_FILE <- file.path(DATADIR, "counts", "featurecounts_raw.txt")
RUNINFO     <- file.path(DATADIR, "metadata", "klepikova_runinfo.csv")
OUTPUT      <- file.path(DATADIR, "klepikova_expression_matrix.tsv")

# Six tissues selected for the default BAT dataset
TISSUES_TO_INCLUDE <- c("seed", "root", "stem", "petal", "carpel", "leaf")

# Tidy up column names for the final output
TISSUE_LABELS <- c(
  seed   = "Seed",
  root   = "Root",
  stem   = "Stem",
  petal  = "Petal",
  carpel = "Carpel",
  leaf   = "Leaf"
)

# ---------------------------------
# Step 1: Read featureCounts output
# ---------------------------------

cat("=== Reading featureCounts output ===\n")

fc_raw <- read.delim(COUNTS_FILE, comment.char = "#", stringsAsFactors = FALSE)

# featureCounts columns: Geneid, Chr, Start, End, Strand, Length, <sample BAM paths>
# Extract gene IDs and count columns (columns 7 onward)
gene_ids   <- fc_raw$Geneid
count_data <- fc_raw[, 7:ncol(fc_raw), drop = FALSE]

# Clean column names: strip path prefixes, keep just SRR accession
colnames(count_data) <- gsub(".*\\/", "", colnames(count_data))        # strip path
colnames(count_data) <- gsub("\\.accepted_hits\\.bam$", "", colnames(count_data))  # strip suffix
colnames(count_data) <- gsub("_trimmed$", "", colnames(count_data))    # strip _trimmed if present

rownames(count_data) <- gene_ids

cat(sprintf("  Genes: %d\n", nrow(count_data)))
cat(sprintf("  Samples: %d\n", ncol(count_data)))

# --------------------------------------------------
# Step 2: Build sample → tissue mapping from runinfo
# --------------------------------------------------

cat("=== Reading sample metadata ===\n")

runinfo <- read_csv(RUNINFO, show_col_types = FALSE)

# The runinfo CSV from NCBI has columns including: Run, SampleName, LibraryName
# Tissue labels may be in SampleName or LibraryName — inspect and adapt
# We try SampleName first; if that doesn't contain tissue info, try LibraryName
# This chunk joins the count matrix (SRR accessions as column names) against the metadata CSV to create a LUT of each SRR to its tissue label
sample_meta <- runinfo %>%
  select(Run, SampleName) %>%
  filter(Run %in% colnames(count_data)) %>%
  distinct(Run, .keep_all = TRUE)

# Assign tissue label by matching TISSUES_TO_INCLUDE against SampleName (case-insensitive)
sample_meta <- sample_meta %>%
  mutate(tissue = NA_character_)

for (tissue in TISSUES_TO_INCLUDE) {
  sample_meta <- sample_meta %>%
    mutate(tissue = ifelse(
      is.na(tissue) & grepl(tissue, SampleName, ignore.case = TRUE),
      tissue,
      tissue
    ))
}

# Filter to samples that matched a target tissue
sample_meta <- sample_meta %>%
  filter(!is.na(tissue))

if (nrow(sample_meta) == 0) {
  stop("ERROR: No samples matched the target tissues. ",
       "Check tissue labels in the runinfo CSV vs. TISSUES_TO_INCLUDE.\n",
       "You may need to adjust the matching logic in this script.")
}

cat(sprintf("  Matched %d samples across %d tissues:\n",
            nrow(sample_meta), n_distinct(sample_meta$tissue)))
print(table(sample_meta$tissue))

# Subset count matrix to matched samples
count_data <- count_data[, sample_meta$Run, drop = FALSE]

# ----------------------------
# Step 3: DESeq2 normalization
# ----------------------------

# DESeq2's size factor method works by computing the ratio for each gene its count to the geometric mean of that gene across all samples. 
# The median of all those ratios for a given sample is its size factor. Dividing each sample's raw counts by its size factor puts all samples on comparable footing. 
# This is opposed to the RPKM normalization for the eFP browser.
# The primary issue with RPKM is compositional bias. RPKM normalizes by the total number of mapped reads in a sample. 
# However, if a specific tissue (e.g., a developing seed) highly expresses a few specific genes, those transcripts will consume a massive proportion of the sequencing reads. 
# In an RPKM calculation, this artificially depresses the read counts and consequently the RPKM values for all other genes in that sample, 
# making them look downregulated compared to a leaf even if their absolute transcript numbers are identical.
# We discuss this in the materials and methods for the paper and provide a supplementary table that shows that the output ratios are comparable to the values 
# given in the eFP browser with the important caveat that this holds for tissues in which the majority of genes are not differentially expressed 
# (AKA the baseline genes being expressed are consistent tissue-to-tissue). In extreme cases like seed tissue, interpreting relative fold changes requires additional caution.

cat("\n=== DESeq2 size-factor normalization ===\n")

col_data <- data.frame(
  row.names = sample_meta$Run,
  tissue    = factor(sample_meta$tissue)
)

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(count_data),
  colData   = col_data,
  design    = ~ tissue
)

dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

cat(sprintf("  Size factors: %s\n",
            paste(round(sizeFactors(dds), 3), collapse = ", ")))

# Lastly, if a huge proportion of your genes are zero in some samples (e.g., you source a low-quality library), 
# the size factors can be estimated from a severely reduced gene set and may be unreliable.

# -------------------------------------
# Step 4: Average replicates per tissue
# -------------------------------------

cat("=== Averaging replicates per tissue ===\n")

tissue_means <- sapply(TISSUES_TO_INCLUDE, function(t) {
  samples <- sample_meta$Run[sample_meta$tissue == t]
  if (length(samples) == 1) {
    return(norm_counts[, samples])
  }
  rowMeans(norm_counts[, samples, drop = FALSE])
})

colnames(tissue_means) <- TISSUE_LABELS[colnames(tissue_means)]

# ----------------------
# Step 5: Log2 transform
# ----------------------

cat("=== Log2(normalized counts + 1) transform ===\n")

# The +1 pseudocount is necessary for genes with a normalized count of 0; otherwise it would give log2(0) = -Inf, which would propagate as NA or cause errors downstream

log2_expr <- log2(tissue_means + 1)

# Round to 2 decimal places for readability
log2_expr <- round(log2_expr, 2)

# --------------------------------------------------
# Step 6: Clean gene IDs to gene-level TAIR10 format
# --------------------------------------------------

cat("=== Cleaning gene IDs ===\n")

# Strip transcript isoform suffixes (e.g., AT1G01010.1 → AT1G01010)
clean_ids <- gsub("\\.[0-9]+$", "", rownames(log2_expr))

# If duplicates after stripping, average them
if (any(duplicated(clean_ids))) {
  cat("  Collapsing isoform duplicates...\n")
  log2_expr <- as.data.frame(log2_expr)
  log2_expr$gene_id <- clean_ids
  log2_expr <- log2_expr %>%
    group_by(gene_id) %>%
    summarise(across(everything(), mean), .groups = "drop")
  clean_ids <- log2_expr$gene_id
  log2_expr <- as.matrix(log2_expr[, -1])
  rownames(log2_expr) <- clean_ids
}

cat(sprintf("  Final gene count: %d\n", nrow(log2_expr)))

# --------------
# Step 7: Export
# --------------

cat(sprintf("=== Exporting to %s ===\n", OUTPUT))

out_df <- data.frame(taxa = rownames(log2_expr), log2_expr,
                     check.names = FALSE, stringsAsFactors = FALSE)

write_tsv(out_df, OUTPUT)

cat(sprintf("  Dimensions: %d genes x %d tissues\n",
            nrow(out_df), ncol(out_df) - 1))
cat("  Done.\n")
cat("\nNext step:\n")
cat("  python3 scripts/format_klepikova_bat.py\n")
