# klepikova_expression_context
We provide a default expression context for the BLAST-align-tree (BAT) by constructing a gene-level Arabidopsis thaliana expression atlas sourced from a subset of Klepikova’s developmental transcriptome atlas data (Klepikova et al. 2016).

The order of events are: 
  1. download_klepikova.sh — fetch SRA metadata to inspect for correct downloads, then download FASTQs (See materials and methods for exact tissue ID's we used)
  2. qc_trim_klepikova.sh — FastQC + Trimmomatic on raw FASTQs (Based on eFP Browser; See Sullivan et al. 2019 (10.1111/tpj.14468)).
  3. align_klepikova.sh — build Bowtie2 index, align trimmed reads with TopHat2; creates BAM files (Based on eFP Browser; See Sullivan et al. 2019 (10.1111/tpj.14468))
  4. quantify_klepikova.sh — count reads per gene from BAMs with featureCounts; creates count matrix to prepare for working in the BAT context
  5. aggregate_klepikova.R — DESeq2 normalization, average replicates per tissue, log2-transform for tissue-to-tissue expression matrix
  6. format_klepikova_bat.py — reformat expression matrix into BAT dataset format and places in folder datasets/
