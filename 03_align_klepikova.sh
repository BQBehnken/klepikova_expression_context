#!/usr/bin/env bash
# align_klepikova.sh
# ==================
# Build Bowtie2 index from TAIR10 genome and align trimmed reads with TopHat2.
# Replicates the Klepikova et al. (2016) Galaxy pipeline with exact parameters:
#   -i / --min-intron-length  50
#   -I / --max-intron-length  5000
#
# Dependencies: Bowtie2 (v2.2.8 used in original), TopHat2 (v2.1.1 used in original)
#   conda install -c bioconda tophat bowtie2
# TopHat2 calls Bowtie2 internally and will crash if it can't find the bowtie2 executable in your $PATH. Both need to be installed. TopHat2 was used for reproducibility with Sullivan's pipeline.
#
# Usage:
#   bash scripts/align_klepikova.sh           # build index + align all
#   bash scripts/align_klepikova.sh index      # build Bowtie2 index only
#   bash scripts/align_klepikova.sh align      # align only (index must exist)
#   bash scripts/align_klepikova.sh align 8    # align with 8 threads

set -euo pipefail

DATADIR="data/klepikova"
TRIMDIR="${DATADIR}/trimmed"
BAMDIR="${DATADIR}/bam"
REFDIR="${DATADIR}/reference"
INDEXDIR="${REFDIR}/bowtie2_index"

THREADS="${2:-40}"

GENOME_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
GTF_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.59.gtf.gz"

# -------------------------------------------------------
# Download reference genome + GTF and build Bowtie2 index
# -------------------------------------------------------
build_index() {
    mkdir -p "${REFDIR}" "${INDEXDIR}"

    # Download TAIR10 genome; -k for retaining the compressed archive and don't have to re-download it if something goes wrong downstream
    GENOME_GZ="${REFDIR}/TAIR10_genome.fa.gz"
    GENOME_FA="${REFDIR}/TAIR10_genome.fa"
    if [[ ! -f "${GENOME_FA}" ]]; then
        echo "=== Downloading TAIR10 genome ==="
        curl -L -o "${GENOME_GZ}" "${GENOME_URL}"
        gunzip -k "${GENOME_GZ}" # 
        echo "  → ${GENOME_FA}"
    else
        echo "  TAIR10 genome already exists: ${GENOME_FA}"
    fi

    # Download TAIR10 GTF annotation; change to specify if you're using a different organism
    GTF_GZ="${REFDIR}/TAIR10.gtf.gz"
    GTF="${REFDIR}/TAIR10.gtf"
    if [[ ! -f "${GTF}" ]]; then
        echo "=== Downloading TAIR10 GTF annotation ==="
        curl -L -o "${GTF_GZ}" "${GTF_URL}"
        gunzip -k "${GTF_GZ}"
        echo "  → ${GTF}"
    else
        echo "  TAIR10 GTF already exists: ${GTF}"
    fi

    # Build Bowtie2 index
    if [[ ! -f "${INDEXDIR}/TAIR10.1.bt2" ]]; then
        echo "=== Building Bowtie2 index ==="
        bowtie2-build --threads "${THREADS}" "${GENOME_FA}" "${INDEXDIR}/TAIR10"
        echo "  → ${INDEXDIR}/TAIR10"
    else
        echo "  Bowtie2 index already exists: ${INDEXDIR}/TAIR10"
    fi

    echo "=== Reference setup complete ==="
}

# --------------------------------
# Align trimmed reads with TopHat2
# --------------------------------
align_reads() {
    if [[ ! -f "${INDEXDIR}/TAIR10.1.bt2" ]]; then
        echo "ERROR: Bowtie2 index not found at ${INDEXDIR}/TAIR10"
        echo "Run: bash scripts/align_klepikova.sh index"
        exit 1
    fi

    GTF="${REFDIR}/TAIR10.gtf"
    if [[ ! -f "${GTF}" ]]; then
        echo "ERROR: GTF not found at ${GTF}"
        exit 1
    fi

    mkdir -p "${BAMDIR}"
    # Parameters from Sullivan et al. 2019
    TRIMMED_FILES=("${TRIMDIR}"/*_trimmed.fastq.gz)
    NFILES=${#TRIMMED_FILES[@]}
    echo "=== Aligning ${NFILES} trimmed FASTQs with TopHat2 ==="
    echo "  Threads: ${THREADS}"
    echo "  Min intron: 50 bp"
    echo "  Max intron: 5000 bp"
    echo ""

    COUNT=0
    for FQ in "${TRIMMED_FILES[@]}"; do
        COUNT=$((COUNT + 1))
        BASENAME=$(basename "${FQ}" _trimmed.fastq.gz)
        OUTDIR="${BAMDIR}/${BASENAME}"
    # If the pipeline dies halfway through, re-running it picks up from the last unfinished sample rather than starting over
        if [[ -f "${OUTDIR}/accepted_hits.bam" ]]; then
            echo "  [${COUNT}/${NFILES}] SKIP ${BASENAME} (BAM exists)"
            continue
        fi

        echo "  [${COUNT}/${NFILES}] Aligning ${BASENAME}..."

        tophat2 \
            -p "${THREADS}" \
            -i 50 \
            -I 5000 \
            -G "${GTF}" \
            -o "${OUTDIR}" \
            "${INDEXDIR}/TAIR10" \
            "${FQ}"

        echo "  [${COUNT}/${NFILES}] Done: ${OUTDIR}/accepted_hits.bam"
    done

    echo ""
    echo "=== Alignment complete ==="
    echo "  BAM files: ${BAMDIR}/*/accepted_hits.bam"
}

# ----
# Main
# ----
case "${1:-all}" in
    index)
        build_index
        ;;
    align)
        align_reads
        ;;
    all)
        build_index
        align_reads
        ;;
    *)
        echo "Usage: bash scripts/align_klepikova.sh {index|align|all} [threads]"
        exit 1
        ;;
esac
