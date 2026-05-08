#!/usr/bin/env bash
# download_klepikova.sh
# =====================
# Two-stage download of Klepikova et al. 2016 RNAseq data from NCBI SRA.
# Stage 1: Fetch metadata for both BioProjects.
# Stage 2: Filter to target tissues and download only those FASTQs.
# Stages can be run one by one, or both at the same time. 
#
# Dependencies: NCBI E-utilities (esearch, efetch), SRA Toolkit (prefetch, fasterq-dump)
# Install E-utilities: https://www.ncbi.nlm.nih.gov/books/NBK179288/
# Install SRA Toolkit: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
#
# How to use: Run each one sequentially in the CLI; open the CSV generated from stage 1, verify the labels look right before committing to a large download.
#   bash scripts/download_klepikova.sh metadata   # Stage 1
#   bash scripts/download_klepikova.sh download    # Stage 2 (after inspecting metadata)
#   bash scripts/download_klepikova.sh all         # Both stages sequentially

set -euo pipefail

DATADIR="data/klepikova"
FASTQDIR="${DATADIR}/fastqs"
METADIR="${DATADIR}/metadata"

BIOPROJECTS=("PRJNA314076" "PRJNA324514") # sourced from Klepikova 2016; change the targets for different desired tissues.

# Five representative tissues selected for the default BAT dataset; change the targets for different desired tissues.
TISSUES=("seed" "root" "stem" "petal" "carpel" "leaf")

# -----------------------
# Stage 1: Fetch metadata
# -----------------------
fetch_metadata() {
    mkdir -p "${METADIR}"

    echo "=== Stage 1: Fetching SRA run metadata ==="

    for PROJ in "${BIOPROJECTS[@]}"; do
        echo "  Querying ${PROJ}..."
        esearch -db sra -query "${PROJ}" \
            | efetch -format runinfo \
            > "${METADIR}/runinfo_${PROJ}.csv"
        echo "  → ${METADIR}/runinfo_${PROJ}.csv"
    done

    # Merge and deduplicate on SRR accession (column 1), keep header from first file
    head -1 "${METADIR}/runinfo_PRJNA314076.csv" > "${METADIR}/klepikova_runinfo.csv"
    for PROJ in "${BIOPROJECTS[@]}"; do
        tail -n +2 "${METADIR}/runinfo_${PROJ}.csv"
    done | sort -t',' -k1,1 -u >> "${METADIR}/klepikova_runinfo.csv"
    # drops any runs that appears in both projects

    TOTAL=$(tail -n +2 "${METADIR}/klepikova_runinfo.csv" | wc -l | tr -d ' ')
    echo ""
    echo "=== Merged metadata: ${TOTAL} total runs ==="
    echo "=== Saved to: ${METADIR}/klepikova_runinfo.csv ==="
    echo ""
    echo "Next step: inspect the CSV to verify tissue labels, then run:"
    echo "  bash scripts/download_klepikova.sh download"
}

# -----------------------------------------------------
# Stage 2: Filter to target tissues and download FASTQs
# -----------------------------------------------------
download_fastqs() {
    if [[ ! -f "${METADIR}/klepikova_runinfo.csv" ]]; then
        echo "ERROR: Metadata not found. Run Stage 1 first:"
        echo "  bash scripts/download_klepikova.sh metadata"
        exit 1
    fi

    mkdir -p "${FASTQDIR}"

    echo "=== Stage 2: Filtering for target tissues ==="

    # Build grep pattern from the tissue list (case-insensitive match on SampleName/LibraryName columns)
    PATTERN=$(IFS='|'; echo "${TISSUES[*]}")

    # Filter runinfo to selected tissues, extract SRR accessions (case-insensitive)
    grep -iE "${PATTERN}" "${METADIR}/klepikova_runinfo.csv" \
        | cut -d',' -f1 \
        | grep -E '^[SDE]RR' \
        > "${METADIR}/selected_srrs.txt" || true

    NRUNS=$(wc -l < "${METADIR}/selected_srrs.txt" | tr -d ' ')
    echo "  Found ${NRUNS} runs matching tissues: ${TISSUES[*]}"

    if [[ "${NRUNS}" -eq 0 ]]; then
        echo "ERROR: No SRR accessions matched the tissue filter."
        echo "Check tissue labels in ${METADIR}/klepikova_runinfo.csv"
        echo "and adjust the TISSUES array in this script if needed."
        exit 1
    fi

    echo ""
    echo "  SRR accessions to download:"
    cat "${METADIR}/selected_srrs.txt" | sed 's/^/    /'
    echo ""

    # Download each SRR
    COUNT=0
    while IFS= read -r SRR; do
        COUNT=$((COUNT + 1))
        echo "  [${COUNT}/${NRUNS}] Downloading ${SRR}..."

        prefetch "${SRR}" --max-size 50G

        fasterq-dump "${SRR}" --outdir "${FASTQDIR}" --progress

        # Compress single-end output
        if [[ -f "${FASTQDIR}/${SRR}.fastq" ]]; then
            gzip "${FASTQDIR}/${SRR}.fastq"
        fi
        # Handle paired files if any exist (Klepikova data is mostly single-end)
        for f in "${FASTQDIR}/${SRR}"_*.fastq; do
            [[ -f "$f" ]] && gzip "$f"
        done

        echo "  [${COUNT}/${NRUNS}] Done: ${SRR}"
    done < "${METADIR}/selected_srrs.txt"

    echo ""
    echo "=== Download complete: ${COUNT} runs ==="
    echo "=== FASTQs in: ${FASTQDIR}/ ==="
    ls -lh "${FASTQDIR}"/*.fastq.gz 2>/dev/null || echo "(no .fastq.gz files found)"
}

# ----
# Main
# ----
case "${1:-all}" in
    metadata)
        fetch_metadata
        ;;
    download)
        download_fastqs
        ;;
    all)
        fetch_metadata
        download_fastqs
        ;;
    *)
        echo "Usage: bash scripts/download_klepikova.sh {metadata|download|all}"
        exit 1
        ;;
esac
