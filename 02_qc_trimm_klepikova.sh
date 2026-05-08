#!/usr/bin/env bash
# qc_trim_klepikova.sh
# ====================
# Quality check (FastQC) and adapter/quality trimming (Trimmomatic-0.32)
# on downloaded Klepikova FASTQs.
#
# Follows the published protocol from Klepikova et al. (2016):
#   - Remove residual adapters
#   - Remove low-quality sequences (Trimmomatic-0.32 default criteria)
#   - Discard reads below 36 bp
#
# Dependencies: FastQC, Trimmomatic (0.32 recommended)
#   conda install -c bioconda fastqc trimmomatic
#
# Usage:
#   bash scripts/qc_trim_klepikova.sh
#   bash scripts/qc_trim_klepikova.sh 8    # use 8 threads

set -euo pipefail

THREADS="${1:-40}"

DATADIR="data/klepikova"
FASTQDIR="${DATADIR}/fastqs"
TRIMDIR="${DATADIR}/trimmed"
QCDIR_PRE="${DATADIR}/fastqc_pretrim"
QCDIR_POST="${DATADIR}/fastqc_posttrim"

mkdir -p "${TRIMDIR}" "${QCDIR_PRE}" "${QCDIR_POST}"

# Locate Trimmomatic adapter file
#  Trimmomatic ships with adapter FASTA files bundled, but their path
#  varies between conda versions and system installs. Rather than hardcoding a path that breaks on other machines, the loop checks four common locations in order and takes the first one that actually exists. If nothing is found, the script warns but does not abort — it runs Trimmomatic without the ILLUMINACLIP step.
ADAPTER_FILE=""
for CANDIDATE in \
    "${CONDA_PREFIX:-/dev/null}/share/trimmomatic/adapters/TruSeq3-SE.fa" \
    "${CONDA_PREFIX:-/dev/null}/share/trimmomatic-0.39-2/adapters/TruSeq3-SE.fa" \
    "/usr/share/trimmomatic/adapters/TruSeq3-SE.fa" \
    "/usr/local/share/trimmomatic/adapters/TruSeq3-SE.fa"; do
    if [[ -f "${CANDIDATE}" ]]; then
        ADAPTER_FILE="${CANDIDATE}"
        break
    fi
done

if [[ -z "${ADAPTER_FILE}" ]]; then
    echo "WARNING: Could not find TruSeq3-SE.fa adapter file."
    echo "Set ADAPTER_FILE manually or install Trimmomatic via conda."
    echo "Proceeding without ILLUMINACLIP step..."
fi

# Count input files
FASTQS=("${FASTQDIR}"/*.fastq.gz)
NFILES=${#FASTQS[@]}
echo "=== QC + Trimming: ${NFILES} files, ${THREADS} threads ==="

COUNT=0
for FQ in "${FASTQS[@]}"; do
    COUNT=$((COUNT + 1))
    BASENAME=$(basename "${FQ}" .fastq.gz)
    echo ""
    echo "--- [${COUNT}/${NFILES}] ${BASENAME} ---"

    # Step 1: Pre-trim FastQC
    echo "  FastQC (pre-trim)..."
    fastqc --quiet --threads "${THREADS}" --outdir "${QCDIR_PRE}" "${FQ}"

    # Step 2: Trimmomatic
    echo "  Trimmomatic trimming..."
    TRIMMED="${TRIMDIR}/${BASENAME}_trimmed.fastq.gz"
    # allow up to 2 mismatches when finding adapters, and use thresholds of 30/10 for the two adapter-finding algorithm; replicating Sullivan et al 2019. This results in many reads being thrown out. Adjust for different parameters if so desired.
    if [[ -n "${ADAPTER_FILE}" ]]; then
        trimmomatic SE -threads "${THREADS}" \
            "${FQ}" "${TRIMMED}" \
            ILLUMINACLIP:"${ADAPTER_FILE}":2:30:10 \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36
    else
        trimmomatic SE -threads "${THREADS}" \
            "${FQ}" "${TRIMMED}" \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36
    fi

    # Step 3: Post-trim FastQC
    echo "  FastQC (post-trim)..."
    fastqc --quiet --threads "${THREADS}" --outdir "${QCDIR_POST}" "${TRIMMED}"

    echo "  Done: ${TRIMMED}"
done

echo ""
echo "=== Trimming complete ==="
echo "  Trimmed FASTQs: ${TRIMDIR}/"
echo "  Pre-trim QC:    ${QCDIR_PRE}/"
echo "  Post-trim QC:   ${QCDIR_POST}/"
