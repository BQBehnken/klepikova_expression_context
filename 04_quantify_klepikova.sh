#!/usr/bin/env bash
# quantify_klepikova.sh
# =====================
# Count reads per gene from TopHat2 BAM files using featureCounts (subread).
# This produces a single raw count matrix across all samples.
#
# Dependencies: featureCounts (part of the subread package)
#   conda install -c bioconda subread
#
# Usage:
#   bash scripts/quantify_klepikova.sh        # default 4 threads
#   bash scripts/quantify_klepikova.sh 8      # use 8 threads

set -euo pipefail

THREADS="${1:-40}"

DATADIR="data/klepikova"
BAMDIR="${DATADIR}/bam"
COUNTDIR="${DATADIR}/counts"
GTF="${DATADIR}/reference/TAIR10.gtf"

mkdir -p "${COUNTDIR}"

# Verify inputs
if [[ ! -f "${GTF}" ]]; then
    echo "ERROR: GTF not found at ${GTF}"
    echo "Run: bash scripts/align_klepikova.sh index"
    exit 1
fi

# Collect all BAM files
BAM_FILES=()
for BAMPATH in "${BAMDIR}"/*/accepted_hits.bam; do
    if [[ -f "${BAMPATH}" ]]; then
        BAM_FILES+=("${BAMPATH}")
    fi
done

NBAMS=${#BAM_FILES[@]}
if [[ "${NBAMS}" -eq 0 ]]; then
    echo "ERROR: No BAM files found in ${BAMDIR}/*/accepted_hits.bam"
    echo "Run: bash scripts/align_klepikova.sh"
    exit 1
fi

echo "=== featureCounts: ${NBAMS} BAM files, ${THREADS} threads ==="
echo "  GTF: ${GTF}"
echo "  Output: ${COUNTDIR}/featurecounts_raw.txt"
echo ""

# Run featureCounts
# -s 0 : Klepikova used a TruSeq RNA Sample Prep Kits v2 kit, which produces unstranded libraries
# -t exon : ignores introns
# -g gene_id : summarize at gene level
featureCounts \
    -T "${THREADS}" \
    -s 0 \
    -t exon \
    -g gene_id \
    -a "${GTF}" \
    -o "${COUNTDIR}/featurecounts_raw.txt" \
    "${BAM_FILES[@]}"

echo ""
echo "=== Quantification complete ==="
echo "  Raw counts: ${COUNTDIR}/featurecounts_raw.txt"
echo "  Summary:    ${COUNTDIR}/featurecounts_raw.txt.summary"
echo ""
echo "Next step: bash scripts/aggregate_klepikova.R"
