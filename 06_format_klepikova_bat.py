#!/usr/bin/env python3
"""
format_klepikova_bat.py
=======================
Convert the aggregated Klepikova expression matrix into the final BAT
dataset format and install it into datasets/.

Reads:  data/klepikova/klepikova_expression_matrix.tsv
Writes: datasets/klepikova_tair10_expression.txt

The output file is tab-separated with:
  - Column 1: "taxa" (TAIR10 gene-level IDs, e.g. AT1G01010)
  - Columns 2+: log2(normalized counts + 1) per tissue
  - "NA" for genes below detection (value == 0.0)

Usage:
    python3 scripts/format_klepikova_bat.py
"""

import csv
import os
import re
import sys

DATADIR = "data/klepikova"
INPUT_FILE = os.path.join(DATADIR, "klepikova_expression_matrix.tsv")
OUTPUT_FILE = os.path.join("datasets", "klepikova_tair10_expression.txt")

# TAIR10 gene-level ID pattern: AT[1-5]G[0-9]{5}
TAIR10_PATTERN = re.compile(r"^AT[1-5]G\d{5}$")


def main():
    if not os.path.isfile(INPUT_FILE):
        print(f"ERROR: Input file not found: {INPUT_FILE}", file=sys.stderr)
        print("Run scripts/aggregate_klepikova.R first.", file=sys.stderr)
        sys.exit(1)

    print(f"=== Reading {INPUT_FILE} ===")

    rows = []
    with open(INPUT_FILE, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)

        # Ensure first column is "taxa"
        header[0] = "taxa"

        for row in reader:
            gene_id = row[0]

            # Strip isoform suffix if present (e.g., AT1G01010.1 -> AT1G01010)
            gene_id = re.sub(r"\.\d+$", "", gene_id)

            # Validate TAIR10 gene ID format
            if not TAIR10_PATTERN.match(gene_id):
                continue  # skip non-TAIR10 entries (organellar, etc.)

            # Replace 0.0 values with NA
            values = []
            for val in row[1:]:
                try:
                    num = float(val)
                    if num == 0.0:
                        values.append("NA")
                    else:
                        values.append(val)
                except ValueError:
                    values.append("NA")

            rows.append([gene_id] + values)

    print(f"  Genes passing filter: {len(rows)}")
    print(f"  Tissues: {', '.join(header[1:])}")

    # Check for duplicate gene IDs (shouldn't happen after R aggregation, but to be safe)
    seen = set()
    unique_rows = []
    for row in rows:
        if row[0] not in seen:
            seen.add(row[0])
            unique_rows.append(row)
    if len(unique_rows) < len(rows):
        print(f"  Deduplicated: {len(rows)} -> {len(unique_rows)} genes")
    rows = unique_rows

    # Write output
    os.makedirs(os.path.dirname(OUTPUT_FILE), exist_ok=True)
    print(f"=== Writing {OUTPUT_FILE} ===")

    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)

    print(f"  {len(rows)} genes x {len(header) - 1} tissues")
    print(f"  Output: {OUTPUT_FILE}")
    print("  Done.")


if __name__ == "__main__":
    main()
