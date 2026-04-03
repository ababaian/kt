#!/bin/bash

# Fast family identifier counter using compressed input and optimized pipeline

# Check if filename is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <filename.zst>"
    echo "Example: $0 kt1-diamond-mar26_26-00001.pro.zst"
    exit 1
fi

input_file="$1"

# Check if file exists
if [ ! -f "$input_file" ]; then
    echo "Error: File '$input_file' not found!"
    exit 1
fi

echo "Fast counting family identifiers in $input_file..."
echo "Family_ID\tCount"
echo "-------------------"

# Optimized pipeline:
# 1. Decompress with zstd
# 2. Grep only lines with kt1.kt in column 6 (much faster pre-filter)
# 3. Cut to extract column 6 directly
# 4. Grep again to ensure kt1.kt; prefix (double-check after cut)
# 5. Cut to extract family ID (2nd field after split by ;)
# 6. Sort and count with uniq -c (highly optimized)
# 7. Format output with awk

zstd -dc "$input_file" | \
grep -F $'\tkt1.kt;' | \
cut -f6 | \
grep '^kt1\.kt;' | \
cut -d';' -f2 | \
sort | \
uniq -c | \
awk '{printf "%s\t%d\n", $2, $1; total += $1} END {print "-------------------"; printf "Total:\t%d\n", total}' | \
sort -k1,1

echo "Processing complete!"