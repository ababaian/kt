#!/bin/bash

# =============================================================================
# Simple HMMER-based Protein Sequence Search Pipeline
# =============================================================================
# This script uses HMMER to build an HMM, search a database, and create
# an output MSA with significant matches (limited for performance)
# =============================================================================

set -e

# Configuration
QUERY_FILE="M28.it1.alpha_trim.fa"
DATABASE_FILE="kt1.M28.ORF50.id99s.fa"
E_VALUE_CUTOFF=0.0001
MAX_HITS_TO_ALIGN=100000
CPU_CORES=16

# Directory setup
BASE_DIR="."
INPUT_DIR="$BASE_DIR/input"
OUTPUT_DIR="$BASE_DIR/output"
RESULTS_DIR="$OUTPUT_DIR/hmmer_simple_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS_DIR"

log() {
    echo "[$(date +'%H:%M:%S')] $1" | tee -a "$RESULTS_DIR/pipeline.log"
}

log "Starting simple HMMER pipeline"
log "Results directory: $RESULTS_DIR"

# Step 1: Build HMM from query MSA
log "Building HMM from query MSA..."
cp "$INPUT_DIR/$QUERY_FILE" "$RESULTS_DIR/"
cp "$INPUT_DIR/$DATABASE_FILE" "$RESULTS_DIR/"

hmmbuild --amino --cpu "$CPU_CORES" "$RESULTS_DIR/query.hmm" "$RESULTS_DIR/$QUERY_FILE" > "$RESULTS_DIR/hmmbuild.log"

model_name=$(grep "^NAME" "$RESULTS_DIR/query.hmm" | awk '{print $2}')
model_length=$(grep "^LENG" "$RESULTS_DIR/query.hmm" | awk '{print $2}')
log "✓ HMM built: $model_name (length: $model_length)"

# Step 2: Search database
log "Searching database..."
hmmsearch \
    --cpu "$CPU_CORES" \
    -E "$E_VALUE_CUTOFF" \
    --tblout "$RESULTS_DIR/hits.tbl" \
    -o "$RESULTS_DIR/search.out" \
    "$RESULTS_DIR/query.hmm" \
    "$RESULTS_DIR/$DATABASE_FILE" > "$RESULTS_DIR/hmmsearch.log"

total_hits=$(grep -c -v "^#" "$RESULTS_DIR/hits.tbl" 2>/dev/null || echo 0)
significant_hits=$(awk -v cutoff="$E_VALUE_CUTOFF" '!/^#/ && $5 < cutoff' "$RESULTS_DIR/hits.tbl" | wc -l)

log "✓ Search complete: $total_hits hits found, $significant_hits significant"

# Step 3: Extract top significant hits
log "Extracting top $MAX_HITS_TO_ALIGN significant hits..."
awk -v cutoff="$E_VALUE_CUTOFF" '!/^#/ && $5 < cutoff {print $1}' "$RESULTS_DIR/hits.tbl" | head -"$MAX_HITS_TO_ALIGN" > "$RESULTS_DIR/top_hit_ids.txt"

hits_to_align=$(wc -l < "$RESULTS_DIR/top_hit_ids.txt")
log "Aligning top $hits_to_align hits"

if [[ $hits_to_align -eq 0 ]]; then
    log "No significant hits found"
    exit 0
fi

# Extract hit sequences
awk '
BEGIN {
    while ((getline id < "'$RESULTS_DIR'/top_hit_ids.txt") > 0) {
        hit_ids[id] = 1
    }
    close("'$RESULTS_DIR'/top_hit_ids.txt")
}
/^>/ {
    seq_id = substr($1, 2)
    if (seq_id in hit_ids) {
        print_seq = 1
        print $0
    } else {
        print_seq = 0
    }
}
!/^>/ && print_seq {
    print $0
}
' "$RESULTS_DIR/$DATABASE_FILE" > "$RESULTS_DIR/hit_sequences.fa"

extracted=$(grep -c "^>" "$RESULTS_DIR/hit_sequences.fa" 2>/dev/null || echo 0)
log "✓ Extracted $extracted hit sequences"

# Step 4: Prepare sequences for alignment (remove gaps from query)
log "Creating alignment-ready sequences..."

awk '
/^>/ { print $0 }
!/^>/ { gsub(/-/, ""); if (length($0) > 0) print $0 }
' "$RESULTS_DIR/$QUERY_FILE" > "$RESULTS_DIR/query_ungapped.fa"

cat "$RESULTS_DIR/query_ungapped.fa" "$RESULTS_DIR/hit_sequences.fa" > "$RESULTS_DIR/all_sequences.fa"

seq_count=$(grep -c "^>" "$RESULTS_DIR/all_sequences.fa")
log "Combined $seq_count sequences for alignment"

# Step 5: Create output MSA
log "Creating output MSA..."
hmmalign --amino --outformat Stockholm "$RESULTS_DIR/query.hmm" "$RESULTS_DIR/all_sequences.fa" > "$RESULTS_DIR/output_msa.sto" 2> "$RESULTS_DIR/hmmalign.log"

hmmalign --amino --outformat FASTA "$RESULTS_DIR/query.hmm" "$RESULTS_DIR/all_sequences.fa" > "$RESULTS_DIR/output_msa.fa" 2>> "$RESULTS_DIR/hmmalign.log"

msa_seqs=$(grep -c "^>" "$RESULTS_DIR/output_msa.fa" 2>/dev/null || echo 0)
log "✓ Output MSA created with $msa_seqs sequences"

# Step 6: Create summary
log "Creating summary..."
cat > "$RESULTS_DIR/SUMMARY.txt" << EOF
HMMER Search Summary
===================
Date: $(date)
Query: $QUERY_FILE
Database: $DATABASE_FILE
E-value cutoff: $E_VALUE_CUTOFF
Max hits aligned: $MAX_HITS_TO_ALIGN

Results:
- Total hits: $total_hits
- Significant hits: $significant_hits
- Hits in final MSA: $msa_seqs

Top 10 Hits:
$(awk '!/^#/ {print NR ". " $1 " (E=" $5 ", Score=" $6 ")"}' "$RESULTS_DIR/hits.tbl" | head -10)

Files Generated:
- query.hmm           HMM model built from query MSA
- hits.tbl            All search hits in tabular format
- search.out          Full HMMER search output
- hit_sequences.fa    Extracted significant hit sequences
- output_msa.fa       Final multiple sequence alignment (FASTA)
- output_msa.sto      Final multiple sequence alignment (Stockholm)

Quick inspection commands:
  cat SUMMARY.txt
  less output_msa.fa
  less search.out
  awk '\$5 < $E_VALUE_CUTOFF' hits.tbl
EOF

log "Pipeline complete!"
log "Summary: cat $RESULTS_DIR/SUMMARY.txt"
log "Output MSA: $RESULTS_DIR/output_msa.fa"