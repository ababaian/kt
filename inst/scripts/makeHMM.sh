#!/bin/bash

# makeHMM.sh - Iterative HMM building and refinement script
# Usage: ./makeHMM.sh [options]
# This script builds HMM models iteratively using hhblits and allows manual inspection

set -e  # Exit on any error

# Configuration
INITIAL_MSA="K1.abg.msa.it0.fa"
TARGET_DB="kt0.preclust.fa"
MAX_ITERATIONS=5
E_VALUE_THRESHOLD=1e-3
MIN_COVERAGE=0.3
OUTPUT_DIR="hmm_iterations"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_step() {
    echo -e "${BLUE}[STEP]${NC} $1"
}

# Function to check required files and tools
check_requirements() {
    print_status "Checking requirements..."

    # Check input files
    if [[ ! -f "$INITIAL_MSA" ]]; then
        print_error "Initial MSA file not found: $INITIAL_MSA"
        exit 1
    fi

    if [[ ! -f "$TARGET_DB" ]]; then
        print_error "Target database file not found: $TARGET_DB"
        exit 1
    fi

    # Check required tools
    for tool in hhblits hhmake hhconsensus muscle mafft; do
        if ! command -v "$tool" &> /dev/null; then
            print_warning "$tool not found. Some features may not work."
        fi
    done

    print_status "Requirements check completed."
}

# Function to setup directory structure
setup_directories() {
    print_status "Setting up directory structure..."

    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR/iteration_logs"
    mkdir -p "$OUTPUT_DIR/msas"
    mkdir -p "$OUTPUT_DIR/hmms"
    mkdir -p "$OUTPUT_DIR/search_results"

    # Copy initial MSA to iteration 0
    cp "$INITIAL_MSA" "$OUTPUT_DIR/msas/K1.msa.it0.fa"

    print_status "Directory structure created in $OUTPUT_DIR"
}

# Function to build HMM from MSA
build_hmm() {
    local iteration=$1
    local msa_file="$OUTPUT_DIR/msas/K1.msa.it${iteration}.fa"
    local hmm_file="$OUTPUT_DIR/hmms/K1.hmm.it${iteration}"

    print_step "Building HMM for iteration $iteration..."

    if [[ ! -f "$msa_file" ]]; then
        print_error "MSA file not found: $msa_file"
        return 1
    fi

    # Build HMM using hhmake (HHsuite format)
    hhmake -M first -i "$msa_file" -o "$hmm_file.hhm" -name "K1_iteration_$iteration" 2>&1 | \
        tee "$OUTPUT_DIR/iteration_logs/hhmake_it${iteration}.log"

    # Also build HMMER format HMM
    if command -v hmmbuild &> /dev/null; then
        hmmbuild "$hmm_file.hmmer" "$msa_file" 2>&1 | \
            tee -a "$OUTPUT_DIR/iteration_logs/hmmbuild_it${iteration}.log"
    fi

    print_status "HMM built successfully: $hmm_file.hhm"

    # Generate HMM statistics
    echo "=== HMM Statistics for Iteration $iteration ===" >> "$OUTPUT_DIR/iteration_logs/hmm_stats.log"
    echo "MSA file: $msa_file" >> "$OUTPUT_DIR/iteration_logs/hmm_stats.log"
    echo "Number of sequences in MSA: $(grep -c '^>' "$msa_file")" >> "$OUTPUT_DIR/iteration_logs/hmm_stats.log"
    echo "MSA length: $(head -2 "$msa_file" | tail -1 | wc -c)" >> "$OUTPUT_DIR/iteration_logs/hmm_stats.log"
    echo "Timestamp: $(date)" >> "$OUTPUT_DIR/iteration_logs/hmm_stats.log"
    echo "" >> "$OUTPUT_DIR/iteration_logs/hmm_stats.log"
}

# Function to create HHblits database properly
create_hhblits_db() {
    local db_basename="${TARGET_DB%.*}_hhblits"
    local db_ffdata="${db_basename}_db.ffdata"
    local db_ffindex="${db_basename}_db.ffindex"
    local db_cs219_ffdata="${db_basename}_db_cs219.ffdata"
    local db_cs219_ffindex="${db_basename}_db_cs219.ffindex"

    # Set global variable for return value
    HHBLITS_DB_BASENAME="$db_basename"

    # Check if database already exists
    if [[ -f "$db_ffdata" && -f "$db_ffindex" && -f "$db_cs219_ffdata" && -f "$db_cs219_ffindex" ]]; then
        print_status "HHblits database already exists" >&2
        return 0
    fi

    print_status "Creating HHblits database from $TARGET_DB..." >&2

    # Step 1: Create ffindex database from FASTA
    if [[ ! -f "$db_ffdata" || ! -f "$db_ffindex" ]]; then
        print_status "Step 1: Creating ffindex database..." >&2
        ffindex_from_fasta -s "$db_ffdata" "$db_ffindex" "$TARGET_DB" \
            &> "$OUTPUT_DIR/iteration_logs/ffindex_build.log"

        if [[ ! -f "$db_ffdata" || ! -f "$db_ffindex" ]]; then
            print_error "Failed to create ffindex database" >&2
            return 1
        fi
        print_status "ffindex database created successfully" >&2
    fi

    # Step 2: Create cs219 translated database for prefiltering
    if [[ ! -f "$db_cs219_ffdata" || ! -f "$db_cs219_ffindex" ]]; then
        print_status "Step 2: Creating CS219 context state database for prefiltering..." >&2

        # Use cstranslate with correct parameters - need to use basename, not full paths
        local db_basename_only="${db_basename}_db"
        local cs219_basename_only="${db_basename}_db_cs219"

        print_status "Running cstranslate: $db_basename_only -> $cs219_basename_only" >&2

        /usr/libexec/hhsuite/cstranslate-plain \
            -f \
            -I fas \
            -x 0.3 \
            -c 4 \
            -i "$db_basename_only" \
            -o "$cs219_basename_only" \
            &> "$OUTPUT_DIR/iteration_logs/cstranslate.log"

        # Check if cstranslate succeeded
        if [[ ! -f "$db_cs219_ffdata" || ! -f "$db_cs219_ffindex" ]]; then
            print_warning "CS219 creation failed with initial parameters, trying simpler approach..." >&2

            # Try with minimal parameters
            /usr/libexec/hhsuite/cstranslate-plain \
                -f \
                -i "$db_basename_only" \
                -o "$cs219_basename_only" \
                &>> "$OUTPUT_DIR/iteration_logs/cstranslate.log"
        fi

        if [[ ! -f "$db_cs219_ffdata" || ! -f "$db_cs219_ffindex" ]]; then
            print_error "Failed to create CS219 database even with alternative approach" >&2
            return 1
        fi
        print_status "CS219 database created successfully" >&2
    fi

    print_status "HHblits database created successfully" >&2
    return 0
}

# Function to run hhblits search
run_hhblits_search() {
    local iteration=$1
    local msa_file="$OUTPUT_DIR/msas/K1.msa.it${iteration}.fa"
    local hmm_file="$OUTPUT_DIR/hmms/K1.hmm.it${iteration}.hhm"
    local search_out="$OUTPUT_DIR/search_results/hhblits_it${iteration}"

    print_step "Running hhblits search for iteration $iteration..."

    if [[ ! -f "$msa_file" ]]; then
        print_error "MSA file not found: $msa_file"
        return 1
    fi

    if [[ ! -f "$hmm_file" ]]; then
        print_error "HMM file not found: $hmm_file"
        return 1
    fi

    # Create database if needed
    if ! create_hhblits_db; then
        return 1
    fi
    local db_basename="$HHBLITS_DB_BASENAME"

    # Convert MSA to A3M format for hhblits input
    local query_a3m="$search_out.query.a3m"

    # Try to use reformat.pl for proper conversion
    if [[ -f "/usr/share/hhsuite/scripts/reformat.pl" ]]; then
        /usr/share/hhsuite/scripts/reformat.pl fas a3m "$msa_file" "$query_a3m" 2>/dev/null || {
            print_warning "reformat.pl failed, using simple copy"
            cp "$msa_file" "$query_a3m"
        }
    else
        # Simple copy as fallback
        cp "$msa_file" "$query_a3m"
    fi

    print_status "Running hhblits with query MSA against custom database..."

    # Run hhblits search using 1 iteration (recommended for custom genome databases)
    hhblits \
        -i "$query_a3m" \
        -d "${db_basename}_db" \
        -o "$search_out.hhr" \
        -oa3m "$search_out.a3m" \
        -e "$E_VALUE_THRESHOLD" \
        -cov "$MIN_COVERAGE" \
        -n 1 \
        -cpu 4 \
        -v 2 \
        2>&1 | tee "$OUTPUT_DIR/iteration_logs/hhblits_it${iteration}.log"

    # Clean up temporary query file
    rm -f "$query_a3m"

    # Check if search was successful
    if [[ ! -f "$search_out.hhr" ]]; then
        print_error "hhblits search failed - no output file generated"
        return 1
    fi

    # Check for common error patterns
    if grep -q "ERROR.*FFindexDatabase" "$search_out.hhr" 2>/dev/null; then
        print_error "Database format error detected. Please check database creation."
        return 1
    fi

    # Extract significant hits from the results
    if [[ -f "$search_out.hhr" && -s "$search_out.hhr" ]]; then
        # Look for the hits section in hhblits output
        if grep -q "No 1" "$search_out.hhr" 2>/dev/null; then
            # Found hits - extract them
            # hhblits output has numbered hits like "No 1", "No 2", etc.
            grep -E "^No [0-9]+" "$search_out.hhr" | head -100 > "$search_out.hits.txt"

            # If that doesn't work, try alternative patterns
            if [[ ! -s "$search_out.hits.txt" ]]; then
                grep -E "^[[:space:]]*[0-9]+[[:space:]]+" "$search_out.hhr" | head -100 > "$search_out.hits.txt"
            fi
        else
            echo "0" > "$search_out.hits.txt"
            print_warning "No significant hits found in hhblits output"
            # Create minimal a3m file for next iteration
            head -2 "$msa_file" > "$search_out.a3m"
        fi
    else
        print_error "hhblits search failed - empty or missing output"
        return 1
    fi

    print_status "hhblits search completed. Results: $search_out.hhr"

    # Count hits
    local hit_count=$(wc -l < "$search_out.hits.txt")
    print_status "Found $hit_count significant hits"

    # Log some statistics
    if [[ -f "$search_out.a3m" && -s "$search_out.a3m" ]]; then
        local a3m_seqs=$(grep -c "^>" "$search_out.a3m")
        print_status "Output MSA contains $a3m_seqs sequences"
    fi

    return 0
}

# Function to create new MSA from search results
create_new_msa() {
    local iteration=$1
    local search_a3m="$OUTPUT_DIR/search_results/hhblits_it${iteration}.a3m"
    local new_msa="$OUTPUT_DIR/msas/K1.msa.it$((iteration+1)).fa"

    print_step "Creating new MSA for iteration $((iteration+1))..."

    if [[ ! -f "$search_a3m" ]]; then
        print_error "Search results file not found: $search_a3m"
        return 1
    fi

    # Check if a3m file has content
    if [[ ! -s "$search_a3m" ]]; then
        print_warning "Empty search results file. No new MSA will be created."
        return 1
    fi

    # Convert a3m to fasta format
    # a3m format uses lowercase letters for insertions relative to consensus
    # Remove lowercase letters (insertions) and convert to standard FASTA
    perl -pe '
        if (/^>/) {
            s/^>(.+)/>iteration_'$((iteration+1))'_$1/;
        } else {
            s/[a-z]//g;  # Remove lowercase insertions
        }
    ' "$search_a3m" > "$new_msa.tmp"

    # Filter sequences by quality
    awk '
    BEGIN { RS=">"; ORS="" }
    NR==1 { next }
    {
        # Extract header and sequence
        split($0, parts, "\n")
        header = parts[1]
        seq = ""
        for(i=2; i<=length(parts); i++) {
            if(parts[i] != "") seq = seq parts[i]
        }

        # Remove any remaining whitespace
        gsub(/[ \t\r]/, "", seq)

        # Calculate gap statistics
        gap_count = gsub(/-/, "-", seq)
        total_len = length(seq)

        # Filter: minimum length 50, maximum 80% gaps
        if(total_len > 50 && (total_len == 0 || gap_count/total_len < 0.8)) {
            print ">" header "\n" seq "\n"
        }
    }' "$new_msa.tmp" > "$new_msa"

    rm -f "$new_msa.tmp"

    # Check if we have any sequences
    if [[ ! -s "$new_msa" ]]; then
        print_warning "No sequences passed quality filters."

        # Try with more lenient filters
        print_status "Trying with more lenient filters..."
        perl -pe '
            if (/^>/) {
                s/^>(.+)/>iteration_'$((iteration+1))'_$1/;
            } else {
                s/[a-z]//g;
            }
        ' "$search_a3m" | awk '
        BEGIN { RS=">"; ORS="" }
        NR==1 { next }
        {
            split($0, parts, "\n")
            header = parts[1]
            seq = ""
            for(i=2; i<=length(parts); i++) {
                if(parts[i] != "") seq = seq parts[i]
            }
            gsub(/[ \t\r]/, "", seq)

            if(length(seq) > 20) {  # Very lenient filter
                print ">" header "\n" seq "\n"
            }
        }' > "$new_msa"
    fi

    if [[ ! -s "$new_msa" ]]; then
        print_error "Could not create valid MSA from search results"
        return 1
    fi

    local seq_count=$(grep -c '^>' "$new_msa")
    print_status "New MSA created with $seq_count sequences: $new_msa"

    # Optionally refine MSA with muscle if available and reasonable number of sequences
    if command -v muscle &> /dev/null && [[ $seq_count -le 500 && $seq_count -ge 3 ]]; then
        print_status "Refining MSA with MUSCLE..."
        if muscle -in "$new_msa" -out "$new_msa.refined" -maxiters 2 &>/dev/null; then
            if [[ -s "$new_msa.refined" ]]; then
                mv "$new_msa.refined" "$new_msa"
                print_status "MSA refined successfully"
            fi
        else
            print_warning "MUSCLE refinement failed, using original MSA"
        fi
    fi

    return 0
}

# Function to display iteration summary
display_summary() {
    local iteration=$1
    local msa_file="$OUTPUT_DIR/msas/K1.msa.it${iteration}.fa"
    local search_results="$OUTPUT_DIR/search_results/hhblits_it${iteration}.hits.txt"
    local search_hhr="$OUTPUT_DIR/search_results/hhblits_it${iteration}.hhr"

    echo ""
    echo "=================================="
    echo "    ITERATION $iteration SUMMARY"
    echo "=================================="
    echo "MSA file: $msa_file"
    if [[ -f "$msa_file" ]]; then
        local seq_count=$(grep -c '^>' "$msa_file")
        echo "Sequences in MSA: $seq_count"

        if [[ $seq_count -gt 0 ]]; then
            local avg_len=$(awk '
            /^>/ {
                if(seq) {
                    print length(seq)
                }
                seq = ""
                next
            }
            {
                gsub(/[ \t]/, "", $0)
                seq = seq $0
            }
            END {
                if(seq) {
                    print length(seq)
                }
            }' "$msa_file" | awk '{sum+=$1; n++} END{if(n>0) print int(sum/n); else print 0}')
            echo "Average sequence length: $avg_len"
        fi
    fi

    # Display search results
    if [[ -f "$search_hhr" ]]; then
        local hit_count=$(wc -l < "$search_results" 2>/dev/null || echo "0")
        echo "Significant hits found: $hit_count"

        if [[ $hit_count -gt 0 && -s "$search_results" ]]; then
            echo ""
            echo "Top 5 hits from hhblits:"
            head -5 "$search_results" | nl
        fi

        # Show search statistics from hhblits output
        if grep -q "Done\|Searching" "$search_hhr" 2>/dev/null; then
            echo ""
            echo "Search completed successfully"

            # Try to extract some basic statistics
            if grep -q "Query" "$search_hhr"; then
                echo "Query info:"
                grep "Query" "$search_hhr" | head -2
            fi
        fi
    fi

    # Show HMM statistics
    local hmm_file="$OUTPUT_DIR/hmms/K1.hmm.it${iteration}.hhm"
    if [[ -f "$hmm_file" ]]; then
        echo ""
        echo "HMM file: $hmm_file"
        if grep -q "LENG" "$hmm_file"; then
            echo "HMM length: $(grep "LENG" "$hmm_file" | awk '{print $2}')"
        fi
        if grep -q "NEFF" "$hmm_file"; then
            echo "Effective sequences: $(grep "NEFF" "$hmm_file" | awk '{print $2}')"
        fi
    fi

    echo "=================================="
    echo ""
}

# Function to prompt for continuation
prompt_continue() {
    local iteration=$1

    echo ""
    print_warning "Please inspect the MSA file: $OUTPUT_DIR/msas/K1.msa.it$((iteration+1)).fa"
    echo ""
    echo "You can:"
    echo "1. View the MSA: less $OUTPUT_DIR/msas/K1.msa.it$((iteration+1)).fa"
    echo "2. Check alignment quality with: awk '/^>/{print NR, \$0}' $OUTPUT_DIR/msas/K1.msa.it$((iteration+1)).fa | head"
    echo "3. Compare with previous iteration: diff $OUTPUT_DIR/msas/K1.msa.it${iteration}.fa $OUTPUT_DIR/msas/K1.msa.it$((iteration+1)).fa"
    echo ""

    while true; do
        read -p "Continue to next iteration? (y/n/s=skip refinement/q=quit): " choice
        case $choice in
            [Yy]* ) return 0;;
            [Ss]* ) return 2;;
            [Nn]* ) return 1;;
            [Qq]* ) exit 0;;
            * ) echo "Please answer y(es), n(o), s(kip), or q(uit).";;
        esac
    done
}

# Main iteration function
run_iteration() {
    local iteration=$1

    echo ""
    print_step "Starting iteration $iteration"
    echo "==============================="

    # Build HMM
    if ! build_hmm "$iteration"; then
        return 1
    fi

    # Run search
    if ! run_hhblits_search "$iteration"; then
        return 1
    fi

    # Create new MSA (unless it's the last allowed iteration)
    if [[ $iteration -lt $MAX_ITERATIONS ]]; then
        if ! create_new_msa "$iteration"; then
            return 1
        fi
    fi

    # Display summary
    display_summary "$iteration"

    return 0
}

# Main function
main() {
    print_status "Starting iterative HMM building process"
    print_status "Initial MSA: $INITIAL_MSA"
    print_status "Target database: $TARGET_DB"
    print_status "Max iterations: $MAX_ITERATIONS"
    print_status "E-value threshold: $E_VALUE_THRESHOLD"
    echo ""

    # Check requirements and setup
    check_requirements
    setup_directories

    # Record run parameters
    {
        echo "Run started: $(date)"
        echo "Initial MSA: $INITIAL_MSA"
        echo "Target DB: $TARGET_DB"
        echo "Max iterations: $MAX_ITERATIONS"
        echo "E-value threshold: $E_VALUE_THRESHOLD"
        echo "Min coverage: $MIN_COVERAGE"
        echo ""
    } > "$OUTPUT_DIR/iteration_logs/run_parameters.log"

    # Run iterations
    for iteration in $(seq 0 $MAX_ITERATIONS); do
        if ! run_iteration "$iteration"; then
            print_error "Iteration $iteration failed"
            exit 1
        fi

        # Check if this was the last iteration
        if [[ $iteration -eq $MAX_ITERATIONS ]]; then
            print_status "Maximum iterations reached ($MAX_ITERATIONS)"
            break
        fi

        # Check if new MSA exists for next iteration
        next_msa="$OUTPUT_DIR/msas/K1.msa.it$((iteration+1)).fa"
        if [[ ! -f "$next_msa" ]]; then
            print_warning "No MSA generated for next iteration. Stopping."
            break
        fi

        # Prompt user for continuation
        prompt_continue "$iteration"
        result=$?

        if [[ $result -eq 1 ]]; then
            print_status "Stopping at user request after iteration $iteration"
            break
        elif [[ $result -eq 2 ]]; then
            print_status "Skipping MSA refinement for next iteration"
        fi

        # Check convergence (optional)
        if [[ $iteration -gt 0 ]]; then
            prev_hmm="$OUTPUT_DIR/hmms/K1.hmm.it$((iteration-1)).hhm"
            curr_hmm="$OUTPUT_DIR/hmms/K1.hmm.it${iteration}.hhm"
            if [[ -f "$prev_hmm" && -f "$curr_hmm" ]]; then
                # Simple convergence check could be added here
                print_status "Convergence check could be implemented here"
            fi
        fi
    done

    # Final summary
    echo ""
    print_status "HMM building process completed!"
    print_status "Results are in: $OUTPUT_DIR/"
    print_status "Final HMM files are in: $OUTPUT_DIR/hmms/"
    print_status "Final MSAs are in: $OUTPUT_DIR/msas/"

    # Create final summary
    {
        echo "Run completed: $(date)"
        echo "Total iterations completed: $(ls $OUTPUT_DIR/hmms/*.hhm | wc -l)"
        echo ""
        echo "Final file listing:"
        ls -la "$OUTPUT_DIR"/hmms/
        echo ""
        ls -la "$OUTPUT_DIR"/msas/
    } >> "$OUTPUT_DIR/iteration_logs/run_parameters.log"
}

# Handle command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -e|--evalue)
            E_VALUE_THRESHOLD="$2"
            shift 2
            ;;
        -c|--coverage)
            MIN_COVERAGE="$2"
            shift 2
            ;;
        -i|--iterations)
            MAX_ITERATIONS="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  -e, --evalue FLOAT     E-value threshold (default: $E_VALUE_THRESHOLD)"
            echo "  -c, --coverage FLOAT   Minimum coverage (default: $MIN_COVERAGE)"
            echo "  -i, --iterations INT   Maximum iterations (default: $MAX_ITERATIONS)"
            echo "  -o, --output DIR       Output directory (default: $OUTPUT_DIR)"
            echo "  -h, --help            Show this help message"
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

# Run main function
main