#!/bin/bash
#SBATCH --job-name=AF3_homo_controller
#SBATCH --output=/home/kamounlab_gmail_com/slurm_logs/af3_homo_controller_%j.out
#SBATCH --error=/home/kamounlab_gmail_com/slurm_logs/af3_homo_controller_%j.err
#SBATCH --partition=datac3
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00

# ─────────────────────────────────────────────────────────────────────────────
# AF3 Homomeric Controller Job
# ─────────────────────────────────────────────────────────────────────────────
# This job runs on a CPU node and submits chunks one at a time.
# It waits for each chunk to complete before submitting the next.
# Only ~1000 jobs in queue at any time - avoids MaxSubmitJobs limits.
#
# Usage:
#   sbatch AF3_homomeric_controller.sh --manifest <file> --input <dir> \
#          --output <dir> --batch <n> --num-protomers <n> [options]
#
# Required:
#   --manifest <file>     Manifest TSV file
#   --input <dir>         Input directory
#   --output <dir>        Output directory
#   --batch <n>           Job name
#   --num-protomers <n>   Number of protomers
#
# Optional:
#   --ligand <spec>       Ligand specification
#   --num-ligands <n>     Number of ligands
#   --seeds <n...>        Seeds (default: 1)
#   --start-chunk <n>     Start from chunk (default: 1)
#   --end-chunk <n>       End at chunk (default: all)
#   --run-ipsae <bool>    Run ipSAE (default: true)
#   --run-scores <bool>   Extract scores (default: true)
#
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# Defaults
MAX_ARRAY_SIZE=1001
SEEDS="1"
START_CHUNK=1
END_CHUNK=""
CHECK_INTERVAL=30  # seconds between completion checks
SCRIPT_DIR="/home/kamounlab_gmail_com/inference_homomer_scripts"
RUN_IPSAE="true"
RUN_SCORES="true"

# ─────────────────────────────────────────────────────────────────────────────
# Parse arguments
# ─────────────────────────────────────────────────────────────────────────────
MANIFEST_FILE=""
INPUT_DIR=""
OUTPUT_DIR=""
BATCH_NAME=""
NUM_PROTOMERS=""
LIGAND=""
NUM_LIGANDS=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --manifest)
            MANIFEST_FILE="$2"
            shift 2
            ;;
        --input)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --batch)
            BATCH_NAME="$2"
            shift 2
            ;;
        --num-protomers)
            NUM_PROTOMERS="$2"
            shift 2
            ;;
        --ligand)
            LIGAND="$2"
            shift 2
            ;;
        --num-ligands)
            NUM_LIGANDS="$2"
            shift 2
            ;;
        --seeds)
            SEEDS=""
            shift
            while [[ $# -gt 0 && "$1" =~ ^[0-9]+$ ]]; do
                SEEDS="$SEEDS $1"
                shift
            done
            SEEDS="${SEEDS# }"
            ;;
        --start-chunk)
            START_CHUNK="$2"
            shift 2
            ;;
        --end-chunk)
            END_CHUNK="$2"
            shift 2
            ;;
        --script-dir)
            SCRIPT_DIR="$2"
            shift 2
            ;;
        --run-ipsae)
            RUN_IPSAE="$2"
            shift 2
            ;;
        --run-scores)
            RUN_SCORES="$2"
            shift 2
            ;;
        --help|-h)
            head -40 "$0" | tail -35
            exit 0
            ;;
        -*)
            echo "ERROR: Unknown option: $1"
            exit 1
            ;;
        *)
            echo "ERROR: Unexpected argument: $1"
            exit 1
            ;;
    esac
done

# ─────────────────────────────────────────────────────────────────────────────
# Set script paths (after argument parsing so --script-dir can override)
# ─────────────────────────────────────────────────────────────────────────────
JOB_SCRIPT="${SCRIPT_DIR}/AF3_homomeric_job.sh"
GENERATE_SCRIPT="${SCRIPT_DIR}/generate_af3_homomeric.py"
IPSAE_SCRIPT="${SCRIPT_DIR}/ipsae.py"
SCORE_SCRIPT="${SCRIPT_DIR}/extract_af3_scores.py"

# ─────────────────────────────────────────────────────────────────────────────
# Validate
# ─────────────────────────────────────────────────────────────────────────────
MISSING=""
[[ -z "$MANIFEST_FILE" ]] && MISSING="$MISSING --manifest"
[[ -z "$INPUT_DIR" ]] && MISSING="$MISSING --input"
[[ -z "$OUTPUT_DIR" ]] && MISSING="$MISSING --output"
[[ -z "$BATCH_NAME" ]] && MISSING="$MISSING --batch"
[[ -z "$NUM_PROTOMERS" ]] && MISSING="$MISSING --num-protomers"

if [[ -n "$MISSING" ]]; then
    echo "ERROR: Missing required arguments:$MISSING"
    exit 1
fi

[[ ! -f "$MANIFEST_FILE" ]] && echo "ERROR: Manifest not found: $MANIFEST_FILE" && exit 1
[[ ! -d "$INPUT_DIR" ]] && echo "ERROR: Input dir not found: $INPUT_DIR" && exit 1
[[ ! -f "$JOB_SCRIPT" ]] && echo "ERROR: Job script not found: $JOB_SCRIPT" && exit 1
[[ ! -f "$GENERATE_SCRIPT" ]] && echo "ERROR: Generate script not found: $GENERATE_SCRIPT" && exit 1

# Convert to absolute paths
INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"
OUTPUT_DIR="$(mkdir -p "$OUTPUT_DIR" && cd "$OUTPUT_DIR" && pwd)"
MANIFEST_FILE="$(cd "$(dirname "$MANIFEST_FILE")" && pwd)/$(basename "$MANIFEST_FILE")"

IPSAE_OUTPUT="${OUTPUT_DIR}/ipsae_scores"
SCORE_OUTPUT="${OUTPUT_DIR}/scores"

# Calculate chunks
TOTAL_JOBS=$(($(wc -l < "$MANIFEST_FILE") - 1))
TOTAL_CHUNKS=$(( (TOTAL_JOBS + MAX_ARRAY_SIZE - 1) / MAX_ARRAY_SIZE ))
[[ -z "$END_CHUNK" ]] && END_CHUNK=$TOTAL_CHUNKS

# ─────────────────────────────────────────────────────────────────────────────
# Display info
# ─────────────────────────────────────────────────────────────────────────────
echo "==============================================================================="
echo "AF3 Homomeric Controller Job Started"
echo "==============================================================================="
echo ""
echo "Controller Job ID: $SLURM_JOB_ID"
echo "Started at:        $(date)"
echo ""
echo "Manifest:          $MANIFEST_FILE"
echo "Total jobs:        $TOTAL_JOBS"
echo "Total chunks:      $TOTAL_CHUNKS"
echo "Processing:        chunks $START_CHUNK to $END_CHUNK"
echo ""
echo "Configuration:"
echo "  Batch:           $BATCH_NAME"
echo "  Input:           $INPUT_DIR"
echo "  Output:          $OUTPUT_DIR"
echo "  Num protomers:   $NUM_PROTOMERS"
if [[ -n "$LIGAND" ]]; then
    echo "  Ligand:          $LIGAND (x${NUM_LIGANDS:-1})"
fi
echo "  Seeds:           $SEEDS"
echo "  Run ipSAE:       $RUN_IPSAE"
echo "  Run scores:      $RUN_SCORES"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Build base sbatch arguments
# ─────────────────────────────────────────────────────────────────────────────
SBATCH_BASE_ARGS=(
    --job-name="${BATCH_NAME}"
    "$JOB_SCRIPT"
    --input "$INPUT_DIR"
    --output "$OUTPUT_DIR"
    --num-protomers "$NUM_PROTOMERS"
    --generate-script "$GENERATE_SCRIPT"
    --run-ipsae "$RUN_IPSAE"
    --run-scores "$RUN_SCORES"
)

# Add ligand if specified
if [[ -n "$LIGAND" ]]; then
    SBATCH_BASE_ARGS+=(--ligand "$LIGAND")
    if [[ -n "$NUM_LIGANDS" ]]; then
        SBATCH_BASE_ARGS+=(--num-ligands "$NUM_LIGANDS")
    fi
fi

# Add seeds (handle multiple seeds properly)
SBATCH_BASE_ARGS+=(--seeds)
NUM_SEEDS=0
for seed in $SEEDS; do
    SBATCH_BASE_ARGS+=("$seed")
    NUM_SEEDS=$((NUM_SEEDS + 1))
done

# Add num-seeds if more than one seed (for seed expansion mode)
if [[ $NUM_SEEDS -gt 1 ]]; then
    SBATCH_BASE_ARGS+=(--num-seeds "$NUM_SEEDS")
fi

# Add optional scripts
if [[ "$RUN_IPSAE" == "true" && -f "$IPSAE_SCRIPT" ]]; then
    SBATCH_BASE_ARGS+=(--ipsae-script "$IPSAE_SCRIPT" --ipsae-output "$IPSAE_OUTPUT")
fi

if [[ "$RUN_SCORES" == "true" && -f "$SCORE_SCRIPT" ]]; then
    SBATCH_BASE_ARGS+=(--score-script "$SCORE_SCRIPT" --score-output "$SCORE_OUTPUT")
fi

# ─────────────────────────────────────────────────────────────────────────────
# Helper: wait for a job array to complete
# ─────────────────────────────────────────────────────────────────────────────
wait_for_job() {
    local job_id=$1
    local chunk_num=$2
    local chunk_size=$3
    
    echo "  Waiting for job $job_id to complete..."
    
    while true; do
        # Check if any tasks are still running or pending
        local remaining=$(squeue -j "$job_id" -h 2>/dev/null | wc -l || echo "0")
        
        if [[ $remaining -eq 0 ]]; then
            # Job finished - check how many succeeded
            local completed=$(sacct -j "$job_id" --format=State -n 2>/dev/null | grep -c "COMPLETED" || echo "0")
            local failed=$(sacct -j "$job_id" --format=State -n 2>/dev/null | grep -c "FAILED" || echo "0")
            echo "  Chunk $chunk_num finished: $completed completed, $failed failed"
            return 0
        fi
        
        # Show progress
        local done_count=$((chunk_size - remaining))
        printf "\r  Chunk %d progress: %d/%d tasks running/pending...   " "$chunk_num" "$remaining" "$chunk_size"
        
        sleep $CHECK_INTERVAL
    done
}

# ─────────────────────────────────────────────────────────────────────────────
# Process chunks one at a time
# ─────────────────────────────────────────────────────────────────────────────
echo "─────────────────────────────────────────────────────────────────────────────"
echo "Processing chunks sequentially (one at a time)..."
echo "─────────────────────────────────────────────────────────────────────────────"
echo ""

PROGRESS_FILE="${OUTPUT_DIR}/.controller_progress"
START_TIME=$(date +%s)

for CHUNK_NUM in $(seq $START_CHUNK $END_CHUNK); do
    # Calculate job range
    START_IDX=$(( (CHUNK_NUM - 1) * MAX_ARRAY_SIZE ))
    END_IDX=$(( START_IDX + MAX_ARRAY_SIZE - 1 ))
    [[ $END_IDX -ge $TOTAL_JOBS ]] && END_IDX=$((TOTAL_JOBS - 1))
    
    CHUNK_SIZE=$((END_IDX - START_IDX + 1))
    ARRAY_RANGE="0-$((CHUNK_SIZE - 1))"
    
    # Calculate elapsed time
    ELAPSED=$(($(date +%s) - START_TIME))
    ELAPSED_H=$((ELAPSED / 3600))
    ELAPSED_M=$(( (ELAPSED % 3600) / 60 ))
    
    echo ""
    echo "[${ELAPSED_H}h ${ELAPSED_M}m] ═══════════════════════════════════════════════════"
    echo "Chunk $CHUNK_NUM/$TOTAL_CHUNKS: jobs $START_IDX-$END_IDX ($CHUNK_SIZE jobs)"
    echo "═══════════════════════════════════════════════════════════════════════════"
    
    # Build command
    SBATCH_CMD=(sbatch --array="$ARRAY_RANGE" "${SBATCH_BASE_ARGS[@]}")
    [[ $START_IDX -gt 0 ]] && SBATCH_CMD+=(--array-offset "$START_IDX")
    
    # Submit
    echo "  Submitting..."
    JOB_OUTPUT=$("${SBATCH_CMD[@]}" 2>&1)
    SBATCH_EXIT=$?
    
    if [[ $SBATCH_EXIT -ne 0 ]]; then
        echo "  ERROR: Submission failed: $JOB_OUTPUT"
        echo "  Retrying in 60 seconds..."
        sleep 60
        JOB_OUTPUT=$("${SBATCH_CMD[@]}" 2>&1)
        SBATCH_EXIT=$?
        if [[ $SBATCH_EXIT -ne 0 ]]; then
            echo "  FATAL: Submission failed again. Exiting."
            echo "  To resume: --start-chunk $CHUNK_NUM"
            exit 1
        fi
    fi
    
    JOB_ID=$(echo "$JOB_OUTPUT" | grep -oE '[0-9]+' | tail -1)
    echo "  Submitted: Job $JOB_ID"
    
    # Save progress
    echo "$CHUNK_NUM $JOB_ID $(date +%s)" >> "$PROGRESS_FILE"
    
    # Wait for this chunk to complete
    wait_for_job "$JOB_ID" "$CHUNK_NUM" "$CHUNK_SIZE"
    
    # Brief pause before next submission
    sleep 5
done

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
TOTAL_TIME=$(($(date +%s) - START_TIME))
TOTAL_H=$((TOTAL_TIME / 3600))
TOTAL_M=$(( (TOTAL_TIME % 3600) / 60 ))

echo ""
echo "==============================================================================="
echo "Controller Job Complete!"
echo "==============================================================================="
echo ""
echo "Processed:    $((END_CHUNK - START_CHUNK + 1)) chunks"
echo "Total time:   ${TOTAL_H}h ${TOTAL_M}m"
echo "Finished at:  $(date)"
echo ""
echo "Check outputs:"
echo "  ls $OUTPUT_DIR/*.tar.gz | wc -l"
if [[ "$RUN_SCORES" == "true" ]]; then
    echo "  ls $SCORE_OUTPUT/*.tsv | wc -l"
fi
echo ""
