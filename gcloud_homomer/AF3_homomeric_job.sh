#!/bin/bash
#===============================================================================
# AF3_homomeric_job.sh
#
# SLURM array job script for AlphaFold 3 homomeric predictions.
# Generates input JSON on-the-fly, runs AlphaFold, then cleans up.
#
# This script is called by AF3_homomeric_submit.sh - do not run directly.
#
# The script computes which input file to process based on SLURM_ARRAY_TASK_ID:
#   file_index = SLURM_ARRAY_TASK_ID / NUM_SEEDS (if seed expansion)
#   file_index = SLURM_ARRAY_TASK_ID (otherwise)
#===============================================================================

#SBATCH --job-name=alphafold3_homomeric
#SBATCH --partition=infa2u
#SBATCH --time=01:00:00
#SBATCH --mem=160G
#SBATCH --cpus-per-task=12
#SBATCH --output=/home/kamounlab_gmail_com/slurm_logs/job_%j/out.txt
#SBATCH --error=/home/kamounlab_gmail_com/slurm_logs/job_%j/err.txt
#SBATCH --gres=gpu:1

# Print all arguments for debugging
echo "==============================================================================="
echo "AF3 Homomeric Job Script Starting"
echo "==============================================================================="
echo "Date: $(date)"
echo "Hostname: $(hostname)"
echo "Working directory: $(pwd)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-not set}"
echo "SLURM_ARRAY_JOB_ID: ${SLURM_ARRAY_JOB_ID:-not set}"
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID:-not set}"
echo ""
echo "Arguments received ($# total):"
for arg in "$@"; do
    echo "  $arg"
done
echo "==============================================================================="
echo ""

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Parse arguments
# ─────────────────────────────────────────────────────────────────────────────
INPUT_DIR=""
OUTPUT_DIR=""
NUM_PROTOMERS=""
LIGAND=""
NUM_LIGANDS=""
NUM_SEEDS=""
GENERATE_SCRIPT=""
SEEDS=""
IPSAE_SCRIPT=""
IPSAE_OUTPUT_DIR=""
SCORE_SCRIPT=""
SCORE_OUTPUT_DIR=""
RUN_IPSAE="true"
RUN_SCORES="true"
ARRAY_OFFSET=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --input)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output)
            OUTPUT_DIR="$2"
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
        --num-seeds)
            NUM_SEEDS="$2"
            shift 2
            ;;
        --generate-script)
            GENERATE_SCRIPT="$2"
            shift 2
            ;;
        --seeds)
            shift
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                SEEDS="$SEEDS $1"
                shift
            done
            SEEDS="${SEEDS# }"  # Trim leading space
            ;;
        --ipsae-script)
            IPSAE_SCRIPT="$2"
            shift 2
            ;;
        --ipsae-output)
            IPSAE_OUTPUT_DIR="$2"
            shift 2
            ;;
        --score-script)
            SCORE_SCRIPT="$2"
            shift 2
            ;;
        --score-output)
            SCORE_OUTPUT_DIR="$2"
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
        --array-offset)
            ARRAY_OFFSET="$2"
            shift 2
            ;;
        *)
            echo "ERROR: Unknown argument: $1"
            exit 1
            ;;
    esac
done

# ─────────────────────────────────────────────────────────────────────────────
# Validate required arguments
# ─────────────────────────────────────────────────────────────────────────────
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" || -z "$NUM_PROTOMERS" ]]; then
    echo "ERROR: Missing required arguments"
    echo "Required: --input, --output, --num-protomers"
    exit 1
fi

if [[ -z "$GENERATE_SCRIPT" ]]; then
    echo "ERROR: --generate-script is required"
    exit 1
fi

if [[ ! -f "$GENERATE_SCRIPT" ]]; then
    echo "ERROR: Generate script not found: $GENERATE_SCRIPT"
    exit 1
fi

# ─────────────────────────────────────────────────────────────────────────────
# Load environment
# ─────────────────────────────────────────────────────────────────────────────
# NOTE: Add your environment setup here if needed
# source package e8edb411-7374-4342-b9f1-408da41fc197

# ─────────────────────────────────────────────────────────────────────────────
# Calculate indices from SLURM_ARRAY_TASK_ID
# ─────────────────────────────────────────────────────────────────────────────
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
JOB_ID=${SLURM_ARRAY_JOB_ID:-$$}

# Apply array offset if provided (for chunked submissions beyond MaxArraySize)
if [[ -n "$ARRAY_OFFSET" ]]; then
    TASK_ID=$((TASK_ID + ARRAY_OFFSET))
fi

# Handle seed expansion mode
if [[ -n "$NUM_SEEDS" && "$NUM_SEEDS" -gt 1 ]]; then
    # Seed expansion: task_id encodes (file_index, seed_index)
    FILE_INDEX=$((TASK_ID / NUM_SEEDS))
    SEED_INDEX=$((TASK_ID % NUM_SEEDS))
    
    # Extract the specific seed from SEEDS string
    SEEDS_ARRAY=($SEEDS)
    CURRENT_SEED="${SEEDS_ARRAY[$SEED_INDEX]}"
else
    # Standard mode: task_id = file_index
    FILE_INDEX=$TASK_ID
    CURRENT_SEED=""
fi

echo "==============================================================================="
echo "AF3 Homomeric Pipeline Job"
echo "==============================================================================="
echo "Job ID:           ${JOB_ID}_${TASK_ID}"
echo "Task ID:          $TASK_ID"
echo "File Index:       $FILE_INDEX"
echo "Input dir:        $INPUT_DIR"
echo "Output dir:       $OUTPUT_DIR"
echo "Num protomers:    $NUM_PROTOMERS"
if [[ -n "$LIGAND" ]]; then
    echo "Ligand:           $LIGAND (x${NUM_LIGANDS:-1})"
fi
if [[ -n "$CURRENT_SEED" ]]; then
    echo "Seed (expanded):  $CURRENT_SEED"
elif [[ -n "$SEEDS" ]]; then
    echo "Seeds:            $SEEDS"
fi
echo "Run ipSAE:        $RUN_IPSAE"
echo "Run scores:       $RUN_SCORES"
echo "==============================================================================="
echo ""

# Set up log directory for tracking failures
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "$LOG_DIR"
FAILED_LOG="${LOG_DIR}/failed_jobs.tsv"

# Create header for failed log if it doesn't exist
if [[ ! -f "$FAILED_LOG" ]]; then
    echo -e "job_id\ttask_id\tfile_index\tprotein\toutput_name\texit_code\terror_message\ttimestamp" > "$FAILED_LOG"
fi

# ─────────────────────────────────────────────────────────────────────────────
# Create temporary JSON file
# ─────────────────────────────────────────────────────────────────────────────
TEMP_JSON="/tmp/af3_homo_${JOB_ID}_${TASK_ID}.json"

echo "Generating input JSON..."
echo "  Temp file: $TEMP_JSON"

# Build generate command
GENERATE_CMD=(
    python3 "$GENERATE_SCRIPT"
    "$INPUT_DIR"
    --index "$FILE_INDEX"
    --protomer-count "$NUM_PROTOMERS"
    --output-file "$TEMP_JSON"
    --job-number "$TASK_ID"
    --quiet
)

# Add ligand if specified
if [[ -n "$LIGAND" ]]; then
    GENERATE_CMD+=(--ligand "$LIGAND")
    if [[ -n "$NUM_LIGANDS" ]]; then
        GENERATE_CMD+=(--ligand-count "$NUM_LIGANDS")
    fi
fi

# Add seed(s) - use single seed if expanded, all seeds otherwise
if [[ -n "$CURRENT_SEED" ]]; then
    GENERATE_CMD+=(--seeds "$CURRENT_SEED")
elif [[ -n "$SEEDS" ]]; then
    GENERATE_CMD+=(--seeds $SEEDS)
fi

# Run generate and capture output
GENERATE_OUTPUT=$("${GENERATE_CMD[@]}")
echo "  Generate output: $GENERATE_OUTPUT"

# Parse the JSON output to get names
PROTEIN_NAME=$(echo "$GENERATE_OUTPUT" | python3 -c "import sys,json; print(json.load(sys.stdin)['protein'])")
OUTPUT_NAME=$(echo "$GENERATE_OUTPUT" | python3 -c "import sys,json; print(json.load(sys.stdin)['output_name'])")

echo ""
echo "Processing: $PROTEIN_NAME (${NUM_PROTOMERS} protomers)"
echo "Output name: $OUTPUT_NAME"

# Verify temp file was created
if [[ ! -f "$TEMP_JSON" ]]; then
    echo "ERROR: Failed to create input JSON: $TEMP_JSON"
    exit 1
fi

# ─────────────────────────────────────────────────────────────────────────────
# Set up local scratch directory (fast SSD storage on compute node)
# ─────────────────────────────────────────────────────────────────────────────
LOCAL_SCRATCH="/tmp/af3_scratch_${JOB_ID}_${TASK_ID}"
mkdir -p "$LOCAL_SCRATCH"

echo "Local scratch:    $LOCAL_SCRATCH"
echo "Final output:     $OUTPUT_DIR/$OUTPUT_NAME.tar.gz"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Run AlphaFold 3
# ─────────────────────────────────────────────────────────────────────────────
echo "─────────────────────────────────────────────────────────────────────────────"
echo "Running AlphaFold 3..."
echo "─────────────────────────────────────────────────────────────────────────────"

apptainer run --nv \
    --bind "$INPUT_DIR" \
    --bind "$OUTPUT_DIR" \
    --bind /opt/apps/af3/models \
    --bind "$LOCAL_SCRATCH" \
    /opt/apps/af3/containers/af3.sif python3 /app/alphafold/run_alphafold.py \
    --json_path="$TEMP_JSON" \
    --model_dir=/opt/apps/af3/models \
    --db_dir=/dev/shm/public_databases \
    --pdb_database_path=/mnt/databases/v3.0/uncompressed/mmcif_files \
    --run_data_pipeline=false \
    --run_inference=True \
    --output_dir="$LOCAL_SCRATCH"

AF_EXIT_CODE=$?

# Update AF3_OUTPUT to point to local scratch location
AF3_LOCAL_OUTPUT="${LOCAL_SCRATCH}/${OUTPUT_NAME}"

# ─────────────────────────────────────────────────────────────────────────────
# Extract confidence scores (if enabled and AF3 succeeded)
# ─────────────────────────────────────────────────────────────────────────────
SCORES_EXIT_CODE=0
if [[ $AF_EXIT_CODE -eq 0 && "$RUN_SCORES" == "true" && -n "$SCORE_SCRIPT" && -f "$SCORE_SCRIPT" ]]; then
    echo ""
    echo "─────────────────────────────────────────────────────────────────────────────"
    echo "Extracting confidence scores..."
    echo "─────────────────────────────────────────────────────────────────────────────"
    
    # Find the summary_confidences.json file
    SUMMARY_JSON=$(find "$AF3_LOCAL_OUTPUT" -maxdepth 1 -name "*_summary_confidences.json" -type f | head -1)
    
    if [[ -z "$SUMMARY_JSON" ]]; then
        echo "WARNING: Could not find summary_confidences.json in $AF3_LOCAL_OUTPUT"
        echo "  Skipping score extraction"
    else
        echo "  Summary file: $(basename "$SUMMARY_JSON")"
        
        # Determine output directory for scores
        if [[ -n "$SCORE_OUTPUT_DIR" ]]; then
            SCORE_DEST="$SCORE_OUTPUT_DIR"
        else
            SCORE_DEST="${OUTPUT_DIR}/scores"
        fi
        mkdir -p "$SCORE_DEST"
        
        # Output base path (script creates _scores.tsv and _matrix.tsv)
        SCORE_OUTPUT_BASE="${SCORE_DEST}/${OUTPUT_NAME}"
        
        set +e
        python3 "$SCORE_SCRIPT" "$SUMMARY_JSON" --output "$SCORE_OUTPUT_BASE" --model-name "$OUTPUT_NAME" --quiet
        SCORES_EXIT_CODE=$?
        set -e
        
        if [[ $SCORES_EXIT_CODE -eq 0 ]]; then
            echo "  Overview:  ${OUTPUT_NAME}_scores.tsv"
            echo "  Matrix:    ${OUTPUT_NAME}_matrix.tsv"
        else
            echo "  WARNING: Score extraction failed with exit code: $SCORES_EXIT_CODE"
        fi
    fi
fi

# ─────────────────────────────────────────────────────────────────────────────
# Run ipSAE scoring (if enabled and AF3 succeeded)
# ─────────────────────────────────────────────────────────────────────────────
IPSAE_EXIT_CODE=0
if [[ $AF_EXIT_CODE -eq 0 && "$RUN_IPSAE" == "true" && -n "$IPSAE_SCRIPT" ]]; then
    echo ""
    echo "─────────────────────────────────────────────────────────────────────────────"
    echo "Running ipSAE scoring on top model..."
    echo "─────────────────────────────────────────────────────────────────────────────"
    
    # Find the *_model.cif file (top model, not in seed subdirectories)
    TOP_MODEL_CIF=$(find "$AF3_LOCAL_OUTPUT" -maxdepth 1 -name "*_model.cif" -type f | head -1)
    
    # Find the *_confidences.json file (not summary_confidences.json)
    TOP_CONFIDENCES=$(find "$AF3_LOCAL_OUTPUT" -maxdepth 1 -name "*_confidences.json" ! -name "*summary*" -type f | head -1)
    
    if [[ -z "$TOP_MODEL_CIF" ]]; then
        echo "WARNING: Could not find top model .cif file in $AF3_LOCAL_OUTPUT"
        echo "  Skipping ipSAE scoring"
        # Log to failed TSV
        IPSAE_FAILED_FILE="${OUTPUT_DIR}/ipsae_failed.tsv"
        if [[ ! -f "$IPSAE_FAILED_FILE" ]]; then
            echo -e "job_num\toutput_name\tprotein\texit_code" > "$IPSAE_FAILED_FILE"
        fi
        (
            flock -x 200
            echo -e "${TASK_ID}\t${OUTPUT_NAME}\t${PROTEIN_NAME}\tno_cif" >> "$IPSAE_FAILED_FILE"
        ) 200>"${IPSAE_FAILED_FILE}.lock"
    elif [[ -z "$TOP_CONFIDENCES" ]]; then
        echo "WARNING: Could not find confidences.json file in $AF3_LOCAL_OUTPUT"
        echo "  Skipping ipSAE scoring"
        # Log to failed TSV
        IPSAE_FAILED_FILE="${OUTPUT_DIR}/ipsae_failed.tsv"
        if [[ ! -f "$IPSAE_FAILED_FILE" ]]; then
            echo -e "job_num\toutput_name\tprotein\texit_code" > "$IPSAE_FAILED_FILE"
        fi
        (
            flock -x 200
            echo -e "${TASK_ID}\t${OUTPUT_NAME}\t${PROTEIN_NAME}\tno_json" >> "$IPSAE_FAILED_FILE"
        ) 200>"${IPSAE_FAILED_FILE}.lock"
    elif [[ ! -f "$IPSAE_SCRIPT" ]]; then
        echo "WARNING: ipSAE script not found: $IPSAE_SCRIPT"
        echo "  Skipping ipSAE scoring"
    else
        echo "  Model file: $(basename "$TOP_MODEL_CIF")"
        echo "  Confidences: $(basename "$TOP_CONFIDENCES")"
        
        # Determine output directory for ipsae results
        if [[ -n "$IPSAE_OUTPUT_DIR" ]]; then
            IPSAE_DEST="$IPSAE_OUTPUT_DIR"
        else
            IPSAE_DEST="${OUTPUT_DIR}/ipsae_scores"
        fi
        mkdir -p "$IPSAE_DEST"
        
        # Get the directory containing the ipsae script for binding
        IPSAE_SCRIPT_DIR=$(dirname "$IPSAE_SCRIPT")
        IPSAE_SCRIPT_NAME=$(basename "$IPSAE_SCRIPT")
        
        # Run ipsae.py inside the apptainer container (has numpy)
        echo "  Running ipSAE inside container..."
        
        # Disable exit-on-error for ipSAE so we can capture failures
        set +e
        apptainer exec \
            --bind "$IPSAE_SCRIPT_DIR" \
            --bind "$AF3_LOCAL_OUTPUT" \
            --bind "$IPSAE_DEST" \
            /opt/apps/af3/containers/af3.sif \
            python3 "${IPSAE_SCRIPT_DIR}/${IPSAE_SCRIPT_NAME}" \
            "$TOP_CONFIDENCES" "$TOP_MODEL_CIF" 10 10
        
        IPSAE_EXIT_CODE=$?
        set -e
        
        if [[ $IPSAE_EXIT_CODE -eq 0 ]]; then
            # Move ipsae output files to destination
            MODEL_STEM=$(basename "$TOP_MODEL_CIF" .cif)
            
            for ext in ".txt" "_byres.txt" ".pml"; do
                IPSAE_FILE="${AF3_LOCAL_OUTPUT}/${MODEL_STEM}_10_10${ext}"
                if [[ -f "$IPSAE_FILE" ]]; then
                    mv "$IPSAE_FILE" "$IPSAE_DEST/"
                    echo "  Moved: $(basename "$IPSAE_FILE") -> $IPSAE_DEST/"
                fi
            done
            echo "  ipSAE scoring completed successfully"
        else
            echo "  WARNING: ipSAE scoring failed with exit code: $IPSAE_EXIT_CODE"
            
            # Log failed ipSAE to TSV file
            IPSAE_FAILED_FILE="${OUTPUT_DIR}/ipsae_failed.tsv"
            
            if [[ ! -f "$IPSAE_FAILED_FILE" ]]; then
                echo -e "job_num\toutput_name\tprotein\texit_code" > "$IPSAE_FAILED_FILE"
            fi
            
            (
                flock -x 200
                echo -e "${TASK_ID}\t${OUTPUT_NAME}\t${PROTEIN_NAME}\t${IPSAE_EXIT_CODE}" >> "$IPSAE_FAILED_FILE"
            ) 200>"${IPSAE_FAILED_FILE}.lock"
            
            echo "  Logged to: $IPSAE_FAILED_FILE"
        fi
    fi
fi

# ─────────────────────────────────────────────────────────────────────────────
# Compress and transfer AF3 output to final storage
# ─────────────────────────────────────────────────────────────────────────────
if [[ $AF_EXIT_CODE -eq 0 && -d "$AF3_LOCAL_OUTPUT" ]]; then
    echo ""
    echo "─────────────────────────────────────────────────────────────────────────────"
    echo "Compressing and transferring output to final storage..."
    echo "─────────────────────────────────────────────────────────────────────────────"
    
    # Ensure final output directory exists
    mkdir -p "$OUTPUT_DIR"
    
    # Compress the AF3 output directory
    ARCHIVE_NAME="${OUTPUT_NAME}.tar.gz"
    ARCHIVE_PATH="${LOCAL_SCRATCH}/${ARCHIVE_NAME}"
    
    echo "  Compressing: $AF3_LOCAL_OUTPUT"
    echo "  Archive:     $ARCHIVE_NAME"
    
    # Create tar.gz archive (from parent directory to preserve folder name)
    cd "$LOCAL_SCRATCH"
    tar -czf "$ARCHIVE_NAME" "$OUTPUT_NAME"
    COMPRESS_EXIT_CODE=$?
    
    if [[ $COMPRESS_EXIT_CODE -eq 0 ]]; then
        # Get archive size for logging
        ARCHIVE_SIZE=$(du -h "$ARCHIVE_PATH" | cut -f1)
        echo "  Archive size: $ARCHIVE_SIZE"
        
        # Move archive to final storage
        echo "  Transferring to: $OUTPUT_DIR/"
        mv "$ARCHIVE_PATH" "$OUTPUT_DIR/"
        TRANSFER_EXIT_CODE=$?
        
        if [[ $TRANSFER_EXIT_CODE -eq 0 ]]; then
            echo "  Transfer completed successfully"
        else
            echo "  WARNING: Transfer failed with exit code: $TRANSFER_EXIT_CODE"
        fi
    else
        echo "  WARNING: Compression failed with exit code: $COMPRESS_EXIT_CODE"
    fi
fi

# ─────────────────────────────────────────────────────────────────────────────
# Cleanup
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "─────────────────────────────────────────────────────────────────────────────"
echo "Cleanup"
echo "─────────────────────────────────────────────────────────────────────────────"

# Remove temp JSON file
if [[ -f "$TEMP_JSON" ]]; then
    rm -f "$TEMP_JSON"
    echo "Removed temp JSON: $TEMP_JSON"
fi

# Remove local scratch directory
if [[ -d "$LOCAL_SCRATCH" ]]; then
    rm -rf "$LOCAL_SCRATCH"
    echo "Removed local scratch: $LOCAL_SCRATCH"
fi

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
echo ""
echo "==============================================================================="
if [[ $AF_EXIT_CODE -eq 0 ]]; then
    echo "Job completed successfully!"
    echo "  Archive: $OUTPUT_DIR/${OUTPUT_NAME}.tar.gz"
    if [[ "$RUN_SCORES" == "true" && -n "$SCORE_SCRIPT" ]]; then
        echo "  Scores:  ${SCORE_DEST}/${OUTPUT_NAME}_scores.tsv"
        echo "  Matrix:  ${SCORE_DEST}/${OUTPUT_NAME}_matrix.tsv"
    fi
    if [[ "$RUN_IPSAE" == "true" && -n "$IPSAE_SCRIPT" ]]; then
        echo "  ipSAE:   $IPSAE_DEST/"
    fi
else
    echo "Job failed with exit code: $AF_EXIT_CODE"
    
    # Log the failure
    ERROR_MSG="AlphaFold3 inference failed"
    TIMESTAMP=$(date -Iseconds)
    
    # Append to failed jobs log (with file locking to handle concurrent writes)
    (
        flock -x 200
        echo -e "${JOB_ID}\t${TASK_ID}\t${FILE_INDEX}\t${PROTEIN_NAME}\t${OUTPUT_NAME}\t${AF_EXIT_CODE}\t${ERROR_MSG}\t${TIMESTAMP}" >> "$FAILED_LOG"
    ) 200>"${FAILED_LOG}.lock"
    
    echo "  Logged to: $FAILED_LOG"
fi
echo "==============================================================================="

exit $AF_EXIT_CODE
