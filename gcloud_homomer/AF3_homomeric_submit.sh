#!/bin/bash
#===============================================================================
# AF3_homomeric_submit.sh
# 
# Submit AlphaFold 3 homomeric modeling array jobs.
#
# This pipeline is for single protein modeling with multiple copies (homomers)
# and optional ligands.
#
# Usage:
#   AF3_homomeric_submit.sh <input_dir> <output_dir> [options]
#
# Options:
#   --batch NAME        Batch name for manifest (default: timestamp)
#   --protomers N       Number of protomers (required)
#   --ligand SPEC       Ligand: CCD code (CCD_OLA, ATP) or SMILES (smiles:CCO)
#   --num-ligands N     Number of ligand copies (default: 1 if ligand specified)
#   --seeds N [M..]     Custom seed numbers
#   --expand-seeds      Create separate jobs for each seed
#   --run-ipsae BOOL    Run ipSAE scoring (default: true)
#   --run-scores BOOL   Extract confidence scores (default: true)
#   --dry-run           Preview combinations without submitting
#   --use-controller    Use controller for large runs (>1000 jobs)
#
# Naming convention: {input}_{X}protomers_{Y}{Ligand}Ligs_seed{N}
#
# Example:
#   AF3_homomeric_submit.sh ./nlrs ./results --protomers 6 --batch NLR_homo
#   AF3_homomeric_submit.sh ./nlrs ./results --protomers 6 --ligand CCD_OLA --num-ligands 50
#===============================================================================

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Default configuration
# ─────────────────────────────────────────────────────────────────────────────
BATCH_NAME=""
NUM_PROTOMERS=""
LIGAND=""
NUM_LIGANDS=""
SEEDS=""
EXPAND_SEEDS=false
DRY_RUN=false
USE_CONTROLLER=false
RUN_IPSAE="true"
RUN_SCORES="true"
MANIFEST_PATH=""

# Detect script directory for default paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ─────────────────────────────────────────────────────────────────────────────
# Usage function
# ─────────────────────────────────────────────────────────────────────────────
usage() {
    cat << EOF
Usage: $(basename "$0") <input_dir> <output_dir> [options]

Submit AlphaFold 3 homomeric modeling jobs.

Required Arguments:
  input_dir       Directory containing input JSON files (with MSA/templates)
  output_dir      Directory for AlphaFold outputs

Required Options:
  --protomers N       Number of protomers

Optional:
  --batch NAME        Batch name for manifest and job naming
                      (default: af3_homo_YYYYMMDD_HHMMSS)
  --ligand SPEC       Ligand specification:
                        CCD code: CCD_OLA, ATP, OLA
                        SMILES: smiles:CCO
  --num-ligands N     Number of ligand copies (default: 1 if ligand specified)
  --seeds N [M..]     Custom seed numbers (space-separated)
  --expand-seeds      Create separate jobs for each seed value
  --run-ipsae BOOL    Run ipSAE scoring (true/false, default: true)
  --run-scores BOOL   Extract scores from summary_confidences.json
                      (true/false, default: true)
  --dry-run           Preview manifest and job info without submitting
  --use-controller    Use controller job for runs >1000 jobs
  --manifest FILE     Custom manifest output path
  --script-dir PATH   Path to pipeline scripts directory
  -h, --help          Show this help message

Output Organization:
  output_dir/
    ├── *.tar.gz               AF3 output archives
    ├── scores/                Confidence scores (TSV)
    │   ├── *_scores.tsv       Overview metrics
    │   └── *_matrix.tsv       Per-chain pairwise matrices
    ├── ipsae_scores/          ipSAE scoring results
    │   ├── *.txt
    │   └── *_byres.txt
    ├── logs/                  Job tracking
    │   └── failed_jobs.tsv    Failed job log
    └── batch_manifest.tsv     Job manifest

Examples:
  # Basic homomeric modeling with 6 protomers
  $(basename "$0") ./nlrs ./results --protomers 6 --batch NLR_homo

  # Add OLA ligand (50 copies)
  $(basename "$0") ./nlrs ./results --protomers 6 --ligand CCD_OLA --num-ligands 50

  # With multiple seeds (expanded)
  $(basename "$0") ./nlrs ./results --protomers 6 --seeds 1 2 3 --expand-seeds

  # Skip ipSAE but keep score extraction
  $(basename "$0") ./nlrs ./results --protomers 6 --run-ipsae false

  # Dry run to preview
  $(basename "$0") ./nlrs ./results --protomers 6 --dry-run
EOF
    exit "${1:-0}"
}

# ─────────────────────────────────────────────────────────────────────────────
# Parse arguments
# ─────────────────────────────────────────────────────────────────────────────
POSITIONAL_ARGS=()
SEED_ARGS=()

while [[ $# -gt 0 ]]; do
    case $1 in
        --batch)
            BATCH_NAME="$2"
            shift 2
            ;;
        --protomers)
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
            shift
            while [[ $# -gt 0 && "$1" =~ ^[0-9]+$ ]]; do
                SEED_ARGS+=("$1")
                shift
            done
            ;;
        --expand-seeds)
            EXPAND_SEEDS=true
            shift
            ;;
        --run-ipsae)
            RUN_IPSAE="$2"
            shift 2
            ;;
        --run-scores)
            RUN_SCORES="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --use-controller)
            USE_CONTROLLER=true
            shift
            ;;
        --manifest)
            MANIFEST_PATH="$2"
            shift 2
            ;;
        --script-dir)
            SCRIPT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            usage 0
            ;;
        -*)
            echo "ERROR: Unknown option: $1"
            usage 1
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done

# Restore positional arguments
set -- "${POSITIONAL_ARGS[@]}"

# ─────────────────────────────────────────────────────────────────────────────
# Validate arguments
# ─────────────────────────────────────────────────────────────────────────────
if [[ $# -lt 2 ]]; then
    echo "ERROR: Missing required arguments"
    echo ""
    usage 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

if [[ -z "$NUM_PROTOMERS" ]]; then
    echo "ERROR: --protomers is required"
    usage 1
fi

# Set script paths
GENERATE_SCRIPT="${SCRIPT_DIR}/generate_af3_homomeric.py"
JOB_SCRIPT="${SCRIPT_DIR}/AF3_homomeric_job.sh"
CONTROLLER_SCRIPT="${SCRIPT_DIR}/AF3_homomeric_controller.sh"
IPSAE_SCRIPT="${SCRIPT_DIR}/ipsae.py"
SCORE_SCRIPT="${SCRIPT_DIR}/extract_af3_scores.py"

# Validate required scripts exist
if [[ ! -f "$GENERATE_SCRIPT" ]]; then
    echo "ERROR: Generate script not found: $GENERATE_SCRIPT"
    exit 1
fi

if [[ ! -f "$JOB_SCRIPT" ]]; then
    echo "ERROR: Job script not found: $JOB_SCRIPT"
    exit 1
fi

# Convert to absolute paths
INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)"

# Set default batch name
if [[ -z "$BATCH_NAME" ]]; then
    BATCH_NAME="af3_homo_$(date +%Y%m%d_%H%M%S)"
fi

# Set default manifest path
if [[ -z "$MANIFEST_PATH" ]]; then
    MANIFEST_PATH="${OUTPUT_DIR}/${BATCH_NAME}_manifest.tsv"
fi

# Build seeds string
if [[ ${#SEED_ARGS[@]} -gt 0 ]]; then
    SEEDS="${SEED_ARGS[*]}"
else
    SEEDS="1"
fi

# Set default num-ligands if ligand specified
if [[ -n "$LIGAND" && -z "$NUM_LIGANDS" ]]; then
    NUM_LIGANDS="1"
fi

# Verify input directory has content (Python script will handle finding files)
if [[ -z "$(ls -A "$INPUT_DIR")" ]]; then
    echo "ERROR: Input directory is empty: $INPUT_DIR"
    exit 1
fi

# ─────────────────────────────────────────────────────────────────────────────
# Generate manifest
# ─────────────────────────────────────────────────────────────────────────────
echo "==============================================================================="
echo "AF3 Homomeric Pipeline"
echo "==============================================================================="
echo ""
echo "Input directory:  $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Batch name:       $BATCH_NAME"
echo ""
echo "Configuration:"
echo "  Protomers:      $NUM_PROTOMERS"
if [[ -n "$LIGAND" ]]; then
    echo "  Ligand:         $LIGAND"
    echo "  Ligand count:   ${NUM_LIGANDS}"
fi
echo "  Seeds:          $SEEDS"
echo "  Expand seeds:   $EXPAND_SEEDS"
echo "  Run ipSAE:      $RUN_IPSAE"
echo "  Run scores:     $RUN_SCORES"
echo ""

# Build generate command for manifest
GENERATE_CMD=(
    python3 "$GENERATE_SCRIPT"
    "$INPUT_DIR"
    --list-combinations
    --batch "$BATCH_NAME"
    --protomers "$NUM_PROTOMERS"
    --seeds $SEEDS
)

if [[ -n "$LIGAND" ]]; then
    GENERATE_CMD+=(--ligand "$LIGAND" --num-ligands "$NUM_LIGANDS")
fi

if [[ "$EXPAND_SEEDS" == "true" ]]; then
    GENERATE_CMD+=(--expand-seeds)
fi

echo "Generating manifest..."
"${GENERATE_CMD[@]}" > "$MANIFEST_PATH"

# Count jobs
NUM_JOBS=$(($(wc -l < "$MANIFEST_PATH") - 1))  # Subtract header

echo "  Generated: $MANIFEST_PATH"
echo "  Total jobs: $NUM_JOBS"
echo ""

# ─────────────────────────────────────────────────────────────────────────────
# Dry run - just show preview
# ─────────────────────────────────────────────────────────────────────────────
if [[ "$DRY_RUN" == "true" ]]; then
    echo "─────────────────────────────────────────────────────────────────────────────"
    echo "DRY RUN - Manifest Preview (first 20 jobs):"
    echo "─────────────────────────────────────────────────────────────────────────────"
    head -21 "$MANIFEST_PATH" | column -t -s $'\t'
    
    if [[ $NUM_JOBS -gt 20 ]]; then
        echo "... ($((NUM_JOBS - 20)) more jobs)"
    fi
    
    echo ""
    echo "─────────────────────────────────────────────────────────────────────────────"
    echo "To submit, remove --dry-run flag"
    echo "─────────────────────────────────────────────────────────────────────────────"
    exit 0
fi

# ─────────────────────────────────────────────────────────────────────────────
# Check if controller is needed
# ─────────────────────────────────────────────────────────────────────────────
MAX_ARRAY_SIZE=1001

if [[ $NUM_JOBS -gt $MAX_ARRAY_SIZE || "$USE_CONTROLLER" == "true" ]]; then
    echo "─────────────────────────────────────────────────────────────────────────────"
    echo "Large run detected ($NUM_JOBS jobs) - using controller mode"
    echo "─────────────────────────────────────────────────────────────────────────────"
    echo ""
    
    if [[ ! -f "$CONTROLLER_SCRIPT" ]]; then
        echo "ERROR: Controller script not found: $CONTROLLER_SCRIPT"
        exit 1
    fi
    
    # Build controller command
    CONTROLLER_CMD=(
        sbatch "$CONTROLLER_SCRIPT"
        --manifest "$MANIFEST_PATH"
        --input "$INPUT_DIR"
        --output "$OUTPUT_DIR"
        --batch "$BATCH_NAME"
        --num-protomers "$NUM_PROTOMERS"
        --run-ipsae "$RUN_IPSAE"
        --run-scores "$RUN_SCORES"
        --script-dir "$SCRIPT_DIR"
        --seeds $SEEDS
    )
    
    if [[ -n "$LIGAND" ]]; then
        CONTROLLER_CMD+=(--ligand "$LIGAND" --num-ligands "$NUM_LIGANDS")
    fi
    
    echo "Submitting controller job..."
    "${CONTROLLER_CMD[@]}"
    
    echo ""
    echo "Controller job submitted. It will manage chunk submissions."
    exit 0
fi

# ─────────────────────────────────────────────────────────────────────────────
# Direct array job submission (for smaller runs)
# ─────────────────────────────────────────────────────────────────────────────
echo "─────────────────────────────────────────────────────────────────────────────"
echo "Submitting array job..."
echo "─────────────────────────────────────────────────────────────────────────────"
echo ""

# Create output subdirectories
mkdir -p "${OUTPUT_DIR}/scores"
mkdir -p "${OUTPUT_DIR}/ipsae_scores"

# Build sbatch command
ARRAY_RANGE="0-$((NUM_JOBS - 1))"

SBATCH_CMD=(
    sbatch
    --array="$ARRAY_RANGE"
    --job-name="$BATCH_NAME"
    "$JOB_SCRIPT"
    --input "$INPUT_DIR"
    --output "$OUTPUT_DIR"
    --num-protomers "$NUM_PROTOMERS"
    --generate-script "$GENERATE_SCRIPT"
    --run-ipsae "$RUN_IPSAE"
    --run-scores "$RUN_SCORES"
    --seeds $SEEDS
)

# Add ligand
if [[ -n "$LIGAND" ]]; then
    SBATCH_CMD+=(--ligand "$LIGAND" --num-ligands "$NUM_LIGANDS")
fi

# Add num-seeds if expanding
if [[ "$EXPAND_SEEDS" == "true" && ${#SEED_ARGS[@]} -gt 1 ]]; then
    SBATCH_CMD+=(--num-seeds "${#SEED_ARGS[@]}")
fi

# Add scoring scripts
if [[ "$RUN_IPSAE" == "true" && -f "$IPSAE_SCRIPT" ]]; then
    SBATCH_CMD+=(--ipsae-script "$IPSAE_SCRIPT" --ipsae-output "${OUTPUT_DIR}/ipsae_scores")
fi

if [[ "$RUN_SCORES" == "true" && -f "$SCORE_SCRIPT" ]]; then
    SBATCH_CMD+=(--score-script "$SCORE_SCRIPT" --score-output "${OUTPUT_DIR}/scores")
fi

echo "Command:"
echo "  ${SBATCH_CMD[*]}"
echo ""

# Submit
echo "Submitting..."
JOB_OUTPUT=$("${SBATCH_CMD[@]}" 2>&1)
SBATCH_EXIT=$?

echo "sbatch exit code: $SBATCH_EXIT"
echo "sbatch output: $JOB_OUTPUT"

if [[ $SBATCH_EXIT -ne 0 ]]; then
    echo "ERROR: Job submission failed:"
    echo "$JOB_OUTPUT"
    exit 1
fi

JOB_ID=$(echo "$JOB_OUTPUT" | grep -oE '[0-9]+' | tail -1)

echo "==============================================================================="
echo "Job Submitted Successfully!"
echo "==============================================================================="
echo ""
echo "Job ID:         $JOB_ID"
echo "Array range:    $ARRAY_RANGE"
echo "Total jobs:     $NUM_JOBS"
echo ""
echo "Manifest:       $MANIFEST_PATH"
echo "Output dir:     $OUTPUT_DIR"
if [[ "$RUN_SCORES" == "true" ]]; then
    echo "Scores dir:     ${OUTPUT_DIR}/scores"
fi
if [[ "$RUN_IPSAE" == "true" ]]; then
    echo "ipSAE dir:      ${OUTPUT_DIR}/ipsae_scores"
fi
echo ""
echo "Monitor with:"
echo "  squeue -j $JOB_ID"
echo "  sacct -j $JOB_ID"
echo ""
