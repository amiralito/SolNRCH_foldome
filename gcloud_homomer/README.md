# AF3 Homomeric Modeling Pipeline v2

Pipeline for AlphaFold 3 homomeric protein modeling with support for:
- Variable protomer counts
- Ligand support (CCD codes)
- Seed expansion for multiple seeds
- Score extraction from summary_confidences.json
- Optional ipSAE scoring
- Job tracking with failure logs

## Input Structure

The pipeline supports inputs in subdirectories (like your existing setup):
```
inference_input/
├── protein1/
│   └── protein1_data.json
├── protein2/
│   └── protein2_data.json
└── protein3/
│   └── protein3_data.json
```

Or direct JSON files in the input directory.

## Quick Start

```bash
# Copy scripts to your inference directory
cp *.sh *.py ~/inference_homomer_scripts/

# Basic usage: 6 protomers, 2 seeds expanded, with OLA ligand
bash ~/inference_homomer_scripts/AF3_homomeric_submit.sh ./inference_input ./output \
    --protomers 6 \
    --ligand CCD_OLA \
    --num-ligands 15 \
    --seeds 1 2 \
    --expand-seeds \
    --batch my_screen

# Dry run to preview
bash ~/inference_homomer_scripts/AF3_homomeric_submit.sh ./inference_input ./output \
    --protomers 6 --ligand CCD_OLA --num-ligands 15 \
    --seeds 1 2 --expand-seeds --dry-run
```

## Output Naming

Names now include job number prefix:
```
job_{N}_{protein}_{X}protomers_{Y}{Ligand}Ligs_seed{S}
```

Examples:
- `job_0_NLR1_6protomers_seed1`
- `job_1_NLR1_6protomers_seed2`
- `job_2_NLR1_6protomers_15OLALigs_seed1`

## Output Organization

```
output_dir/
├── *.tar.gz               AF3 output archives
├── scores/                Confidence scores
│   ├── *_scores.tsv       Overview metrics
│   └── *_matrix.tsv       Per-chain pairwise matrices
├── ipsae_scores/          ipSAE scoring results
├── logs/                  Job tracking
│   └── failed_jobs.tsv    Failed job log
└── batch_manifest.tsv     Job manifest
```

## Score Output Files

The pipeline outputs two score files per model:

### 1. Overview scores (`*_scores.tsv`)
Tabular format with one row per model:
```
model_name    seed    ptm    iptm    ranking_score    fraction_disordered    has_clash    mean_pair_iptm    min_pair_iptm    max_pair_iptm
```

### 2. Per-chain matrix (`*_matrix.tsv`)
Matrix format for pairwise scores:
```
# Model: job_0_protein_6protomers_seed1
# Seed: 1

## chain_ptm
chain   A       B       C       D       E       F
ptm     0.88    0.87    0.86    0.89    0.85    0.84

## chain_iptm
chain   A       B       C       D       E       F
iptm    0.75    0.73    0.71    0.74    0.72    0.70

## chain_pair_iptm
        A       B       C       D       E       F
A       1.00    0.82    0.78    0.75    0.73    0.71
B       0.82    1.00    0.80    0.77    0.74    0.72
...

## chain_pair_pae_min
        A       B       C       D       E       F
A       0.50    2.10    2.80    3.20    3.50    3.80
...
```

## Failed Job Tracking

Failed jobs are logged to `logs/failed_jobs.tsv`:
```
job_id  task_id  file_index  protein  output_name  exit_code  error_message  timestamp
```

This helps identify which models need to be rerun.

## Memory Considerations

AF3 memory usage scales quadratically with token count. For 80GB A100:
- Safe limit: ~3500-4000 tokens
- 6 protomers (~550 AA) + 15 OLA ligands ≈ 3600 tokens ✓
- 6 protomers + 50 OLA ligands ≈ 4200 tokens ✗ (OOM)

## Files

| File | Purpose |
|------|---------|
| `AF3_homomeric_submit.sh` | Main entry point |
| `AF3_homomeric_job.sh` | SLURM job script |
| `AF3_homomeric_controller.sh` | Controller for large runs (>1000 jobs) |
| `generate_af3_homomeric.py` | Generate AF3 input JSONs |
| `extract_af3_scores.py` | Extract scores from outputs |

## Options Reference

```
Required:
  input_dir           Directory with input JSON files (or subdirs)
  output_dir          Output directory
  --protomers N       Number of protomers

Optional:
  --batch NAME        Batch name (default: timestamp)
  --ligand SPEC       Ligand: CCD_OLA, OLA, ccd:OLA
  --num-ligands N     Number of ligand copies
  --seeds N [M..]     Seed numbers
  --expand-seeds      Create separate job per seed
  --run-ipsae BOOL    Run ipSAE (default: true)
  --run-scores BOOL   Extract scores (default: true)
  --dry-run           Preview without submitting
  --use-controller    Force controller mode
```
