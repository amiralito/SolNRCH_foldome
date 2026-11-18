# Resistosome Structure Analysis Pipeline (v2.0.0)

Batch analysis of predicted (AlphaFold 3) or empirical resistosome protein structures. Computes structural quality metrics for oligomeric NLR assemblies and applies statistical penalties (LCB/UCB) to handle variability across replicates and interfaces.

## Installation

```bash
# Python dependencies (pinned versions used for the paper)
pip install -r requirements.txt

# DSSP (system package)
sudo apt-get install dssp          # Ubuntu/Debian
# or: conda install -c salilab dssp

# HMMER3 (optional, for --hmm mode)
conda install -c bioconda hmmer    # or: brew install hmmer
```

## Quick Start

```bash
# Basic usage
python run_pipeline.py \
    --input_dir /path/to/af3_output \
    --output_dir results/ \
    --mhd_csv NRCH_MHD_nlrexp_coords.csv \
    --ploop_csv NRCH_ploop_nlrexp_coords.csv \
    --nbarc_csv NRCH_NBARC_coords.csv

# With HMM-based α1 detection (recommended for CNLs with MADA motif)
python run_pipeline.py \
    --input_dir /path/to/af3_output \
    --output_dir results/ \
    --hmm MADA.hmm \
    --mhd_csv NRCH_MHD_nlrexp_coords.csv \
    --ploop_csv NRCH_ploop_nlrexp_coords.csv \
    --nbarc_csv NRCH_NBARC_coords.csv

# With pre-computed hmmsearch results
python run_pipeline.py \
    --input_dir /path/to/af3_output \
    --output_dir results/ \
    --hmm_domtbl MADA_domtbl.out \
    --hmm_fallback_dssp \
    --protein_names protein_names.txt \
    --workers 12
```

## Command-Line Options

| Flag | Default | Description |
|------|---------|-------------|
| `--input_dir`, `-i` | *(required)* | Directory with structure files |
| `--output_dir`, `-o` | *(required)* | Output directory |
| `--mhd_csv` | None | MHD motif annotation CSV |
| `--ploop_csv` | None | P-loop motif annotation CSV |
| `--nbarc_csv` | None | NB-ARC domain boundary annotation CSV |
| `--protein_names` | None | Text file with known protein names (one per line) for flexible filename matching |
| `--hmm` | None | HMM model for α1-helix boundary detection (runs hmmsearch; requires HMMER3) |
| `--hmm_domtbl` | None | Pre-computed hmmsearch `--domtblout` file (alternative to `--hmm`) |
| `--hmm_fallback_dssp` | off | When HMM has no hit for a protein, fall back to DSSP-based α1 instead of reporting NA |
| `--workers` | auto (min(CPUs, 8)) | Number of parallel workers. Use `--workers 1` for serial/debug mode |
| `--contact_cutoff` | 4.0 | Heavy-atom contact cutoff (Å). This is also the search radius: pairs within it that match no specific contact type are counted as vdW |
| `--pae_cutoff` | 10.0 | PAE threshold for contact filtering (Å) |
| `--k` | 1.96 | Penalty weight k in LCB/UCB = mean ∓ k·SD/√n |
| `--work_dir` | `output_dir/tmp` | Temporary working directory |
| `--log_level` | INFO | Logging verbosity (DEBUG, INFO, WARNING, ERROR) |

## Input Formats

### Structure Files

The pipeline accepts three formats, automatically detected:

1. **AlphaFold 3 `.tar.gz` archives** — standard AF3 output. Extracts `*_model.cif`, `*_summary_confidences.json`, `*_confidences.json`. Subdirectories (`seed-N_sample-M/`) are ignored.

2. **Uncompressed AF3 output directories** — folders containing `*_model.cif` and confidence files directly. The pipeline scans up to 2 levels deep.

3. **Plain PDB/CIF files** — optionally `.gz` compressed. No confidence files required (ipTM and PAE metrics will be NA).

### File Naming Convention

Replicate/seed numbers are detected from `_seedN` in filenames. Structures with the same base ID but different seeds are grouped and aggregated:
```
protein_a_6protomers_25OLALigs_seed1.tar.gz
protein_a_6protomers_25OLALigs_seed2.tar.gz    → grouped as "protein_a"
protein_a_6protomers_25OLALigs_seed3.tar.gz
```

### Protein Names File (optional)

When AF3 output filenames are long or inconsistent, provide a `--protein_names` file for cleaner ID matching:
```
# protein_names.txt
nrc2a
nrc3
nrc4a
```
The pipeline matches filenames against these names (case-insensitive, longest match wins).

### Annotation CSVs (optional)

CSV files with columns: `ID_normalized`, `start`, `end`. IDs are matched against structure base IDs after normalization.

- **MHD CSV**: MHD motif residue boundaries
- **P-loop CSV**: P-loop motif residue boundaries
- **NB-ARC CSV**: NB-ARC domain boundaries (used to define CC domain extent)

## Metrics Reference

### Penalized Metrics (LCB/UCB)

These metrics use confidence bounds to reward consistency across interfaces and replicates.

| Metric | Description | Penalty | Units |
|--------|-------------|---------|-------|
| **ipTM_LCB** | AlphaFold 3 interface predicted TM-score. Measures predicted quality of inter-chain contacts. Higher = better predicted interface. | LCB across replicates | Score (0–1) |
| **Sum_CONTACTS** | Total inter-protomer heavy-atom contacts (< 4.0 Å, PAE < 10 Å if available) between ring-adjacent protomers. Counted per adjacent interface, LCB within structure, then LCB across replicates. Higher = more inter-protomer interactions. | LCB within + LCB across | Count |
| **CONTACTS_HBOND** | Hydrogen bond contacts (N/O donor–acceptor, <3.5 Å heavy-atom proxy). Same penalization as Sum_CONTACTS. | LCB within + LCB across | Count |
| **CONTACTS_SALT_BRIDGE** | Salt bridge contacts (charged atom pairs: Asp/Glu COO⁻ ↔ Arg/Lys/His NH₃⁺, <4.0 Å). | LCB within + LCB across | Count |
| **CONTACTS_HYDROPHOBIC** | Hydrophobic contacts (apolar C–C pairs, <4.0 Å). | LCB within + LCB across | Count |
| **CONTACTS_DISULFIDE** | Disulfide bonds (Cys SG–SG, <2.5 Å). | LCB within + LCB across | Count |
| **CONTACTS_VDW** | Van der Waals contacts (all other atom pairs, <4.0 Å). | LCB within + LCB across | Count |
| **BSA_INT_PROTO** | Buried surface area between ring-adjacent protomers, (SASA₁ + SASA₂ − SASA₁₂)/2 (FreeSASA, Lee–Richards, element radii, 1.4 Å probe). Geometric only: not PAE-filtered. LCB within structure interfaces, LCB across replicates. Higher = larger interface. | LCB within + LCB across | Å² |
| **SD_THETA_ROT** | Standard deviation of inter-protomer rotational angles. Measures rotational symmetry: 0° = perfectly symmetric. | UCB across replicates | Degrees |
| **S_PROTO** | Symmetry RMSD of protomer centroids relative to a regular polygon. Measures spatial symmetry: 0 Å = perfectly symmetric. | UCB across replicates | Å |
| **D_APEX** | Mean pairwise N-terminal Cα distance across all protomers. Indicates pore aperture size at the funnel entrance. When a MADA HMM hit is available, uses the Met at (or closest upstream of) the MADA alignment start as the reference residue, correcting for spurious N-terminal gene model extensions. NA if no N-terminal coordinate is found. | UCB across replicates | Å |

### Non-Penalized Metrics (Mean ± SD)

These metrics are averaged across protomers within each structure, then across replicates.

| Metric | Description | Units |
|--------|-------------|-------|
| **THETA_APEX** | Tilt angle of the α1-helix relative to the central pore axis. 0° = parallel to pore (tight funnel), 90° = perpendicular (splayed). Filtered by α1 linearity (R² ≥ 0.9); reports NA for disordered/kinked α1. | Degrees |
| **L_APEX** | Length of the α1-helix in residue count. | Residues |
| **H_ABS** | Mean Kyte–Doolittle hydrophobicity over all residues of the CC-domain helices (DSSP helices starting before the NB-ARC start, each truncated at the NB-ARC start). Indicates hydrophobic core packing quality. | Hydrophobicity index |
| **MU_H** | Amphipathic moment of the α1-helix (Eisenberg consensus hydrophobicity scale, 100° angular offset). Higher = more amphipathic (one hydrophobic face, one hydrophilic). | Moment (a.u.) |
| **D_MHD_P** | Spatial distance between MHD motif centroid and P-loop centroid (Cα atoms). Measures NB-ARC domain compactness. Requires MHD and P-loop annotation CSVs. | Å |
| **D_CC_NB_LINKER** | Spatial distance between last CC-domain helix end and NB-ARC domain start (Cα). Measures CC-to-NB-ARC linker extension. | Å |
| **L_CC_NB_LINKER** | Residue count of the linker between last CC helix and NB-ARC start. | Residues |

### Quality Indicators

Reported for monitoring but not used as primary scoring metrics.

| Metric | Description | Units |
|--------|-------------|-------|
| **ALPHA1_LINEARITY** | R² of PCA on α1 Cα coordinates. Measures how well α1 fits a straight line: 1.0 = perfect helix, <0.9 = kinked/disordered. When R² < 0.9, THETA_APEX is set to NA for that protomer. | R² (0–1) |
| **ALPHA1_PLDDT** | Mean pLDDT (from Cα B-factors) of the α1-helix region. Indicates AF3 prediction confidence for the α1 region. | pLDDT (0–100) |

### Penalty Formulas

**LCB** (Lower Confidence Bound) = x̄ − k × (s / √n)
- Conservative lower estimate: rewards metrics that are both *high* and *consistent*
- Used for metrics where higher is better (contacts, BSA, ipTM)

**UCB** (Upper Confidence Bound) = x̄ + k × (s / √n)
- Conservative upper estimate: penalizes metrics that are *noisy* even when low
- Used for metrics where lower is better (symmetry deviations, distances)

**Two-level penalization** (for Sum_CONTACTS, BSA, contact types):
1. Within each structure: LCB across all adjacent interfaces
2. Across replicates: LCB of the within-structure LCBs

Default k = 1.96. With n = 3 seeds this is a variance penalty, not a 95 % confidence interval. Single replicate → no penalty applied.

### Contact Type Classification

Contact types are assigned based on heavy-atom distance proxies (suitable for AF3 models without explicit hydrogens):

| Type | Atoms | Cutoff | Notes |
|------|-------|--------|-------|
| Disulfide | Cys SG ↔ SG | <2.5 Å | Checked first (highest priority) |
| Salt bridge | Charged pairs (Asp/Glu O ↔ Arg/Lys/His N) | <4.0 Å | |
| Hydrogen bond | N/O ↔ N/O/S | <3.5 Å | Heavy-atom proxy |
| Hydrophobic | Apolar C ↔ C | <4.0 Å | Excludes backbone Cα, C, O, N |
| VdW | Everything else | < `--contact_cutoff` (4.0 Å) | Default category |

## Chain Ring Order

The pipeline determines true spatial adjacency by computing the angular position of each protomer around the ring, rather than relying on chain letter order. This is critical for PDB/RCSB structures where chain naming is arbitrary (e.g. NRC3 `9RI9` has chains B,D,G,I,J,L; Sr35 has chains A,C,E,G,I in non-adjacent alphabetical order).

**Algorithm**: Compute Cα centroid per chain → fit plane via PCA → project onto plane → sort by `atan2` angle → true ring-adjacent order.

**QC metrics** reported per structure:
- `ring_adjacency_distances`: centroid-to-centroid distances for each adjacent pair (should be uniform)
- `ring_adjacency_cv`: coefficient of variation of adjacent distances (>0.15 warns of unreliable ring detection)

For AlphaFold 3 outputs (chains A–F) the geometric order often, but not always, matches alphabetical order: in the main NRCH set, 138 of the 1,911 models are regular hexamers whose chain letters are not ring-adjacent. For experimental structures it often differs. Each summary row reports `RING_ADJ_CV_max`, the largest `ring_adjacency_cv` across seeds.

## α1-Helix Detection

The pipeline supports three modes for defining the α1-helix boundary, in priority order:

### 1. HMM-Based (recommended)

Provide `--hmm MADA.hmm` or `--hmm_domtbl MADA_domtbl.out`. The pipeline:
1. Extracts chain A sequences from each structure (fast CIF/PDB parser, threaded I/O)
2. Runs `hmmsearch` against the HMM model
3. Parses the domain table: for each protein, takes the best N-terminal hit (lowest `ali_from`, then best E-value)
4. Defines α1 as residues 1 to `ali_to`

This is sequence-based, so it works regardless of structural prediction quality.

**N-terminal correction**: Gene models sometimes have spurious insertions before the real MADA motif. When a MADA HMM hit is found, the pipeline identifies the "true" N-terminal as: (1) the Met at the HMM alignment start position, or (2) the closest Met upstream of the alignment start. This corrected residue is used as the reference for D_APEX distance calculations and as the α1-helix start. Proteins without an HMM hit use residue 1 as usual.

### 2. DSSP-Based with Kink Splitting (fallback)

Active when no HMM is provided, or with `--hmm_fallback_dssp` for proteins without an HMM hit. Uses DSSP secondary structure to find the first CC-domain helix. If this helix is unusually long (DSSP merged α1+α2 due to a short helical linker), **kink detection** splits it at the point of maximum backbone direction change (>30° angular threshold, measured via PCA on sliding windows of 5 Cα atoms).

### 3. None (default when HMM misses)

When HMM is provided but has no hit and `--hmm_fallback_dssp` is off, α1-dependent metrics (THETA_APEX, L_APEX, MU_H) are reported as NA.

## Output

```
output_dir/
├── resistosome_analysis_summary.xlsx   # Main summary table (all metrics; also RING_ADJ_CV_max, PIPELINE_VERSION, GIT_COMMIT)
├── resistosome_analysis_summary.csv    # Same in CSV format
├── pipeline.log                        # Full processing log
├── hmm/                                # HMM analysis (if --hmm used)
│   ├── sequences.fasta                 # Extracted chain A sequences
│   ├── hmmsearch.out                   # hmmsearch full output
│   └── hmmsearch_domtbl.out            # Domain table (reusable via --hmm_domtbl)
└── per_structure_json/                 # Detailed per-structure results
    ├── protein_a_seed1_results.json
    ├── protein_a_seed2_results.json
    ├── protein_a_aggregated.json       # Cross-replicate aggregation
    └── ...
```

### Per-Structure JSON Fields

Each JSON includes:
- `pipeline_version`, `git_commit` and `parameters` (contact cutoff, PAE cutoff, k) used for the run
- `alpha1_source`: "hmm", "dssp", or "none"
- Per-interface: contact counts (total and by type), BSA, PAE statistics
- Per-protomer: THETA_APEX, L_APEX, MU_H, H_ABS, α1 linearity/pLDDT, linker metrics
- Within-structure summaries: means, SDs, LCBs
- `cc_helix_info`: per-chain CC-domain helix metadata including boundaries, sequences, and (for HMM α1) linearity R² and mean pLDDT
- `contact_type_summary`: breakdown by contact type per interface

### Reusing HMM Results

After a first run with `--hmm`, the generated `hmmsearch_domtbl.out` can be reused:
```bash
# First run (extracts sequences, runs hmmsearch)
python run_pipeline.py --hmm MADA.hmm -i input/ -o results/

# Subsequent runs (skip hmmsearch, much faster)
python run_pipeline.py --hmm_domtbl results/hmm/hmmsearch_domtbl.out -i input/ -o results2/
```

## Performance

Structure analysis is parallelized across all available CPU cores.

```bash
# Use all available cores
python run_pipeline.py --workers 14 -i input/ -o results/

# Serial mode (for debugging)
python run_pipeline.py --workers 1 -i input/ -o results/
```

Per-structure processing time: ~8–12s (6-protomer, dominated by DSSP and FreeSASA).

## SLURM Usage

See `submit_slurm.sh` for an HPC submission template.

## Slim AF3 inputs (`tools/`)

Full AF3 outputs are large because of the PAE matrix in `*_confidences.json`. The pipeline only uses the PAE to filter contacts at `--pae_cutoff`, so the per-model inputs can be reduced to about 0.6 MB without changing any result. The v2.0.0 main-set run was done this way.

| Script | Purpose |
|---|---|
| `tools/slim_af3.py OUT.tar 10.0 *.tar.gz` | Keep `model.cif` (gzipped), `summary_confidences.json` and, for every protein-chain pair, a packed boolean mask of min(PAE<sub>ij</sub>, PAE<sub>ji</sub>ᵀ) < cutoff. Needs NumPy only. |
| `tools/unpack_slim.py SLIM_DIR UNPACKED_DIR` | Unpack slim tars into AF3-like folders. |
| `tools/run_pipeline_slim.py <pipeline_dir> <run_pipeline.py args>` | Run the pipeline on unpacked slim folders. It rebuilds the PAE filter from the masks. |
| `tools/aggregate_v2_json.py --pipeline_dir . --out OUT json_dir1 [json_dir2 …]` | Merge per-structure JSONs from several partial runs into one summary, using the pipeline's own aggregation. |

The slim and full inputs gave identical values for every metric and every per-interface contact and BSA value (tested on 3 models). Slim inputs are valid only for the PAE cutoff used when slimming.

## Notes

- **Chain adjacency** is determined geometrically (angular sorting of centroids), not alphabetically — safe for both AF3 and PDB/RCSB structures
- Protomer count is detected automatically (supports 4, 5, 6, 8, etc.)
- Missing annotations → affected metrics are NA (pipeline continues normally)
- Annotation files may contain more entries than structures — only matching IDs are used
- Single replicate → no confidence penalty applied (LCB = UCB = raw value)
- `--contact_cutoff`, `--pae_cutoff` and `--k` are passed explicitly to every worker process (in v1 they did not reach spawned workers)
- The pipeline handles both compressed (.tar.gz) and uncompressed AF3 output directories
- Standalone ring order tool: `python get_ring_order.py structure.cif` (requires gemmi)
