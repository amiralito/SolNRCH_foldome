# Supporting scripts and material for "Discovery of atypical NLR resistosomes"

[![DOI](https://img.shields.io/badge/Zenodo-10.5281/zenodo.XXXXXXX-blue.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) [![DOI](https://img.shields.io/badge/bioRxiv-doi.org/10.64898/2026.05.03.722499-BE2634.svg)](https://doi.org/10.64898/2026.05.03.722499)

AmirAli Toghani<sup>†</sup>, Benjamin A. Seager<sup>†</sup>, Yu Sugihara, Lisa-Marie Roijen, Juan M. Azcue, Maián Garro, Maryam Sargolzaei, Ioanna Morianou, Adeline Harant, Sam Gallop, Jiorgos Kourelis, Dan MacLean, Mauricio P. Contreras, Sophien Kamoun\*, Daniel Lüdke\*

<sup>†</sup> These authors contributed equally to this work.

\*Corresponding authors.

This repository contains the scripts, intermediate data, and detailed methods used for the structural analysis of the Solanaceae NRCH (NRC-helper) clade and the calculation of the Structural Novelty Index (SNI) from AlphaFold 3 hexameric resistosome predictions.

> **Version 2.0.0.** The SNI pipeline was corrected (geometric ring adjacency in the published main-set run, CC-domain helices truncated at the NB-ARC start, contact-cutoff handling) and both datasets were re-run. All tables, supplementary data and figures here are from v2.0.0. See [`CHANGELOG.md`](CHANGELOG.md) for the fixes and their effect on the results.

## Repository layout

```
.
├── CHANGELOG.md                    v1 → v2.0.0 fixes and their effect on the results
├── SNI_methods_detailed.md         Full SNI methods with equations and software versions
├── gcloud_homomer/                 AlphaFold 3 homomeric inference pipeline (Google Cloud)
│   ├── README.md                   Pipeline usage and option reference
│   ├── AF3_homomeric_submit.sh     Main entry point
│   ├── AF3_homomeric_job.sh        SLURM job script
│   ├── AF3_homomeric_controller.sh Controller for large runs (>1000 jobs)
│   ├── generate_af3_homomeric.py   Generate AF3 input JSONs (n protomers, ligands, seeds)
│   └── extract_af3_scores.py       Extract per-model and per-chain scores from AF3 output
├── resistosome_pipeline/           Python pipeline (v2.0.0): per-structure SNI metrics
│   ├── README.md                   Pipeline usage and command-line reference
│   ├── requirements.txt            Pinned Python dependencies
│   ├── run_pipeline.py             Command-line entry point
│   ├── submit_slurm.sh             SLURM submission template (Kronos HPC)
│   ├── test_pipeline.py            Smoke tests
│   ├── get_ring_order.py           Helper: derive ring protomer order
│   ├── resistosome/                Pipeline modules (af3, contacts, geometry, helix, …)
│   └── tools/                      Slim AF3 inputs (model + PAE<10 Å masks) and JSON aggregation
├── structural_analysis/            Pipeline inputs and outputs (v2.0.0)
│   ├── main/                       NRCH dataset: 637 NLRs × 3 seeds = 1,911 hexamers (622 NLRs in 18 clades)
│   │   ├── MADA.hmm                MADA HMM profile (used for α1 detection)
│   │   ├── *_coords.csv            MHD / P-loop / NB-ARC reference coordinates
│   │   └── output/
│   │       ├── hmm/                HMMER 3.4 outputs (hmmsearch.out, domtbl, fasta)
│   │       ├── per_structure_json.zip
│   │       ├── pipeline.log            Warnings from the 7 slim-input batches the main set was run in
│   │       ├── resistosome_analysis_summary.csv
│   │       └── resistosome_analysis_summary.xlsx
│   └── test/                       Benchmark set: 18 characterised NLRs (9 helpers, 9 sensors) × 3 seeds
│       ├── MADA.hmm
│       ├── ids.txt
│       ├── merged_*_cooreds.csv
│       └── output/                 Per-structure JSON + summary tables
├── phylo/                          Phylogenetic analysis of the NRC-H clade
│   ├── NRCH_clade.tree                               Initial NRC-H clade tree
│   ├── NRCH_clade_filtered_len_seq*.fasta            NRC-H sequence sets
│   ├── NRCH_clade_filtered_len_seq_95.fasta.clstr    CD-HIT 95 % cluster table
│   ├── NRCH_clade_filtered_len_seq_95_ids.txt        Final 637-entry ID list
│   ├── NRCH_clade_filtered_len_seq_CCNBARC_filtered_*.{fasta,afa,newick,tree}
│   │                                                 NB-ARC alignments and trees
│   ├── *_rooted_itol.txt                             iTOL annotation file
│   └── clades/                                       18 per-clade Newick subtrees
├── R/                              Data preparation + SNI figures (R)
│   ├── SolNRCH_foldome_v1.rmd      Sequence/metadata preparation pipeline (NLRtracker
│   │                               -> NB-ARC filtering -> CC-NB-ARC extraction for AF3
│   │                               -> 95 % similarity cut -> motif coordinate export)
│   ├── SNI_figures.R               All SNI figures, clade summaries and parameter clustering
│   ├── figures/                    Fig 1A, 2A, 2B, S4, S6A/B, S7A/B (+ *_2way), ParamClust_*,
│   │                               v1_vs_v2_*; each as PDF/PNG/SVG, every heatmap with
│   │                               *_raw_values.csv and *_zscores.csv
│   └── v1_reference/               v1 summary tables (input for the v1_vs_v2 diagnostics)
└── supplementary/                  Manuscript supplementary data tables
    ├── Data_S1.xlsx                SNI metric definitions, calculations and aggregation rules
    ├── Data_S2.xlsx                Benchmark set: per-NLR SNI summary (v2.0.0)
    ├── Data_S3.csv                 NRC-H sequences and NLRtracker annotations (637 entries)
    ├── Data_S4.xlsx                Main NRCH set: per-NLR SNI summary (v2.0.0)
    ├── Data_S5.csv                 Clade means of the 11 SNI parameters (v2.0.0)
    ├── Data_S6.zip                 AI co-scientist hypotheses and reports
    ├── Data_extS1.xlsx             Solanaceae genome/proteome assemblies used (346 + 11 deprecated)
    └── NRC_ref*.fasta(.clstr)      Reference NRC sequences and clusters
```

The full set of unpacked AlphaFold 3 predictions and the per-structure JSON outputs that feed `resistosome_pipeline` are deposited separately (see [Supplementary Data](#supplementary-data)).

## Resources

### Software

| Source |
|--------|
| _Python v3.12_ ([https://www.python.org/](https://www.python.org/)) |
| _BioPython v1.84_ ([https://biopython.org/](https://biopython.org/)) |
| _SciPy v1.15.2_ ([https://scipy.org/](https://scipy.org/)) |
| _NumPy v1.26.3_ ([https://numpy.org/](https://numpy.org/)) |
| _pandas v2.3.3_ ([https://pandas.pydata.org/](https://pandas.pydata.org/)) |
| _openpyxl v3.1.5_ ([https://openpyxl.readthedocs.io/](https://openpyxl.readthedocs.io/)) |
| _FreeSASA v2.2.1_ ([https://freesasa.github.io/](https://freesasa.github.io/)) |
| _HMMER v3.4_ ([http://hmmer.org/](http://hmmer.org/)) |
| _DSSP v4.5.8_ (v4.2.2 for the v2.0.0 rerun) ([https://github.com/PDB-REDO/dssp](https://github.com/PDB-REDO/dssp)) |
| _AlphaFold 3 server_ ([https://alphafoldserver.com/](https://alphafoldserver.com/)) |
| _NLRtracker_ ([https://github.com/slt666666/NLRtracker](https://github.com/slt666666/NLRtracker)) |
| _FAMSA v2.2.2_ ([https://github.com/refresh-bio/FAMSA](https://github.com/refresh-bio/FAMSA)) |
| _FastTree v2.1.10_ ([http://www.microbesonline.org/fasttree/](http://www.microbesonline.org/fasttree/)) |
| _Dendroscope v3.8.8_ ([https://software-ab.cs.uni-tuebingen.de/download/dendroscope3/welcome.html](https://software-ab.cs.uni-tuebingen.de/download/dendroscope3/welcome.html)) |
| _R v4.4.2_ (data preparation); _R v4.3.3_ (SNI figures) ([https://cran.r-project.org/](https://cran.r-project.org/)) |
| _pvclust v2.2_ ([https://github.com/shimo-lab/pvclust](https://github.com/shimo-lab/pvclust)) |

### Python packages

```bash
# Pinned versions (exact reproduction of the analysis):
pip install -r resistosome_pipeline/requirements.txt

# System packages:
# DSSP 4.5.8
sudo apt-get install dssp           # Ubuntu/Debian
# or: conda install -c salilab dssp

# HMMER 3.4
conda install -c bioconda hmmer     # or: brew install hmmer
```

### R packages

```r
# SNI figures (R/SNI_figures.R)
install.packages(c("tidyverse", "ape", "patchwork", "svglite", "pvclust"))

# Data preparation (R/SolNRCH_foldome_v1.rmd)
install.packages(c("readxl", "reshape2"))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("ggtree")
```

[_NLRtracker_](https://github.com/slt666666/NLRtracker) and [_NLRexpress_](https://github.com/eliodelmarre/NLRexpress) are used upstream of `SolNRCH_foldome_v1.rmd` to extract NLR domain annotations and motif coordinates; install per their own instructions.

## Description
**Sequence selection and motif coordinate extraction** &mdash; `R/SolNRCH_foldome_v1.rmd` documents how the NRC-H clade was extracted from NLRtracker output, filtered for length and architecture, deduplicated at 95 % similarity, sliced to CC-NB-ARC domains for AlphaFold 3 input, and how the MHD/P-loop coordinates used downstream were exported from NLRexpress. The intermediate alignments and trees are deposited under `phylo/`.

**AlphaFold 3 homomeric modelling** &mdash; the SLURM-based pipeline under `gcloud_homomer/` was used to predict hexameric resistosomes for every NRC-H entry on a Google Cloud A100 cluster. Per-protein inference inputs (MSAs and templates) were pre-computed for each protomer using the AlphaFold 3 data pipeline ([docs](https://github.com/google-deepmind/alphafold3/blob/main/docs/performance.md#data-pipeline)). Starting from a folder of these single-protomer JSONs, the pipeline generates AF3 input JSONs for the homomeric assembly, expands seeds and (optional) ligand copies, submits one AF3 inference job per (protein, seed) pair, and extracts per-model and per-chain confidence scores from `summary_confidences.json`:

   ```bash
   cd gcloud_homomer
   bash AF3_homomeric_submit.sh ./inference_input ./output \
       --protomers 6 \
       --seeds 1 2 \
       --expand-seeds \
       --batch NRCH_hexamers
   ```

   This produces a per-job tarball, `scores/*_scores.tsv` (overview metrics) and `scores/*_matrix.tsv` (per-chain pairwise matrices) under the chosen `output/` directory. Failed jobs are logged to `logs/failed_jobs.tsv` for re-submission. See `gcloud_homomer/README.md` for the full option reference, output naming convention, and memory considerations (the per-job token budget for an 80 GB A100). The unpacked tarballs from this step are deposited on Zenodo and are the input to the SNI metrics step below.

**Per-structure SNI metrics** &mdash; from a folder of unpacked AlphaFold 3 hexamer predictions (downloaded from the Zenodo deposit, or produced by the AF3 step above), compute the per-structure metric table:

   ```bash
   cd resistosome_pipeline
   python run_pipeline.py \
       --input_dir <path to unpacked AF3 predictions> \
       --output_dir ../structural_analysis/main/output \
       --hmm ../structural_analysis/main/MADA.hmm \
       --mhd_csv  ../structural_analysis/main/NRCH_MHD_nlrexp_coords.csv \
       --ploop_csv ../structural_analysis/main/NRCH_ploop_nlrexp_coords.csv \
       --nbarc_csv ../structural_analysis/main/NRCH_clade_filtered_len_NBARC_filtered_95_coords.csv \
       --workers 12
   ```

   This produces `output/per_structure_json/*.json`, `output/hmm/*`, `pipeline.log`, and the aggregated `resistosome_analysis_summary.{csv,xlsx}`.
   Contacts use a 4.0 Å heavy-atom cutoff (`--contact_cutoff`, default) and PAE < 10 Å. The run takes about 6–10 s per model per worker.
   See `resistosome_pipeline/README.md` for the full command-line reference and `SNI_methods_detailed.md` for the equations, residue-level thresholds, and contact-classification rules.
   To save disk space, the AF3 tarballs can first be reduced to the model, summary confidences and PAE < 10 Å masks (about 0.6 MB per model) with `resistosome_pipeline/tools/`. This gives identical results. See the pipeline README.

**SNI figures, clade summaries and parameter clustering** &mdash; from the repository root:

   ```bash
   Rscript R/SNI_figures.R --base . \
       --main structural_analysis/main/output/resistosome_analysis_summary.csv \
       --benchmark structural_analysis/test/output/resistosome_analysis_summary.csv \
       --trees phylo \
       --v1_main R/v1_reference/main_v1_summary.csv \
       --v1_benchmark R/v1_reference/benchmark_v1_summary.csv \
       --out R/figures --nboot 1000 \
       --font /System/Library/Fonts/HelveticaNeue.ttc
   ```

   This produces the clade and benchmark heatmaps (Fig 1A, 2A, S6, S7, each also as a two-way clustered version), the tree-aligned contact and SNI panels (Fig 2B, S4), and the parameter correlation and bootstrapped parameter dendrograms (`ParamClust_*`). Each heatmap also gets its raw and z-scored values as CSV.
   Clade membership comes from `phylo/clades/*.tree`. Penalised contact, BSA and ipTM bounds are set to max(x, 0) before clade means are taken. `--font` is optional; without it, PDFs reference Helvetica Neue by name.

## Supplementary Data

- AlphaFold 3 hexameric predictions for the main NRCH set and the per-structure JSON outputs (v2.0.0) are deposited on Zenodo: [DOI: 10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX).
- Changes from the first preprint version: [`CHANGELOG.md`](CHANGELOG.md).
- Detailed SNI methods: [`SNI_methods_detailed.md`](SNI_methods_detailed.md).
