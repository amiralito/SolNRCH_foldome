# Supporting scripts and material for "AI-guided discovery of atypical protein assemblies"

[![DOI](https://img.shields.io/badge/Zenodo-10.5281/zenodo.XXXXXXX-blue.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) [![DOI](https://img.shields.io/badge/bioRxiv-doi.org/10.1101/2026.XX.XX.XXXXXX-BE2634.svg)](https://doi.org/10.1101/2026.XX.XX.XXXXXX)

AmirAli Toghani<sup>†</sup>, Benjamin A. Seager<sup>†</sup>, Yu Sugihara, Lisa-Marie Roijen, Juan M. Azcue, Maián Garro, Maryam Sargolzaei, Ioanna Morianou, Adeline Harant, Sam Gallop, Jiorgos Kourelis, Dan MacLean, Mauricio P. Contreras, Sophien Kamoun\*, Daniel Lüdke\*

<sup>†</sup> These authors contributed equally to this work.

\*Corresponding authors.

This repository contains the scripts, intermediate data, and detailed methods used for the structural analysis of the Solanaceae NRCH (NRC-helper) clade and the calculation of the Structural Novelty Index (SNI) from AlphaFold 3 hexameric resistosome predictions.

## Repository layout

```
.
├── SNI_methods_detailed.md         Full SNI methods with equations and software versions
├── gcloud_homomer/                 AlphaFold 3 homomeric inference pipeline (Google Cloud)
│   ├── README.md                   Pipeline usage and option reference
│   ├── AF3_homomeric_submit.sh     Main entry point
│   ├── AF3_homomeric_job.sh        SLURM job script
│   ├── AF3_homomeric_controller.sh Controller for large runs (>1000 jobs)
│   ├── generate_af3_homomeric.py   Generate AF3 input JSONs (n protomers, ligands, seeds)
│   └── extract_af3_scores.py       Extract per-model and per-chain scores from AF3 output
├── resistosome_pipeline/           Python pipeline: per-structure SNI metrics
│   ├── README.md                   Pipeline usage and command-line reference
│   ├── requirements.txt            Pinned Python dependencies
│   ├── run_pipeline.py             Command-line entry point
│   ├── submit_slurm.sh             SLURM submission template (Kronos HPC)
│   ├── test_pipeline.py            Smoke tests
│   ├── get_ring_order.py           Helper: derive ring protomer order
│   └── resistosome/                Pipeline modules (af3, contacts, geometry, helix, …)
├── structural_analysis/            Pipeline inputs and per-structure outputs
│   ├── main/                       NRCH dataset (407 hexamers across 18 clades)
│   │   ├── MADA.hmm                MADA HMM profile (used for α1 detection)
│   │   ├── *_coords.csv            MHD / P-loop / NB-ARC reference coordinates
│   │   └── output/
│   │       ├── hmm/                HMMER 3.4 outputs (hmmsearch.out, domtbl, fasta)
│   │       ├── pipeline.log
│   │       ├── resistosome_analysis_summary.csv
│   │       └── resistosome_analysis_summary.xlsx
│   └── test/                       Reference set (characterised NLRs)
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
├── R/                              Data preparation + heatmap analysis (R)
│   ├── SolNRCH_foldome_v1.rmd      Sequence/metadata preparation pipeline (NLRtracker
│   │                               -> NB-ARC filtering -> CC-NB-ARC extraction for AF3
│   │                               -> 95 % similarity cut -> motif coordinate export)
│   ├── NRC_heatmap_analysis.R      Per-clade SNI heatmaps + parameter-subset heatmaps
│   ├── main/                       Main set: summary xlsx + heatmap outputs
│   │   ├── resistosome_analysis_summary.xlsx
│   │   ├── heatmaps/               main_clade_heatmap{,_zraw}.{pdf,png,svg}
│   │   └── comparison/             main_{scientist,coscientist}_params.{pdf,png,svg}
│   └── test/                       Test set: summary xlsx + heatmap outputs
│       ├── resistosome_analysis_summary.xlsx
│       ├── helpers.txt             Reference helper NLR IDs
│       ├── sensors.txt             Reference sensor NLR IDs
│       ├── heatmaps/
│       └── comparison/
└── supplementary/                  Manuscript supplementary data tables
    ├── Data_S1.xlsx                NLRtracker-derived NLR metadata for the Solanaceae set
    ├── Data_S2.xlsx                Domain-architecture summary
    ├── Data_S3.csv                 Per-entry NB-ARC sequence table
    ├── Data_S4.xlsx                NRCH clade assignments and per-entry annotations
    ├── Data_S5.csv                 Reference NLR coordinates
    └── Data_S6.zip                 Bundled per-clade alignments and trees
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
| _DSSP v4.5.8_ ([https://github.com/PDB-REDO/dssp](https://github.com/PDB-REDO/dssp)) |
| _AlphaFold 3 server_ ([https://alphafoldserver.com/](https://alphafoldserver.com/)) |
| _NLRtracker_ ([https://github.com/slt666666/NLRtracker](https://github.com/slt666666/NLRtracker)) |
| _NLRexpress_ ([https://github.com/eliza-m/NLRexpress](https://github.com/eliza-m/NLRexpress)) |
| _FAMSA v2.2.2_ ([https://github.com/refresh-bio/FAMSA](https://github.com/refresh-bio/FAMSA)) |
| _FastTree v2.1.10_ ([http://www.microbesonline.org/fasttree/](http://www.microbesonline.org/fasttree/)) |
| _Dendroscope v3.8.8_ ([https://software-ab.cs.uni-tuebingen.de/download/dendroscope3/welcome.html](https://software-ab.cs.uni-tuebingen.de/download/dendroscope3/welcome.html)) |
| _R v4.4.2_ ([https://cran.r-project.org/](https://cran.r-project.org/)) |

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
install.packages("tidyverse")
install.packages("readxl")
install.packages("ape")
install.packages("pheatmap")
install.packages("reshape2")
install.packages("svglite")

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
   See `resistosome_pipeline/README.md` for the full command-line reference and `SNI_methods_detailed.md` for the equations, residue-level thresholds, and contact-classification rules.

## Supplementary Data

- AlphaFold 3 hexameric predictions for the main NRCH set and the per-structure JSON outputs are deposited on Zenodo: [DOI: 10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX).
- Detailed SNI methods: [`SNI_methods_detailed.md`](SNI_methods_detailed.md).
