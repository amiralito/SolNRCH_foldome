# Supporting scripts and material for "Discovery of atypical NLR resistosomes"

[![DOI](https://img.shields.io/badge/Zenodo-10.5281/zenodo.XXXXXXX-blue.svg)](https://doi.org/10.5281/zenodo.XXXXXXX) [![DOI](https://img.shields.io/badge/bioRxiv-doi.org/10.1101/2026.XX.XX.XXXXXX-BE2634.svg)](https://doi.org/10.1101/2026.XX.XX.XXXXXX)

AmirAli Toghani\*

\*corresponding author

This repository contains the scripts, intermediate data, and detailed methods used for the structural analysis of the Solanaceae NRCH (NRC-helper) clade and the calculation of the Structural Novelty Index (SNI) from AlphaFold 3 hexameric resistosome predictions.

## Repository layout

```
.
├── SNI_methods_detailed.md         Full SNI methods with equations and software versions
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
└── R/                              Heatmap analysis (R)
    ├── NRC_heatmap_analysis.R      Per-clade SNI heatmaps + parameter-subset heatmaps
    ├── main/                       Main set: summary xlsx + heatmap outputs
    │   ├── resistosome_analysis_summary.xlsx
    │   ├── heatmaps/               main_clade_heatmap{,_zraw}.{pdf,png,svg}
    │   └── comparison/             main_{scientist,coscientist}_params.{pdf,png,svg}
    ├── test/                       Test set: summary xlsx + heatmap outputs
    │   ├── resistosome_analysis_summary.xlsx
    │   ├── helpers.txt             Reference helper NLR IDs
    │   ├── sensors.txt             Reference sensor NLR IDs
    │   ├── heatmaps/
    │   └── comparison/
    └── trees/                      Clade definitions (used to group main-set entries)
        ├── NRCH_clade_filtered_len_seq_CCNBARC_filtered_95_NBARC_rooted_itol.txt
        └── clades/                 18 per-clade Newick tree files
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
install.packages("svglite")
```

## Reproducing the analysis

1. **Per-structure SNI metrics** &mdash; from a folder of unpacked AlphaFold 3 hexamer predictions (downloaded from the Zenodo deposit), compute the per-structure metric table:

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

2. **Per-clade heatmaps** &mdash; from the aggregated summary tables, generate the figure-quality heatmaps:

   ```bash
   cd R
   Rscript NRC_heatmap_analysis.R
   ```

   Two datasets are processed automatically: the *main* set (entries grouped by NRCH clade from `trees/clades/*.tree`) and the *test* set (each entry treated as its own group, with NRC0 and SlNRC0-Sa excluded). Outputs are written as PDF, PNG, and SVG into `R/{main,test}/heatmaps/` and `R/{main,test}/comparison/`.

## Supplementary Data

- Unpacked AlphaFold 3 hexameric predictions for the main NRCH set are deposited on Zenodo: [DOI: 10.5281/zenodo.XXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXX).
- Detailed SNI methods: [`SNI_methods_detailed.md`](SNI_methods_detailed.md).

