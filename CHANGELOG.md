# Changelog

## v2.0.0 (November 2026)

The SNI pipeline was corrected and both datasets were re-run: the main NRCH set (1,911 models, 637 NLRs) and the benchmark set (54 models, 18 NLRs). All summary tables, supplementary data and figures in this repository are from v2.0.0.

### Pipeline fixes (`resistosome_pipeline/`)

| # | Issue | Effect in v1 | Fix in v2 |
|---|---|---|---|
| 1 | **Ring adjacency.** The published main-set summary was produced before geometric ring ordering was added, so protomers were paired by chain letter (A–B, B–C, …, F–A). | 138 models (118 NLRs) are regular hexamers whose AF3 chain letters are not ring-adjacent. For these, contacts and BSA were counted between non-neighbours and σθ ROT was inflated. 126 NLRs change in v2 (the 118 plus irregular assemblies). Median change in the mis-ordered NLRs: Σ_CONTACTS +76, BSA +1,121 Å², σθ ROT −28°. | Protomers are always ordered by their angle around the ring axis. The QC value `RING_ADJ_CV_max` is reported per NLR. |
| 2 | **CC-domain helices ran into the NB-ARC domain.** Helices starting before the NB-ARC boundary were kept whole. | NB-ARC residues were counted in H_ABS and the CC–NB-ARC linker metrics. | Helices are truncated at the NB-ARC start. H_ABS changes by ≤ 0.17 (236 NLRs). |
| 3 | **`--contact_cutoff` had no effect.** The search radius was always max(cutoff, 4.0 Å) and vdW used a fixed 4.0 Å, while the log reported 3.5 Å. | The cutoff actually used was 4.0 Å. | The search radius is the cutoff; the default is now 4.0 Å. Results are unchanged. |
| 4 | **`--contact_cutoff`, `--pae_cutoff` and `--k` did not reach worker processes.** Workers always used the defaults. | — | All three are passed explicitly. |
| 5 | **D_APEX fell back to 0.0** when no N-terminal coordinate was found. | No instances in either dataset. | Returns NA. |
| 6 | **`pipeline.log` location and provenance.** The log could be written outside the output folder. | — | The output folder is created first. The pipeline version, git commit and parameters are written to the log, every JSON and the summary table. |

BSA is unchanged: it is geometric and not PAE-filtered, as in v1.

### Unchanged values

ipTM_LCB, S_PROTO, D_APEX, φ_APEX, L_APEX, μH and D_MHD-P are bit-identical between v1 and v2 for all 637 NLRs.

In the benchmark set only H_ABS changes, because v1 had already been run with ring ordering. `R/figures/v1_vs_v2_*` give the per-NLR comparison.

### Figures (`R/`)

`NRC_heatmap_analysis.R` is replaced by `SNI_figures.R`.

- **H_ABS was forced to 0.** v1 applied `pmax(x, 0)` to every parameter before computing clade means, so the H_ABS rows of Fig 1A, 2A, S6 and S7 (and `H_ABS_mean` in Data_S5) were 0. Now only the penalised contact, BSA and ipTM bounds are clamped at 0.
- **Clade membership.** Tree tips and pipeline IDs are matched after mapping every non-alphanumeric character to "_". This recovers the *S. tuberosum* E86-69 entries, whose tips contain "-". 622 of the 637 NLRs fall in the 18 clades, as in the published Data_S5.
- **Rendering.** All panels are drawn in ggplot2 and assembled with patchwork on shared coordinates (trees, dendrograms, tiles and labels align). Text is set in Helvetica Neue, and parameters are shown with their mathematical symbols.
- **New panels.** Fig 2B and S4 include the NB-ARC tree.
- **New figures.** `*_2way` heatmaps also cluster the parameters. `ParamClust_*` show the parameter correlation heatmap and a bootstrapped parameter dendrogram (pvclust, 1,000 resamples).
- **Diagnostics.** `v1_vs_v2_*` compare the two versions.

### Supplementary data (`supplementary/`)

| File | Change |
|---|---|
| Data_S1.xlsx (metric definitions) | Contacts and σθ ROT use ring-adjacent pairs and the 4.0 Å cutoff. BSA is computed with FreeSASA (Lee–Richards) and is not PAE-filtered. D_APEX uses the MADA-corrected N-terminus and returns NA on failure. L_APEX is in residues (not Å). H_ABS and linker helices are truncated at the NB-ARC start. LCB/UCB use SD/√n. |
| Data_S2.xlsx (benchmark summary) | Regenerated with v2.0.0. |
| Data_S4.xlsx (main summary) | Regenerated with v2.0.0. |
| Data_S5.csv (clade means) | Regenerated: real H_ABS values, corrected Σ_CONTACTS, BSA and σθ ROT. |

Data_S3, Data_S6, Data_extS1 and the `NRC_ref*` files are unchanged.

### Documentation

`SNI_methods_detailed.md` is updated with:

- the BSA formula (÷ 2) and the SASA algorithm (Lee–Richards, the FreeSASA default, not Shrake–Rupley);
- ring ordering;
- helix truncation;
- D_APEX = NA on failure;
- clade-level summaries and clustering (new §10);
- corrected file paths.

### Effect on the results

- **Fig 1 / S6 (benchmark).** Only H_ABS changes. The conclusions stand:
  - the combined SNI separates the 9 helpers from the 9 sensors;
  - scientist parameters alone misplace NbNRC4c;
  - co-scientist parameters alone misplace Rpa1.
- **Fig 2A (clades).**
  - **σθ ROT** in v1 was inflated by mis-ordering; it was up to 23° for NRC4d and 43° for NRC4u. In v2 it is below 3° in 14 of 18 clades. The exceptions are NRC4u (17.1°), NRC9 (7.2°), NRC7 (5.5°) and NRC5u (4.2°).
  - **Ward clustering, k = 2.** NRC4u, NRC5u, NRC6 and NRC7 fall outside the canonical cluster. NRC5u is new: it has 5 members, low ipTM and few contacts.
- **Fig S7 (parameter subsets).**
  - **Scientist parameters.** NRC4u, NRC6 and NRC7 split off together at k = 2, and NRC4u separates from NRC6/NRC7 at k = 4.
  - **Co-scientist parameters.** NRC4u is the only clear outlier. NRC6 and NRC7 group with NRC0 at k = 4, but not on their own.
- **Parameter clustering (622 NLRs).** Three groups are supported:
  - interface: ipTM, Σ_CONTACTS and BSA (AU/BP 100/100);
  - ring geometry: S_PROTO + σθ ROT (100/100) and L_APEX + D_APEX (99/98);
  - CC/α1 descriptors (φ_APEX, μH, H_ABS, D_MHD-P), which are nearly uncorrelated with the rest.

### Reproducibility notes

- The v2.0.0 main-set run used slim AF3 inputs (`resistosome_pipeline/tools/`). Slim and full inputs give identical outputs.
- DSSP 4.2.2 was used for the rerun instead of 4.5.8. Every metric that v2 does not change was bit-identical.
