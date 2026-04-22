# Structural Novelty Index (SNI) — Detailed Methods

This document accompanies the *Discovery of atypical NLR resistosomes* manuscript and the associated code repository. It provides the full per-metric definitions, equations, thresholds, and aggregation rules summarised in the main Methods section. The implementation is distributed under `resistosome_pipeline/` in the project repository (<https://github.com/amiralito/SolNRCH_foldome>). File paths referenced below (e.g. `resistosome_pipeline/README.md`) are relative to the repository root.

All metrics and scores were computed for each structure (i.e. each seed replicate) and then aggregated across replicates as described in *Cross-replicate aggregation* below.

## 1. Interface confidence

**ipTM.** The interface predicted TM-score was extracted directly from the AlphaFold 3 summary confidence JSON.

## 2. Inter-protomer contacts

Contacts were enumerated between all pairs of ring-adjacent protomers (pairs (1,2), (2,3), …, (n,1) in the geometrically determined ring order). For each adjacent pair, a k-d tree spatial index (`scipy.spatial.cKDTree`) was constructed from heavy-atom coordinates of the second chain. Heavy-atom pairs within a search radius of 4.0 Å were identified, and each contact was classified into one of five types by the following priority rules:

- **Disulfide bond:** both atoms are Cys Sγ, distance < 2.5 Å.
- **Salt bridge:** one atom carries a formal positive charge (Arg NH1/NH2/Nε, Lys Nζ, His Nδ1/Nε2) and the other a formal negative charge (Asp Oδ1/Oδ2, Glu Oε1/Oε2), distance < 4.0 Å.
- **Hydrogen bond:** one atom is a potential donor (N or O) and the other a potential acceptor (N, O, or S), distance < 3.5 Å. This heavy-atom distance criterion serves as a proxy for hydrogen-bond geometry without explicit hydrogen placement.
- **Hydrophobic contact:** both atoms are carbon, and at least one belongs to a hydrophobic or aromatic residue (Ala, Val, Leu, Ile, Met, Phe, Trp, Pro, Tyr, His), or both are non-backbone side-chain carbons, distance < 4.0 Å.
- **Van der Waals:** any remaining heavy-atom pair within 4.0 Å not classified above.

When an AlphaFold 3 predicted aligned error (PAE) matrix was available, contacts between residue pairs with PAE ≥ 10 Å were excluded as low-confidence. Σ<sub>CONTACTS</sub> is the total number of contacts across all adjacent interfaces. Per-type totals (hydrogen bonds, salt bridges, hydrophobic contacts, disulfides, van der Waals) are reported separately.

## 3. Buried surface area

**BSA<sub>INTER-PROTO</sub>.** Buried surface area per inter-protomer interface was calculated as

$$\mathrm{BSA}_{ij} = \mathrm{SASA}(i) + \mathrm{SASA}(j) - \mathrm{SASA}(i \cup j)$$

where SASA(*i*) denotes the solvent-accessible surface area of chain *i*, computed with FreeSASA using the Shrake–Rupley algorithm.

## 4. Rotational symmetry

**σ<sub>θ ROT</sub>.** For each ring-adjacent pair, the angle subtended at the ring centroid by the two protomer centroids was computed. For an ideal *n*-mer, each angle equals 360°/*n* (60° for hexamers). σ<sub>θ ROT</sub> is the population standard deviation of these *n* angles:

$$\sigma_{\theta\ \mathrm{ROT}} = \sqrt{\tfrac{1}{n}\sum_{i=1}^{n}\left(\theta_i - \bar{\theta}\right)^2}$$

where θ<sub>i</sub> is the angle subtended by protomer pair (*i*, *i*+1 mod *n*) at the ring centroid.

**S<sub>PROTO</sub>.** A radial-symmetry metric computed as the standard deviation of protomer centroid-to-ring-centre distances. For a perfectly symmetric ring all protomer centroids lie at the same radius and S<sub>PROTO</sub> = 0; higher values indicate increasing deviation from rotational symmetry:

$$S_{\mathrm{PROTO}} = \sqrt{\tfrac{1}{n}\sum_{i=1}^{n}\left(d_i - \bar{d}\right)^2}$$

where d<sub>i</sub> = ‖**c**<sub>i</sub> − **c̄**‖ is the Euclidean distance from protomer centroid **c**<sub>i</sub> to the ring centre **c̄**. Here **c**<sub>i</sub> and **c̄** are three-dimensional position vectors, so ‖·‖ denotes the Euclidean (L2) norm of the vector difference, computed with `numpy.linalg.norm`; d<sub>i</sub> is therefore a scalar, and d̄ is the arithmetic mean of the d<sub>i</sub> across the *n* protomers.

## 5. Pore aperture

**D<sub>APEX</sub>.** The pore aperture was quantified as the mean of all pairwise Euclidean distances between N-terminal Cα atoms across protomers:

$$D_{\mathrm{APEX}} = \frac{2}{n(n-1)} \sum_{i \lt j} \lVert \mathbf{r}^{\mathrm{Nterm}}_i - \mathbf{r}^{\mathrm{Nterm}}_j \rVert$$

where **r**<sub>i</sub><sup>Nterm</sup> is the Cα coordinate of the reference N-terminal residue of protomer *i*. The sum runs over the *n*(*n*−1)/2 unordered pairs of protomers, so the leading factor 2/[*n*(*n*−1)] is the reciprocal of the pair count and yields the arithmetic mean. As above, ‖·‖ denotes the Euclidean (L2) norm of the difference between two three-dimensional position vectors, computed with `scipy.spatial.distance.pdist`.

## 6. CC-domain helix metrics

The pore axis of the resistosome ring was defined as the normal to the best-fit plane of protomer Cα centroids, obtained as the eigenvector corresponding to the smallest singular value from SVD of the mean-centred centroid matrix. In practice, the *n* × 3 matrix of protomer Cα centroids was mean-centred and decomposed using `numpy.linalg.svd`; the right-singular vector associated with the smallest singular value, which is normal to the least-squares plane of the centroids, was taken as the pore axis and normalised to unit length.

**φ<sub>APEX</sub>.** The tilt angle of the α1-helix relative to the pore axis was computed per protomer. The α1 helix axis was estimated as the first principal component (largest singular value) from PCA of the Cα coordinates within the α1 segment. The tilt angle was then

$$\varphi_{\mathrm{APEX}} = \arccos\left(\left|\hat{\mathbf{h}} \cdot \hat{\mathbf{p}}\right|\right)$$

where **ĥ** is the unit helix axis vector and **p̂** is the unit pore axis vector. The absolute value of the dot product makes the result independent of helix directionality (N→C vs C→N). A value of 0° indicates that α1 runs parallel to the pore axis (tight funnel); 90° indicates a perpendicular orientation. φ<sub>APEX</sub> was only reported for protomers where the α1-helix linearity (f<sub>PC1</sub>, defined under *Quality indicators* below) exceeded 0.9, to exclude disordered or poorly predicted α1 segments.

**L<sub>APEX</sub>.** The length of the α1-helix, reported as the number of residues in the α1 segment.

**H<sub>ABS</sub>.** Mean Kyte–Doolittle hydrophobicity of the CC-domain helices:

$$H_{\mathrm{ABS}} = \frac{1}{N}\sum_{i=1}^{N} H_{\mathrm{KD}}(a_i)$$

where H<sub>KD</sub>(a<sub>i</sub>) is the Kyte–Doolittle hydrophobicity of amino acid a<sub>i</sub> and *N* is the total number of residues in the first four CC-domain helices. With NB-ARC domain boundaries annotated, helices with start positions before the NB-ARC start coordinate were included.

**μ<sub>H</sub>.** Eisenberg amphipathic moment of the α1 helix:

$$\mu_H = \frac{1}{N}\sqrt{\left(\sum_{i=1}^{N} H_i \sin(i\delta)\right)^2 + \left(\sum_{i=1}^{N} H_i \cos(i\delta)\right)^2}$$

where H<sub>i</sub> is the Eisenberg hydrophobicity of residue *i*, δ = 100° is the angular increment per residue for an α-helix, and *N* is the number of residues in α1. Higher values indicate stronger segregation of hydrophobic residues to one face of the helix, consistent with membrane-inserting amphipathic helices.

## 7. Motif distance

**D<sub>MHD–P</sub>.** Distance between the MHD motif and the P-loop motif within each protomer, computed as the Euclidean distance between the Cα centroids of each motif:

$$D_{\mathrm{MHD\_P}} = \lVert \bar{\mathbf{r}}_{\mathrm{MHD}} - \bar{\mathbf{r}}_{\mathrm{Ploop}} \rVert$$

where $\bar{\mathbf{r}}_{\mathrm{MHD}}$ and $\bar{\mathbf{r}}_{\mathrm{Ploop}}$ are the mean Cα coordinates over the residue ranges of the MHD and P-loop motifs, respectively. Motif positions were provided as pre-computed residue coordinate ranges per protein.

## 8. Quality indicators

**N-terminal correction for spurious gene-model extensions.** Automated gene prediction can introduce spurious residues upstream of the true N-terminal methionine of the MADA motif. With an HMM-based MADA annotation, the pipeline identified the biologically relevant N-terminus as the methionine residue at the start of the MADA HMM alignment, or the closest upstream methionine if the alignment start position was not itself a methionine. This corrected N-terminal position was used as the reference for pore aperture measurements (D<sub>APEX</sub>) and as the start of the α1-helix segment. Proteins without a MADA HMM hit used the first modelled residue.

**α1-helix boundary definition.** The boundary of the first α-helix (α1) in the CC domain was defined in two steps:

1. **HMM-based definition.** `hmmsearch` was run against a bespoke database of chain A sequences using the MADA HMM with an E-value threshold of 0.01; the α1-helix was defined from the corrected N-terminal methionine to the alignment end coordinate (`ali_to`). The threshold was chosen empirically: because the MADA HMM was searched against a small, bespoke database rather than a large reference proteome, the effective E-value scale is set by this small search space, and more permissive thresholds did not recover additional true-positive α1 segments in reference NRC sequences with experimentally resolved structures. In practice, hits above this threshold were absent from our NRC dataset, so the specific cutoff does not affect any of the reported results. Proteins returning no hit were handled as described below.
2. **No α1 assignment.** When the MADA HMM was provided but yielded no hit, α1-dependent metrics (φ<sub>APEX</sub>, L<sub>APEX</sub>, μ<sub>H</sub>) were reported as NA. A DSSP-based fallback is available via the `--hmm_fallback_dssp` flag of `run_pipeline.py` (see `scripts/resistosome_pipeline/README.md`) but was not used for the primary analysis.

**α1-helix linearity (f<sub>PC1</sub>).** The fraction of variance in Cα positions explained by the first principal component, computed per protomer over the α1 segment (equivalent to the quantity sometimes reported as R² for a PCA fit). Concretely, the *k* × 3 matrix of mean-centred Cα coordinates for the *k* residues of the α1 segment was decomposed with `numpy.linalg.svd`, yielding singular values σ<sub>1</sub> ≥ σ<sub>2</sub> ≥ σ<sub>3</sub>, and

$$f_{\mathrm{PC1}} = \sigma_1^2 \big/ (\sigma_1^2 + \sigma_2^2 + \sigma_3^2)$$

A perfectly straight helix yields f<sub>PC1</sub> ≈ 1.0; kinked or disordered segments yield lower values. This metric was used to gate φ<sub>APEX</sub> reporting: protomers with f<sub>PC1</sub> < 0.9 were assigned φ<sub>APEX</sub> = NA to prevent unreliable tilt measurements from poorly resolved α1 segments.

**α1-helix pLDDT.** Mean predicted local distance difference test (pLDDT) score over the α1 region, extracted from the B-factor column of AlphaFold 3 output structures. Values range from 0 to 100, with higher values indicating greater per-residue confidence.

## 9. Cross-replicate aggregation

Each NLR was modelled across *n* = 3 seed replicates to assess prediction robustness. *n* = 3 was chosen as a pragmatic compromise between cost (each replicate requires a full AlphaFold 3 inference run) and coverage of the large NRC dataset, and is insufficient to support formal statistical inference on individual NLRs. The aggregation scheme below therefore uses variance as a penalty on the mean rather than as the basis of a confidence interval.

Metrics were aggregated in a two-level scheme. First, within each structure, per-interface values (contacts, BSA) or per-protomer values (φ<sub>APEX</sub>, L<sub>APEX</sub>, H<sub>ABS</sub>, μ<sub>H</sub>, D<sub>MHD–P</sub>) were summarised. Second, the within-structure summaries were aggregated across replicates using confidence-bound penalties.

For metrics where higher values indicate better structural prediction quality (ipTM, Σ<sub>CONTACTS</sub>, per-type contacts, BSA<sub>INTER-PROTO</sub>), the Lower Confidence Bound (LCB) was applied:

$$\mathrm{LCB} = \bar{x} - k \cdot \frac{s}{\sqrt{n}}$$

For metrics where lower values indicate better structural prediction quality (σ<sub>θ ROT</sub>, S<sub>PROTO</sub>, D<sub>APEX</sub>), the Upper Confidence Bound (UCB) was applied:

$$\mathrm{UCB} = \bar{x} + k \cdot \frac{s}{\sqrt{n}}$$

where x̄ is the mean across *n* replicates (*n* = 3 in this study), *s* is the sample standard deviation (NumPy `std` with `ddof = 1`, i.e. normalised by *n* − 1), and *k* = 1.96 is a fixed penalty weight. With *n* = 3, LCB and UCB are **not** valid 95 % confidence intervals and should not be interpreted as such: the Gaussian approximation requires far larger *n*, and the value of *k* was chosen as a conventional shorthand rather than a statistically justified interval width. Functionally, the penalty down-weights replicate means by a multiple of the per-mean standard error (s/√n), so that NLRs with equivalent means are ordered by the reproducibility of their replicates. A simpler two-key rank (by mean, then by standard deviation) gives broadly similar results; we retain the LCB/UCB form for ease of downstream aggregation into a single scalar score per NLR. LCB/UCB values were computed directly with `numpy.mean` and `numpy.std(ddof=1)` rather than by calling a dedicated confidence-interval library. This penalisation scheme ensures that NLRs with reproducible predictions across seeds are ranked above those with equivalent mean scores but high cross-replicate variance.

Non-penalised metrics (φ<sub>APEX</sub>, L<sub>APEX</sub>, H<sub>ABS</sub>, μ<sub>H</sub>, D<sub>MHD–P</sub>) were reported as mean ± standard deviation across replicates.

## 10. Software and dependencies

Analyses were run under Python 3.12 with BioPython 1.84, SciPy 1.15.2, NumPy 1.26.3, pandas 2.3.3, openpyxl 3.1.5, FreeSASA 2.2.1, HMMER 3.4 (August 2023 release) and DSSP 4.5.8. Minimum dependency versions are defined in `scripts/resistosome_pipeline/requirements.txt`.

Key references for the external tools used:

- Shrake & Rupley, *J. Mol. Biol.* **79**, 351–371 (1973) — SASA algorithm.
- Mitternacht, *F1000Research* (2016), <https://doi.org/10.12688/f1000research.7931.1> — FreeSASA.
- Kabsch & Sander, *Biopolymers* **22**, 2577–2637 (1983) — DSSP.
- Eddy, *PLOS Comput. Biol.* **7**, e1002195 (2011) — HMMER.
- Cock et al., *Bioinformatics* **25**, 1422–1423 (2009) — BioPython.
- Virtanen et al., *Nat. Methods* **17**, 261–272 (2020) — SciPy.
- Harris et al., *Nature* **585**, 357–362 (2020) — NumPy.
- Kyte & Doolittle, *J. Mol. Biol.* **157**, 105–132 (1982) — Kyte–Doolittle hydrophobicity.
- Eisenberg, Weiss & Terwilliger, *Proc. Natl. Acad. Sci.* **81**, 140–144 (1984) — hydrophobic moment.
- Abramson et al., *Nature* (2024) — AlphaFold 3.
- Adachi et al., *eLife* **8** (2019) — MADA motif.
