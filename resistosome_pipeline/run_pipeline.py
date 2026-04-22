#!/usr/bin/env python3
"""
Resistosome Structure Analysis Pipeline
========================================

Batch analysis of predicted (AlphaFold 3) or empirical resistosome structures.

Computes per-structure and cross-replicate metrics including:
  - ipTM (AF3), contacts, BSA, interface angles, symmetry RMSD
  - N-terminal distances, helix geometry, hydrophobicity, amphipathic moment
  - MHD-to-P-loop distance (with annotation files)

Outputs:
  - Per-structure JSON files with detailed metrics
  - Summary Excel/CSV table with penalized (LCB/UCB) scores

Usage:
    python run_pipeline.py \\
        --input_dir /path/to/structures \\
        --output_dir /path/to/output \\
        --mhd_csv /path/to/NRCH_MHD_nlrexp_coords.csv \\
        --ploop_csv /path/to/NRCH_ploop_nlrexp_coords.csv \\
        [--contact_cutoff 3.5] \\
        [--pae_cutoff 10.0] \\
        [--k 1.96] \\
        [--log_level INFO]
"""

import argparse
import logging
import sys
import os


def main():
    parser = argparse.ArgumentParser(
        description='Resistosome Structure Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--input_dir', '-i', required=True,
                        help='Directory containing structure files (PDB/CIF/tar.gz)')
    parser.add_argument('--output_dir', '-o', required=True,
                        help='Output directory for results')
    parser.add_argument('--mhd_csv', default=None,
                        help='Path to MHD motif annotation CSV')
    parser.add_argument('--ploop_csv', default=None,
                        help='Path to P-loop motif annotation CSV')
    parser.add_argument('--nbarc_csv', default=None,
                        help='Path to NB-ARC domain annotation CSV (defines CC/NB-ARC boundary)')
    parser.add_argument('--protein_names', default=None,
                        help='Path to text file with protein names (one per line). '
                             'Used to match filenames to protein IDs instead of regex parsing.')
    parser.add_argument('--hmm', default=None,
                        help='Path to HMM model for α1-helix boundary detection. '
                             'Extracts sequences from structures, runs hmmsearch, '
                             'and uses hit coordinates to define α1. Requires HMMER3.')
    parser.add_argument('--hmm_domtbl', default=None,
                        help='Path to pre-computed hmmsearch --domtblout file. '
                             'Alternative to --hmm when you already have the results.')
    parser.add_argument('--hmm_fallback_dssp', action='store_true', default=False,
                        help='When --hmm is provided, fall back to DSSP for α1 if '
                             'no HMM hit is found. Default: off (report NA).')
    parser.add_argument('--workers', type=int, default=None,
                        help='Number of parallel workers for structure analysis. '
                             'Default: min(CPU count, 8). Use 1 for serial mode.')
    parser.add_argument('--contact_cutoff', type=float, default=3.5,
                        help='Contact distance cutoff in Å (default: 3.5)')
    parser.add_argument('--pae_cutoff', type=float, default=10.0,
                        help='PAE threshold for contact filtering (default: 10.0)')
    parser.add_argument('--k', type=float, default=1.96,
                        help='Confidence bound k-value (default: 1.96 for 95%% CI)')
    parser.add_argument('--work_dir', default=None,
                        help='Temporary working directory (default: output_dir/tmp)')
    parser.add_argument('--log_level', default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        help='Logging level (default: INFO)')

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s [%(levelname)s] %(name)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(os.path.join(args.output_dir, 'pipeline.log')
                                if os.path.isdir(args.output_dir)
                                else 'pipeline.log'),
        ]
    )

    # Override config with CLI arguments
    from resistosome import config
    config.CONTACT_CUTOFF_A = args.contact_cutoff
    config.PAE_CUTOFF_A = args.pae_cutoff
    config.K_CONFIDENCE = args.k

    # Ensure output dir exists (needed for log file)
    os.makedirs(args.output_dir, exist_ok=True)

    from resistosome.pipeline import run_pipeline

    logger = logging.getLogger(__name__)
    logger.info("=" * 70)
    logger.info("Resistosome Structure Analysis Pipeline")
    logger.info("=" * 70)
    logger.info(f"Input directory : {args.input_dir}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"MHD annotations : {args.mhd_csv or 'None'}")
    logger.info(f"P-loop annot.   : {args.ploop_csv or 'None'}")
    logger.info(f"NB-ARC annot.   : {args.nbarc_csv or 'None'}")
    logger.info(f"Protein names   : {args.protein_names or 'None (regex mode)'}")
    logger.info(f"HMM model       : {args.hmm or 'None'}")
    logger.info(f"HMM domtbl      : {args.hmm_domtbl or 'None'}")
    logger.info(f"Workers         : {args.workers or 'auto'}")
    logger.info(f"Contact cutoff  : {args.contact_cutoff} Å")
    logger.info(f"PAE cutoff      : {args.pae_cutoff} Å")
    logger.info(f"k (CI)          : {args.k}")
    logger.info("=" * 70)

    summary_path = run_pipeline(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        mhd_csv=args.mhd_csv,
        ploop_csv=args.ploop_csv,
        nbarc_csv=args.nbarc_csv,
        protein_names_file=args.protein_names,
        hmm_path=args.hmm,
        hmm_domtbl=args.hmm_domtbl,
        hmm_fallback_dssp=args.hmm_fallback_dssp,
        workers=args.workers,
        work_dir=args.work_dir,
    )

    if summary_path:
        logger.info(f"\n✓ Pipeline complete. Summary: {summary_path}")
    else:
        logger.warning("\n✗ Pipeline finished with no output.")


if __name__ == '__main__':
    main()
