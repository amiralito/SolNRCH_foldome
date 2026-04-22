#!/usr/bin/env python3
"""
Extract confidence scores from AlphaFold 3 summary_confidences.json files.

Parses the summary_confidences.json output and extracts key metrics:
- pTM (predicted TM-score)
- ipTM (interface pTM)
- ranking_score
- fraction_disordered
- has_clash
- chain_ptm (per-chain pTM)
- chain_iptm (per-chain interface pTM)
- chain_pair_iptm (pairwise interface pTM)
- chain_pair_pae_min (pairwise PAE minimum)

Outputs two TSV files:
1. Overview scores (scalar metrics + summary stats)
2. Per-chain pairwise scores as matrices
"""

import json
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Any, Optional
import re


def load_json(filepath: Path) -> Dict[str, Any]:
    """Load JSON file and return parsed data."""
    with open(filepath, 'r') as f:
        return json.load(f)


def extract_model_name(filepath: Path) -> str:
    """
    Extract model name from summary_confidences.json path.
    
    Expected format: {model_name}_summary_confidences.json
    """
    name = filepath.stem
    # Remove _summary_confidences suffix
    name = re.sub(r'_summary_confidences$', '', name)
    return name


def extract_seed_from_name(name: str) -> Optional[int]:
    """Extract seed number from model name if present."""
    match = re.search(r'seed(\d+)', name, re.IGNORECASE)
    if match:
        return int(match.group(1))
    return None


def extract_scores(filepath: Path) -> Dict[str, Any]:
    """
    Extract all relevant scores from a summary_confidences.json file.
    
    Returns a dictionary with all extracted metrics.
    """
    data = load_json(filepath)
    model_name = extract_model_name(filepath)
    
    # Initialize scores dictionary
    scores = {
        'model_name': model_name,
        'seed': extract_seed_from_name(model_name),
        'ptm': None,
        'iptm': None,
        'ranking_score': None,
        'fraction_disordered': None,
        'has_clash': None,
        # Raw data for matrix output
        '_chain_ids': [],
        '_chain_ptm': [],
        '_chain_iptm': [],
        '_chain_pair_iptm': [],
        '_chain_pair_pae_min': [],
    }
    
    # Extract basic scalar scores
    for key in ['ptm', 'iptm', 'ranking_score', 'fraction_disordered', 'has_clash']:
        if key in data:
            scores[key] = data[key]
    
    # Extract chain IDs if available
    if 'chain_ids' in data:
        scores['_chain_ids'] = data['chain_ids']
    
    # Extract chain-level scores
    if 'chain_ptm' in data:
        scores['_chain_ptm'] = data['chain_ptm']
    
    if 'chain_iptm' in data:
        scores['_chain_iptm'] = data['chain_iptm']
    
    # Extract pairwise matrices
    if 'chain_pair_iptm' in data:
        scores['_chain_pair_iptm'] = data['chain_pair_iptm']
    
    if 'chain_pair_pae_min' in data:
        scores['_chain_pair_pae_min'] = data['chain_pair_pae_min']
    
    # Calculate summary statistics for pair ipTM
    if 'chain_pair_iptm' in data:
        pair_vals = []
        chain_pair_iptm = data['chain_pair_iptm']
        if isinstance(chain_pair_iptm, list):
            for i, row in enumerate(chain_pair_iptm):
                if isinstance(row, list):
                    for j, val in enumerate(row):
                        if i < j and val is not None:  # Upper triangle only
                            pair_vals.append(val)
        if pair_vals:
            scores['mean_pair_iptm'] = sum(pair_vals) / len(pair_vals)
            scores['min_pair_iptm'] = min(pair_vals)
            scores['max_pair_iptm'] = max(pair_vals)
    
    return scores


def find_summary_confidence_files(directory: Path, recursive: bool = True) -> List[Path]:
    """Find all summary_confidences.json files in a directory."""
    if recursive:
        return sorted(directory.rglob('*_summary_confidences.json'))
    else:
        return sorted(directory.glob('*_summary_confidences.json'))


def write_overview_tsv(scores_list: List[Dict[str, Any]], output_file: Path) -> None:
    """
    Write overview scores to a TSV file.
    
    Contains scalar metrics and summary statistics only.
    """
    if not scores_list:
        return
    
    # Fixed columns for overview
    columns = [
        'model_name', 'seed', 'ptm', 'iptm', 'ranking_score',
        'fraction_disordered', 'has_clash', 'mean_pair_iptm',
        'min_pair_iptm', 'max_pair_iptm'
    ]
    
    # Filter to only columns that exist in data
    all_keys = set()
    for scores in scores_list:
        all_keys.update(scores.keys())
    columns = [c for c in columns if c in all_keys]
    
    # Ensure output directory exists
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join(columns) + '\n')
        
        # Write data rows
        for scores in scores_list:
            row = []
            for col in columns:
                val = scores.get(col, '')
                if val is None:
                    row.append('')
                elif isinstance(val, float):
                    row.append(f'{val:.6f}')
                elif isinstance(val, bool):
                    row.append(str(val).lower())
                else:
                    row.append(str(val))
            f.write('\t'.join(row) + '\n')


def write_matrix_tsv(scores: Dict[str, Any], output_file: Path) -> None:
    """
    Write per-chain pairwise scores as matrices to a TSV file.
    
    Format:
    ## chain_pair_iptm
    	A	B	C	D	E	F
    A	1.0	0.8	0.7	...
    B	0.8	1.0	0.75	...
    ...
    
    ## chain_pair_pae_min
    ...
    """
    chain_ids = scores.get('_chain_ids', [])
    chain_pair_iptm = scores.get('_chain_pair_iptm', [])
    chain_pair_pae_min = scores.get('_chain_pair_pae_min', [])
    chain_ptm = scores.get('_chain_ptm', [])
    chain_iptm = scores.get('_chain_iptm', [])
    
    # If no chain IDs provided, generate them
    if not chain_ids and chain_pair_iptm:
        n_chains = len(chain_pair_iptm)
        chain_ids = [chr(65 + i) if i < 26 else f"C{i}" for i in range(n_chains)]
    
    if not chain_ids:
        return
    
    # Ensure output directory exists
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        # Write model name header
        f.write(f"# Model: {scores.get('model_name', 'unknown')}\n")
        f.write(f"# Seed: {scores.get('seed', 'N/A')}\n")
        f.write("\n")
        
        # Write per-chain scores (single row for each metric)
        if chain_ptm:
            f.write("## chain_ptm\n")
            f.write('\t'.join(['chain'] + chain_ids) + '\n')
            row = ['ptm'] + [f'{v:.4f}' if v is not None else '' for v in chain_ptm]
            f.write('\t'.join(row) + '\n')
            f.write("\n")
        
        if chain_iptm:
            f.write("## chain_iptm\n")
            f.write('\t'.join(['chain'] + chain_ids) + '\n')
            row = ['iptm'] + [f'{v:.4f}' if v is not None else '' for v in chain_iptm]
            f.write('\t'.join(row) + '\n')
            f.write("\n")
        
        # Write chain_pair_iptm matrix
        if chain_pair_iptm:
            f.write("## chain_pair_iptm\n")
            # Header row
            f.write('\t'.join([''] + chain_ids) + '\n')
            # Data rows
            for i, row_data in enumerate(chain_pair_iptm):
                if i < len(chain_ids):
                    row_label = chain_ids[i]
                    if isinstance(row_data, list):
                        row_vals = [f'{v:.4f}' if v is not None else '' for v in row_data]
                    else:
                        row_vals = [str(row_data) if row_data is not None else '']
                    f.write('\t'.join([row_label] + row_vals) + '\n')
            f.write("\n")
        
        # Write chain_pair_pae_min matrix
        if chain_pair_pae_min:
            f.write("## chain_pair_pae_min\n")
            # Header row
            f.write('\t'.join([''] + chain_ids) + '\n')
            # Data rows
            for i, row_data in enumerate(chain_pair_pae_min):
                if i < len(chain_ids):
                    row_label = chain_ids[i]
                    if isinstance(row_data, list):
                        row_vals = [f'{v:.4f}' if v is not None else '' for v in row_data]
                    else:
                        row_vals = [str(row_data) if row_data is not None else '']
                    f.write('\t'.join([row_label] + row_vals) + '\n')


def process_single_file(
    filepath: Path,
    output_file: Optional[Path] = None,
    model_name: Optional[str] = None,
    quiet: bool = False
) -> Dict[str, Any]:
    """
    Process a single summary_confidences.json file.
    
    Returns extracted scores and writes to two output files:
    - {output}_scores.tsv (overview)
    - {output}_matrix.tsv (per-chain matrices)
    """
    if not filepath.exists():
        print(f"Error: File not found: {filepath}", file=sys.stderr)
        sys.exit(1)
    
    scores = extract_scores(filepath)
    
    if model_name:
        scores['model_name'] = model_name
    
    if not quiet:
        print(f"Processed: {filepath.name}")
        print(f"  Model: {scores['model_name']}")
        print(f"  pTM: {scores.get('ptm', 'N/A')}")
        print(f"  ipTM: {scores.get('iptm', 'N/A')}")
        print(f"  Ranking score: {scores.get('ranking_score', 'N/A')}")
    
    if output_file:
        # Determine output file paths
        output_stem = output_file.stem
        # Remove common suffixes to get base name
        for suffix in ['_score_tbl', '_scores', '_matrix']:
            if output_stem.endswith(suffix):
                output_stem = output_stem[:-len(suffix)]
                break
        
        overview_file = output_file.parent / f"{output_stem}_scores.tsv"
        matrix_file = output_file.parent / f"{output_stem}_matrix.tsv"
        
        write_overview_tsv([scores], overview_file)
        write_matrix_tsv(scores, matrix_file)
        
        if not quiet:
            print(f"  Overview: {overview_file}")
            print(f"  Matrix:   {matrix_file}")
    
    return scores


def process_directory(
    directory: Path,
    output_file: Path,
    recursive: bool = True,
    quiet: bool = False
) -> List[Dict[str, Any]]:
    """
    Process all summary_confidences.json files in a directory.
    
    Writes combined overview TSV and individual matrix files.
    """
    files = find_summary_confidence_files(directory, recursive=recursive)
    
    if not files:
        print(f"No summary_confidences.json files found in {directory}", file=sys.stderr)
        return []
    
    if not quiet:
        print(f"Found {len(files)} summary_confidences.json files")
    
    all_scores = []
    for filepath in files:
        try:
            scores = extract_scores(filepath)
            all_scores.append(scores)
            
            # Write individual matrix file
            model_name = scores['model_name']
            matrix_file = output_file.parent / f"{model_name}_matrix.tsv"
            write_matrix_tsv(scores, matrix_file)
            
            if not quiet:
                print(f"  Processed: {model_name}")
        except Exception as e:
            print(f"  Warning: Failed to process {filepath}: {e}", file=sys.stderr)
    
    if all_scores:
        # Write combined overview
        write_overview_tsv(all_scores, output_file)
        if not quiet:
            print(f"\nWrote {len(all_scores)} records to {output_file}")
    
    return all_scores


def write_scores_json(scores: Dict[str, Any], output_file: Path) -> None:
    """Write scores to a JSON file."""
    # Remove internal keys
    output_scores = {k: v for k, v in scores.items() if not k.startswith('_')}
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(output_scores, f, indent=2)


def main():
    """Main function to handle command line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract confidence scores from AlphaFold 3 summary_confidences.json files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process a single file (outputs two files: *_scores.tsv and *_matrix.tsv)
  python extract_af3_scores.py model_summary_confidences.json -o model
  
  # Process all files in a directory
  python extract_af3_scores.py ./outputs -o all_scores.tsv
  
  # Process recursively
  python extract_af3_scores.py ./outputs -o ./scores/all_scores.tsv --recursive

Output files:
  For single file mode:
    - {output}_scores.tsv  : Overview metrics (ptm, iptm, ranking_score, etc.)
    - {output}_matrix.tsv  : Per-chain pairwise scores as matrices
  
  For directory mode:
    - {output}             : Combined overview for all models
    - {model}_matrix.tsv   : Individual matrix file per model
        """
    )
    
    parser.add_argument(
        'input',
        type=Path,
        help='Input file or directory containing summary_confidences.json files'
    )
    parser.add_argument(
        '--output', '-o',
        type=Path,
        required=True,
        help='Output file path (base name for single file, TSV for directory)'
    )
    parser.add_argument(
        '--recursive', '-r',
        action='store_true',
        help='Search recursively in directories'
    )
    parser.add_argument(
        '--json',
        action='store_true',
        help='Output overview as JSON instead of TSV'
    )
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Minimize output'
    )
    parser.add_argument(
        '--model-name',
        type=str,
        default=None,
        help='Override model name (for single file mode)'
    )
    
    args = parser.parse_args()
    
    if not args.input.exists():
        print(f"Error: Input not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    if args.input.is_file():
        # Single file mode
        scores = process_single_file(
            args.input,
            output_file=args.output,
            model_name=args.model_name,
            quiet=args.quiet
        )
        
        if args.json:
            json_output = args.output.parent / f"{args.output.stem}_scores.json"
            write_scores_json(scores, json_output)
        
        if args.quiet:
            # Output just the key scores as JSON for parsing
            print(json.dumps({
                'model_name': scores['model_name'],
                'ptm': scores.get('ptm'),
                'iptm': scores.get('iptm'),
                'ranking_score': scores.get('ranking_score')
            }))
    else:
        # Directory mode
        all_scores = process_directory(
            args.input,
            args.output,
            recursive=args.recursive,
            quiet=args.quiet
        )
        
        if args.quiet and all_scores:
            print(json.dumps({'processed': len(all_scores), 'output': str(args.output)}))


if __name__ == "__main__":
    main()
