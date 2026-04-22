#!/usr/bin/env python3
"""
Generate AlphaFold 3 input JSON files for homomeric modeling.

This script takes precomputed MSA/template JSON files and generates inputs
with specified numbers of protomers and optional ligands.

Features:
- Single input directory (homomeric modeling)
- Variable number of protomers (multiple values → separate jobs)
- Ligand support (CCD codes or SMILES strings)
- Variable number of ligands per job
- Job naming: {input_name}_{X}protomer_{Y}Ligand_seed{N}

Supports:
- Single file processing
- Batch processing of all files in a directory
- Single-file mode for array job integration (--index)
- Manifest generation (--list-combinations)
"""

import json
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from itertools import product
import string
import re


def load_json(filepath: Path) -> Dict[str, Any]:
    """Load JSON file and return parsed data."""
    with open(filepath, 'r') as f:
        return json.load(f)


def save_json(data: Dict[str, Any], filepath: Path, indent: int = 2) -> None:
    """Save data to JSON file with proper formatting."""
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=indent)


def get_json_files(path: Path) -> List[Path]:
    """Get all JSON files from a path (file or directory), recursively."""
    if path.is_file():
        if path.suffix.lower() == '.json':
            return [path]
        else:
            print(f"Warning: {path} is not a JSON file, skipping.", file=sys.stderr)
            return []
    elif path.is_dir():
        json_files = sorted(path.rglob('*.json'))
        return json_files
    else:
        print(f"Error: {path} is neither a file nor a directory", file=sys.stderr)
        return []


def extract_protein_name(filepath: Path) -> str:
    """Extract clean protein name from filepath, removing common suffixes."""
    name = filepath.stem
    # Remove common suffixes
    for suffix in ['_data', '_msa', '_template', '_input', '_precomputed']:
        name = name.replace(suffix, '')
    return name


def parse_ligand_input(ligand_input: str) -> Dict[str, Any]:
    """
    Parse ligand input string or file path.
    
    Returns a dictionary with ligand info for AF3 JSON (alphafold3 dialect):
    - CCD codes: {"type": "ccd", "value": "OLA"} -> {"ligand": {"id": "X", "ccdCodes": ["OLA"]}}
    - SMILES strings: {"type": "smiles", "value": "..."} -> {"ligand": {"id": "X", "smiles": "..."}}
    
    Input formats:
    - CCD code: "CCD_OLA", "OLA", "ccd:OLA" -> stored as "OLA"
    - SMILES: "smiles:CCO" or auto-detected from special chars
    - File path: @/path/to/ligand.txt (reads first line)
    """
    # Check if it's a file reference
    if ligand_input.startswith('@'):
        filepath = Path(ligand_input[1:])
        if not filepath.exists():
            print(f"Error: Ligand file not found: {filepath}", file=sys.stderr)
            sys.exit(1)
        with open(filepath, 'r') as f:
            ligand_input = f.readline().strip()
    
    # Parse the ligand specification
    ligand_input = ligand_input.strip()
    
    # Strip CCD_ prefix if present
    if ligand_input.upper().startswith('CCD_'):
        ccd_code = ligand_input[4:].strip().upper()
        return {'type': 'ccd', 'value': ccd_code}
    
    # Explicit ccd: prefix
    if ligand_input.lower().startswith('ccd:'):
        ccd_code = ligand_input[4:].strip().upper()
        return {'type': 'ccd', 'value': ccd_code}
    
    # Explicit smiles: prefix
    if ligand_input.lower().startswith('smiles:'):
        smiles = ligand_input[7:].strip()
        return {'type': 'smiles', 'value': smiles}
    
    # Auto-detect based on content
    # SMILES contain special characters like = ( ) # @ [ ] + - / \ .
    smiles_chars = set('=()#@[]+-/\\.')
    
    if any(c in ligand_input for c in smiles_chars):
        # Contains SMILES-specific characters
        return {'type': 'smiles', 'value': ligand_input}
    elif re.match(r'^[A-Za-z0-9]{1,10}$', ligand_input):
        # Pure alphanumeric, 1-10 chars -> treat as CCD code
        return {'type': 'ccd', 'value': ligand_input.upper()}
    else:
        # Default to SMILES for anything else
        return {'type': 'smiles', 'value': ligand_input}


def create_ligand_entry(ligand_info: Dict[str, str], chain_id: str) -> Dict[str, Any]:
    """
    Create a ligand entry for the AF3 sequences array (alphafold3 dialect).
    
    Format: {"ligand": {"id": "F", "ccdCodes": ["OLA"]}}
    or:     {"ligand": {"id": "F", "smiles": "CCO"}}
    """
    if ligand_info['type'] == 'ccd':
        return {
            'ligand': {
                'id': chain_id,
                'ccdCodes': [ligand_info['value']]
            }
        }
    else:  # smiles
        return {
            'ligand': {
                'id': chain_id,
                'smiles': ligand_info['value']
            }
        }


def get_ligand_name(ligand_info: Dict[str, str]) -> str:
    """Get a short name for the ligand for job naming."""
    if ligand_info['type'] == 'ccd':
        # Value is like "OLA" - use directly
        return ligand_info['value']
    else:
        # Use first 8 chars of SMILES hash for naming
        import hashlib
        smiles_hash = hashlib.md5(ligand_info['value'].encode()).hexdigest()[:8]
        return f"smi{smiles_hash}"


def generate_output_name(
    protein_name: str,
    num_protomers: int,
    ligand_info: Optional[Dict[str, str]] = None,
    num_ligands: int = 0,
    seed: Optional[int] = None,
    job_number: Optional[int] = None
) -> str:
    """
    Generate output name following the convention:
    job_{N}_{input_name}_{X}protomers_{Y}{Ligand}Ligs_seed{S}
    
    Args:
        protein_name: Base protein name
        num_protomers: Number of protomers
        ligand_info: Ligand information (optional)
        num_ligands: Number of ligands (0 if no ligand)
        seed: Seed number (optional)
        job_number: Job number for prefixing (optional)
    """
    parts = []
    
    # Add job number prefix if provided
    if job_number is not None:
        parts.append(f"job_{job_number}")
    
    # Protein name
    parts.append(protein_name)
    
    # Protomer count
    protomer_suffix = "protomer" if num_protomers == 1 else "protomers"
    parts.append(f"{num_protomers}{protomer_suffix}")
    
    # Ligand info
    if ligand_info and num_ligands > 0:
        ligand_name = get_ligand_name(ligand_info)
        ligand_suffix = "Lig" if num_ligands == 1 else "Ligs"
        parts.append(f"{num_ligands}{ligand_name}{ligand_suffix}")
    
    # Seed
    if seed is not None:
        parts.append(f"seed{seed}")
    
    return "_".join(parts)


def generate_homomeric_input(
    template_data: Dict[str, Any],
    template_path: Path,
    num_protomers: int,
    ligand_info: Optional[Dict[str, str]] = None,
    num_ligands: int = 0,
    custom_seeds: Optional[List[int]] = None,
    output_name: Optional[str] = None,
    verbose: bool = True
) -> Dict[str, Any]:
    """
    Generate AF3 input JSON for homomeric modeling (alphafold3 dialect).
    
    Uses the alphafold3 format with individual entries for each protomer
    and ligand copy, each with a unique chain ID (A, B, C, ...).
    
    Args:
        template_data: Source JSON data with MSA/templates
        template_path: Path to source file (for naming)
        num_protomers: Number of protein copies
        ligand_info: Optional ligand information
        num_ligands: Number of ligand copies
        custom_seeds: Optional seed numbers
        output_name: Optional name override
        verbose: Print progress messages
    
    Returns:
        Complete AF3 input dictionary in alphafold3 dialect
    """
    # Generate chain IDs - extend beyond 26 if needed
    # A-Z (26), then AA, BA, CA, ..., ZA, AB, BB, ..., ZZ
    def generate_chain_ids(n: int) -> List[str]:
        """Generate n unique chain IDs."""
        ids = []
        # First pass: A-Z
        for c in string.ascii_uppercase:
            ids.append(c)
            if len(ids) >= n:
                return ids
        # Extended: AA, BA, CA, ..., ZA, AB, BB, ..., ZZ
        for c2 in string.ascii_uppercase:
            for c1 in string.ascii_uppercase:
                ids.append(f"{c1}{c2}")
                if len(ids) >= n:
                    return ids
        return ids
    
    total_chains = num_protomers + (num_ligands if ligand_info else 0)
    all_chain_ids = generate_chain_ids(total_chains)
    
    # Get the protein sequence from template - handle both formats
    protein_seq_entry = None
    protein_sequence_str = None
    
    for seq in template_data.get('sequences', []):
        if 'protein' in seq:
            protein_seq_entry = seq['protein']
            protein_sequence_str = protein_seq_entry.get('sequence', '')
            break
        elif 'proteinChain' in seq:
            # Convert from alphafoldserver format
            protein_sequence_str = seq['proteinChain'].get('sequence', '')
            protein_seq_entry = {
                'sequence': protein_sequence_str
            }
            # Copy MSA/template info if present
            if 'unpairedMsa' in seq['proteinChain']:
                protein_seq_entry['unpairedMsa'] = seq['proteinChain']['unpairedMsa']
            if 'pairedMsa' in seq['proteinChain']:
                protein_seq_entry['pairedMsa'] = seq['proteinChain']['pairedMsa']
            if 'templates' in seq['proteinChain']:
                protein_seq_entry['templates'] = seq['proteinChain']['templates']
            break
    
    if protein_seq_entry is None or not protein_sequence_str:
        raise ValueError(f"No protein sequence found in {template_path}")
    
    # Build output data structure using alphafold3 dialect
    output_data = {
        'name': '',
        'modelSeeds': [],
        'sequences': [],
        'dialect': 'alphafold3',
        'version': template_data.get('version', 1)
    }
    
    # Set name
    if output_name:
        output_data['name'] = output_name
    else:
        protein_name = extract_protein_name(template_path)
        output_data['name'] = generate_output_name(
            protein_name, num_protomers, ligand_info, num_ligands
        )
    
    if verbose:
        print(f"  Generating: {output_data['name']}")
        print(f"    Protomers: {num_protomers}")
    
    # Assign chain IDs
    protein_chain_ids = all_chain_ids[:num_protomers]
    ligand_chain_ids = all_chain_ids[num_protomers:num_protomers + num_ligands] if num_ligands > 0 else []
    
    # Create protein entry with list of IDs for homomeric copies
    protein_entry = {
        'protein': {
            'id': protein_chain_ids if num_protomers > 1 else protein_chain_ids[0],
            'sequence': protein_sequence_str
        }
    }
    
    # Copy MSA and template information from source
    if 'unpairedMsa' in protein_seq_entry:
        protein_entry['protein']['unpairedMsa'] = protein_seq_entry['unpairedMsa']
    if 'unpairedMsaPath' in protein_seq_entry:
        protein_entry['protein']['unpairedMsaPath'] = protein_seq_entry['unpairedMsaPath']
    if 'pairedMsa' in protein_seq_entry:
        protein_entry['protein']['pairedMsa'] = protein_seq_entry['pairedMsa']
    if 'pairedMsaPath' in protein_seq_entry:
        protein_entry['protein']['pairedMsaPath'] = protein_seq_entry['pairedMsaPath']
    if 'templates' in protein_seq_entry:
        protein_entry['protein']['templates'] = protein_seq_entry['templates']
    if 'modifications' in protein_seq_entry:
        protein_entry['protein']['modifications'] = protein_seq_entry['modifications']
    
    output_data['sequences'].append(protein_entry)
    
    if verbose:
        first_id = protein_chain_ids[0]
        last_id = protein_chain_ids[-1] if num_protomers > 1 else first_id
        print(f"    Protein chains: {first_id}-{last_id} ({len(protein_sequence_str)} residues each)")
    
    # Add ligand with list of IDs if specified
    if ligand_info and num_ligands > 0:
        ligand_entry = {
            'ligand': {
                'id': ligand_chain_ids if num_ligands > 1 else ligand_chain_ids[0]
            }
        }
        
        if ligand_info['type'] == 'ccd':
            ligand_entry['ligand']['ccdCodes'] = [ligand_info['value']]
        else:  # smiles
            ligand_entry['ligand']['smiles'] = ligand_info['value']
        
        output_data['sequences'].append(ligand_entry)
        
        if verbose:
            ligand_name = get_ligand_name(ligand_info)
            first_lig_id = ligand_chain_ids[0]
            last_lig_id = ligand_chain_ids[-1] if num_ligands > 1 else first_lig_id
            print(f"    Ligand chains: {first_lig_id}-{last_lig_id} ({ligand_name}, {num_ligands} copies)")
    
    # Handle model seeds
    if custom_seeds is not None:
        output_data['modelSeeds'] = custom_seeds
        if verbose:
            print(f"    Seeds: {custom_seeds}")
    elif 'modelSeeds' in template_data:
        output_data['modelSeeds'] = template_data['modelSeeds']
    else:
        # Default seed if none specified
        output_data['modelSeeds'] = [1]
    
    # Copy bonded atom pairs if present
    if 'bondedAtomPairs' in template_data:
        output_data['bondedAtomPairs'] = template_data['bondedAtomPairs']
    
    # Copy user CCD if present
    if 'userCCD' in template_data:
        output_data['userCCD'] = template_data['userCCD']
    if 'userCCDPath' in template_data:
        output_data['userCCDPath'] = template_data['userCCDPath']
    
    return output_data


def generate_single_job(
    input_dir: Path,
    index: int,
    num_protomers: int,
    output_file: Path,
    ligand_info: Optional[Dict[str, str]] = None,
    num_ligands: int = 0,
    custom_seeds: Optional[List[int]] = None,
    job_number: Optional[int] = None,
    quiet: bool = False
) -> Dict[str, Any]:
    """
    Generate a single job based on index.
    
    For SLURM array job integration.
    
    Returns:
        Dictionary with metadata about what was generated.
    """
    files = get_json_files(input_dir)
    
    if not files:
        print(f"Error: No JSON files found in {input_dir}", file=sys.stderr)
        sys.exit(1)
    
    if index < 0 or index >= len(files):
        print(f"Error: index={index} out of range [0, {len(files)-1}]", file=sys.stderr)
        sys.exit(1)
    
    # Select the file
    template_file = files[index]
    protein_name = extract_protein_name(template_file)
    
    if not quiet:
        print(f"Single-job mode:")
        print(f"  Input [{index}]: {template_file.name} -> {protein_name}")
    
    # Load template
    template_data = load_json(template_file)
    
    # Generate output name
    output_name = generate_output_name(
        protein_name, num_protomers, ligand_info, num_ligands,
        seed=custom_seeds[0] if custom_seeds and len(custom_seeds) == 1 else None,
        job_number=job_number
    )
    
    # Generate the input JSON
    result_data = generate_homomeric_input(
        template_data, template_file,
        num_protomers=num_protomers,
        ligand_info=ligand_info,
        num_ligands=num_ligands,
        custom_seeds=custom_seeds,
        output_name=output_name,
        verbose=not quiet
    )
    
    # Ensure output directory exists
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Save the JSON
    save_json(result_data, output_file)
    
    if not quiet:
        print(f"  Output: {output_file}")
    
    # Return metadata
    return {
        'protein': protein_name,
        'num_protomers': num_protomers,
        'ligand': get_ligand_name(ligand_info) if ligand_info else None,
        'num_ligands': num_ligands,
        'output_file': str(output_file),
        'output_name': output_name
    }


def list_combinations(
    input_dir: Path,
    protomer_counts: List[int],
    ligand_info: Optional[Dict[str, str]] = None,
    ligand_counts: List[int] = None,
    custom_seeds: Optional[List[int]] = None,
    expand_seeds: bool = False,
    batch_name: Optional[str] = None
) -> List[Dict[str, Any]]:
    """
    List all combinations that would be generated.
    
    Args:
        input_dir: Directory containing input JSON files
        protomer_counts: List of protomer numbers to generate
        ligand_info: Optional ligand information
        ligand_counts: List of ligand counts (default: [0] or [1] if ligand specified)
        custom_seeds: Optional seed numbers
        expand_seeds: Create separate entries per seed
        batch_name: Optional batch name for manifest
    
    Returns:
        List of dictionaries with combination metadata
    """
    files = get_json_files(input_dir)
    
    if not files:
        return []
    
    # Set defaults
    if ligand_counts is None:
        ligand_counts = [1] if ligand_info else [0]
    
    # Filter out 0 ligands if no ligand specified
    if not ligand_info:
        ligand_counts = [0]
    
    # Determine seeds to iterate over
    if expand_seeds and custom_seeds and len(custom_seeds) > 1:
        seeds_to_expand = custom_seeds
    else:
        seeds_to_expand = [None]  # None means no seed suffix
    
    combinations = []
    job_number = 0
    
    for file_idx, template_file in enumerate(files):
        protein_name = extract_protein_name(template_file)
        
        for num_protomers in protomer_counts:
            for num_ligands in ligand_counts:
                # Skip invalid combinations
                if num_ligands > 0 and not ligand_info:
                    continue
                
                for seed in seeds_to_expand:
                    output_name = generate_output_name(
                        protein_name, num_protomers, ligand_info, num_ligands,
                        seed=seed, job_number=job_number
                    )
                    
                    combo = {
                        'batch': batch_name or '',
                        'job_number': job_number,
                        'file_index': file_idx,
                        'protein': protein_name,
                        'file_path': str(template_file),
                        'num_protomers': num_protomers,
                        'ligand': get_ligand_name(ligand_info) if ligand_info else '',
                        'num_ligands': num_ligands,
                        'seed': seed,
                        'output_name': output_name
                    }
                    combinations.append(combo)
                    job_number += 1
    
    return combinations


def print_manifest(combinations: List[Dict[str, Any]], file=sys.stdout) -> None:
    """Print combinations as a TSV manifest."""
    if not combinations:
        print("No combinations to display.", file=file)
        return
    
    # Header
    has_seeds = any(combo.get('seed') is not None for combo in combinations)
    has_ligand = any(combo.get('num_ligands', 0) > 0 for combo in combinations)
    
    columns = ['batch', 'job_number', 'file_index', 'protein', 'num_protomers']
    if has_ligand:
        columns.extend(['ligand', 'num_ligands'])
    if has_seeds:
        columns.append('seed')
    columns.append('output_name')
    
    print('\t'.join(columns), file=file)
    
    for combo in combinations:
        row = [str(combo.get(col, '')) for col in columns]
        print('\t'.join(row), file=file)


def main():
    """Main function to handle command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate AlphaFold 3 input JSON files for homomeric modeling",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate input for 2 protomers
  python generate_af3_homomeric.py ./inputs --protomers 2 --output ./outputs
  
  # Generate inputs for multiple protomer counts (creates separate jobs)
  python generate_af3_homomeric.py ./inputs --protomers 2 3 4 --output ./outputs
  
  # Add ATP ligand
  python generate_af3_homomeric.py ./inputs --protomers 2 --ligand ATP --output ./outputs
  
  # Add ligand by SMILES
  python generate_af3_homomeric.py ./inputs --protomers 2 --ligand "smiles:CCO" --output ./outputs
  
  # Add ligand from file
  python generate_af3_homomeric.py ./inputs --protomers 2 --ligand @ligand.txt --output ./outputs
  
  # Multiple ligand copies
  python generate_af3_homomeric.py ./inputs --protomers 2 --ligand ATP --num-ligands 2 --output ./outputs
  
  # Custom seeds
  python generate_af3_homomeric.py ./inputs --protomers 2 --seeds 1 2 3 --output ./outputs

=== Array Job Integration ===

  # Single-job mode (for SLURM array jobs)
  python generate_af3_homomeric.py ./inputs --index 0 --protomers 2 --output-file /tmp/job.json
  
  # List all combinations (for manifest generation)
  python generate_af3_homomeric.py ./inputs --protomers 2 3 --list-combinations --batch my_screen
        """
    )
    
    parser.add_argument(
        'input_dir',
        type=Path,
        help='Path to directory containing input JSON files with MSA/templates'
    )
    parser.add_argument(
        '--output', '-o',
        type=Path,
        default=None,
        help='Output directory for generated JSON files'
    )
    parser.add_argument(
        '--protomers', '-p',
        nargs='+',
        type=int,
        default=[1],
        help='Number of protomers (can specify multiple values for separate jobs)'
    )
    parser.add_argument(
        '--ligand', '-l',
        type=str,
        default=None,
        help='Ligand specification: CCD code (ATP), SMILES (smiles:CCO), or file (@ligand.txt)'
    )
    parser.add_argument(
        '--num-ligands', '-nl',
        nargs='+',
        type=int,
        default=None,
        help='Number of ligand copies (can specify multiple values)'
    )
    parser.add_argument(
        '--seeds', '-s',
        nargs='+',
        type=int,
        default=None,
        help='Custom seed numbers'
    )
    parser.add_argument(
        '--expand-seeds',
        action='store_true',
        help='Create separate jobs for each seed value'
    )
    
    # Array job integration
    array_group = parser.add_argument_group('Array Job Integration')
    array_group.add_argument(
        '--index', '-i',
        type=int,
        default=None,
        help='Index of file to process (0-based, for single-job mode)'
    )
    array_group.add_argument(
        '--output-file',
        type=Path,
        default=None,
        help='Exact output file path for single-job mode'
    )
    array_group.add_argument(
        '--job-number',
        type=int,
        default=None,
        help='Job number for naming'
    )
    array_group.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Minimize output'
    )
    
    # For array job: need to specify exact protomer count and ligand count
    array_group.add_argument(
        '--protomer-count',
        type=int,
        default=None,
        help='Exact protomer count for single-job mode'
    )
    array_group.add_argument(
        '--ligand-count',
        type=int,
        default=None,
        help='Exact ligand count for single-job mode'
    )
    
    # Manifest options
    manifest_group = parser.add_argument_group('Manifest Generation')
    manifest_group.add_argument(
        '--list-combinations',
        action='store_true',
        help='List all combinations as TSV without processing'
    )
    manifest_group.add_argument(
        '--batch',
        type=str,
        default=None,
        help='Batch name for manifest output'
    )
    
    args = parser.parse_args()
    
    # Validate input
    if not args.input_dir.exists():
        print(f"Error: Input path not found: {args.input_dir}", file=sys.stderr)
        sys.exit(1)
    
    # Parse ligand if specified
    ligand_info = None
    if args.ligand:
        ligand_info = parse_ligand_input(args.ligand)
    
    # Set default ligand counts
    ligand_counts = args.num_ligands
    if ligand_counts is None and ligand_info:
        ligand_counts = [1]
    elif ligand_counts is None:
        ligand_counts = [0]
    
    # === Mode: List combinations ===
    if args.list_combinations:
        combinations = list_combinations(
            args.input_dir,
            protomer_counts=args.protomers,
            ligand_info=ligand_info,
            ligand_counts=ligand_counts,
            custom_seeds=args.seeds,
            expand_seeds=args.expand_seeds,
            batch_name=args.batch
        )
        print_manifest(combinations)
        return
    
    # === Mode: Single-job (array job) ===
    if args.index is not None:
        if args.output_file is None:
            print("Error: --output-file is required when using --index", file=sys.stderr)
            sys.exit(1)
        
        # Use exact counts for single-job mode
        protomer_count = args.protomer_count if args.protomer_count else args.protomers[0]
        ligand_count = args.ligand_count if args.ligand_count is not None else (ligand_counts[0] if ligand_counts else 0)
        
        result = generate_single_job(
            args.input_dir,
            args.index,
            num_protomers=protomer_count,
            output_file=args.output_file,
            ligand_info=ligand_info,
            num_ligands=ligand_count,
            custom_seeds=args.seeds,
            job_number=args.job_number,
            quiet=args.quiet
        )
        
        # Print result as JSON for easy parsing
        if args.quiet:
            print(json.dumps(result))
        
        return
    
    # === Mode: Batch processing ===
    if args.output is None:
        print("Error: --output is required for batch processing", file=sys.stderr)
        sys.exit(1)
    
    # Generate all combinations
    combinations = list_combinations(
        args.input_dir,
        protomer_counts=args.protomers,
        ligand_info=ligand_info,
        ligand_counts=ligand_counts,
        custom_seeds=args.seeds,
        expand_seeds=args.expand_seeds,
        batch_name=args.batch
    )
    
    if not combinations:
        print("No combinations to process.", file=sys.stderr)
        sys.exit(1)
    
    print(f"Processing {len(combinations)} combinations...")
    args.output.mkdir(parents=True, exist_ok=True)
    
    files = get_json_files(args.input_dir)
    
    for combo in combinations:
        template_file = files[combo['file_index']]
        template_data = load_json(template_file)
        
        # Determine seeds for this job
        if combo['seed'] is not None:
            job_seeds = [combo['seed']]
        else:
            job_seeds = args.seeds
        
        # Generate input
        result_data = generate_homomeric_input(
            template_data, template_file,
            num_protomers=combo['num_protomers'],
            ligand_info=ligand_info if combo['num_ligands'] > 0 else None,
            num_ligands=combo['num_ligands'],
            custom_seeds=job_seeds,
            output_name=combo['output_name'],
            verbose=True
        )
        
        # Save
        output_file = args.output / f"{combo['output_name']}.json"
        save_json(result_data, output_file)
        print(f"  -> {output_file}")
    
    print(f"\nGenerated {len(combinations)} input files in {args.output}")


if __name__ == "__main__":
    main()
