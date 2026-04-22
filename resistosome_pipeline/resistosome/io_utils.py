"""
File discovery, compression handling, and replicate grouping.
"""

import os
import re
import gzip
import tarfile
import tempfile
import shutil
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

# ── Patterns ──
SEED_PATTERN = re.compile(r'_seed(\d+)')
STRUCTURE_EXTS = {'.pdb', '.cif', '.pdb.gz', '.cif.gz'}
AF3_TAR_EXTS = {'.tar.gz', '.tar'}


def detect_seed(filename: str) -> Optional[int]:
    """Extract seed/replicate number from filename."""
    m = SEED_PATTERN.search(filename)
    return int(m.group(1)) if m else None


def strip_seed(filename: str) -> str:
    """Remove seed suffix to get the base model name."""
    return SEED_PATTERN.sub('', filename)


def normalize_id_from_filename(filename: str,
                               protein_names: Optional[List[str]] = None) -> str:
    """
    Extract a normalized protein ID from a filename.

    If protein_names is provided, matches the filename against known names
    (longest match wins). This is the recommended approach for flexible naming.

    Otherwise falls back to regex-based stripping of common AF3 suffixes.
    Always returns lowercase.
    """
    if protein_names:
        match = match_protein_name(filename, protein_names)
        if match:
            return match

    # Fallback: regex-based extraction
    name = Path(filename).stem
    # Remove common extensions that might remain
    for ext in ['.tar', '.gz', '.pdb', '.cif']:
        if name.endswith(ext):
            name = name[:-len(ext)]
    # Remove leading job_N_ prefix
    name = re.sub(r'^job_\d+_', '', name)
    # Remove trailing _seedN
    name = SEED_PATTERN.sub('', name)
    # Remove protomer/ligand specs like _6protomers_25OLALigs
    name = re.sub(r'_\d+protomers.*$', '', name, flags=re.IGNORECASE)
    return name.lower()


def match_protein_name(filename: str,
                       protein_names: List[str]) -> Optional[str]:
    """
    Match a filename against a list of known protein names.
    Returns the longest matching name (lowercase), or None if no match.
    Case-insensitive substring matching.
    """
    fname_lower = filename.lower()
    best = None
    for name in protein_names:
        nl = name.lower()
        if nl in fname_lower:
            if best is None or len(nl) > len(best):
                best = nl
    return best


def load_protein_names(path: str) -> List[str]:
    """
    Load protein names from a text file (one name per line).
    Strips whitespace, skips empty lines and comments (#).
    """
    names = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                names.append(line)
    logger.info(f"Loaded {len(names)} protein names from {path}")
    return names


def is_af3_tar(filepath: str) -> bool:
    """Check if a file is an AF3 tar(.gz) archive."""
    fp = filepath.lower()
    return fp.endswith('.tar.gz') or fp.endswith('.tar')


def is_structure_file(filepath: str) -> bool:
    """Check if a file is a PDB or CIF structure file (possibly gzipped)."""
    fp = filepath.lower()
    return (fp.endswith('.pdb') or fp.endswith('.cif') or
            fp.endswith('.pdb.gz') or fp.endswith('.cif.gz'))


class StructureInput:
    """Represents a single structure input (one seed/replicate)."""

    def __init__(self, source_path: str, seed: Optional[int] = None,
                 protein_names: Optional[List[str]] = None):
        self.source_path = source_path
        self.seed = seed
        self.is_af3 = is_af3_tar(source_path)
        self._is_af3_dir = os.path.isdir(source_path) and _dir_has_model(source_path)
        self.base_id = normalize_id_from_filename(
            os.path.basename(source_path), protein_names=protein_names)
        self._tmpdir = None

        # These will be populated after extraction
        self.model_path: Optional[str] = None         # path to CIF/PDB file
        self.confidences_path: Optional[str] = None   # path to full confidences JSON
        self.summary_conf_path: Optional[str] = None  # path to summary confidences JSON
        self.ranking_path: Optional[str] = None        # path to ranking_scores CSV
        self.model_format: Optional[str] = None        # 'cif' or 'pdb'

    def prepare(self, work_dir: str) -> 'StructureInput':
        """
        Extract/decompress files into a working directory.
        For AF3 tar: extract root-level files only.
        For AF3 directory: scan for model/confidence files.
        For plain gz: decompress.
        For plain files: just reference them.
        """
        self._tmpdir = tempfile.mkdtemp(dir=work_dir, prefix='struct_')

        if self.is_af3:
            self._extract_af3_tar()
        elif self._is_af3_dir:
            self._scan_af3_dir()
        elif self.source_path.lower().endswith('.gz'):
            self._decompress_gz()
        else:
            # Plain PDB/CIF file
            self.model_path = self.source_path
            self.model_format = 'cif' if self.source_path.lower().endswith('.cif') else 'pdb'

        return self

    def _scan_af3_dir(self):
        """Scan an AF3 output directory for model and confidence files."""
        for fname in os.listdir(self.source_path):
            fpath = os.path.join(self.source_path, fname)
            if not os.path.isfile(fpath):
                continue
            fl = fname.lower()
            if fl.endswith('_model.cif') or (fl.endswith('.cif') and self.model_path is None):
                self.model_path = fpath
                self.model_format = 'cif'
            elif fl.endswith('_model.pdb') or (fl.endswith('.pdb') and self.model_path is None):
                self.model_path = fpath
                self.model_format = 'pdb'
            elif 'summary_confidences' in fl and fl.endswith('.json'):
                self.summary_conf_path = fpath
            elif 'confidences' in fl and fl.endswith('.json'):
                self.confidences_path = fpath
            elif 'ranking_scores' in fl and fl.endswith('.csv'):
                self.ranking_path = fpath

    def _extract_af3_tar(self):
        """Extract root-level files from AF3 tar archive."""
        open_mode = 'r:gz' if self.source_path.endswith('.gz') else 'r:'
        try:
            with tarfile.open(self.source_path, open_mode) as tar:
                for member in tar.getmembers():
                    # Only root-level files (no directory separators in name
                    # after stripping the top-level directory if present)
                    parts = Path(member.name).parts
                    # AF3 tars often have a top-level dir; we want files at depth 0 or 1
                    if member.isfile() and len(parts) <= 2:
                        # Check if it's in a subdirectory like seed-1_sample-0/
                        if len(parts) == 2:
                            parent = parts[0]
                            if re.match(r'seed-\d+_sample-\d+', parent):
                                continue  # skip sample subdirectories
                        fname = parts[-1].lower()
                        # Extract to tmpdir with just the filename
                        member.name = parts[-1]
                        tar.extract(member, self._tmpdir)
                        extracted_path = os.path.join(self._tmpdir, parts[-1])

                        if fname.endswith('_model.cif') or fname.endswith('.cif'):
                            self.model_path = extracted_path
                            self.model_format = 'cif'
                        elif fname.endswith('_model.pdb') or fname.endswith('.pdb'):
                            self.model_path = extracted_path
                            self.model_format = 'pdb'
                        elif 'summary_confidences' in fname and fname.endswith('.json'):
                            self.summary_conf_path = extracted_path
                        elif 'confidences' in fname and fname.endswith('.json'):
                            self.confidences_path = extracted_path
                        elif 'ranking_scores' in fname and fname.endswith('.csv'):
                            self.ranking_path = extracted_path
        except Exception as e:
            logger.error(f"Failed to extract AF3 tar {self.source_path}: {e}")

    def _decompress_gz(self):
        """Decompress a gzipped structure file."""
        basename = os.path.basename(self.source_path)
        if basename.endswith('.gz'):
            basename = basename[:-3]
        out_path = os.path.join(self._tmpdir, basename)
        try:
            with gzip.open(self.source_path, 'rb') as f_in:
                with open(out_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            self.model_path = out_path
            self.model_format = 'cif' if out_path.lower().endswith('.cif') else 'pdb'
        except Exception as e:
            logger.error(f"Failed to decompress {self.source_path}: {e}")

    def cleanup(self):
        """Remove temporary files."""
        if self._tmpdir and os.path.exists(self._tmpdir):
            shutil.rmtree(self._tmpdir, ignore_errors=True)


def _dir_has_model(dirpath: str) -> bool:
    """Check if a directory contains an AF3 model file (CIF or PDB)."""
    try:
        for fname in os.listdir(dirpath):
            fl = fname.lower()
            if fl.endswith('_model.cif') or fl.endswith('_model.pdb'):
                return True
            # Also accept plain .cif/.pdb if no _model variant
            if fl.endswith('.cif') or fl.endswith('.pdb'):
                return True
    except OSError:
        pass
    return False


def discover_inputs(input_dir: str,
                    protein_names: Optional[List[str]] = None) -> Dict[str, List[StructureInput]]:
    """
    Scan a directory for structure files and group by base model ID.

    Handles:
      - Plain .pdb/.cif files (and .gz compressed)
      - AF3 tar.gz archives
      - AF3 uncompressed output directories (containing *_model.cif etc.)
      - Nested directories (walks up to 2 levels deep to find AF3 folders)

    Parameters
    ----------
    input_dir : str
        Directory to scan.
    protein_names : list of str, optional
        Known protein names. If provided, filenames are matched against
        these names instead of regex-based ID extraction. Longest match wins.

    Returns {base_id: [StructureInput, ...]} grouped by replicate.
    """
    groups: Dict[str, List[StructureInput]] = {}
    found_items = []

    # Walk input directory (up to 2 levels: input_dir/level1/level2)
    for root, dirs, files in os.walk(input_dir):
        depth = root.replace(input_dir, '').count(os.sep)
        if depth > 1:
            dirs.clear()  # don't recurse deeper
            continue

        # Skip AF3 sample subdirectories (seed-N_sample-N)
        dirs[:] = [d for d in dirs
                   if not re.match(r'seed-\d+_sample-\d+', d)]

        for fname in sorted(files):
            fpath = os.path.join(root, fname)

            if is_af3_tar(fname):
                seed = detect_seed(fname)
                si = StructureInput(fpath, seed=seed, protein_names=protein_names)
                groups.setdefault(si.base_id, []).append(si)
                found_items.append(('tar', fname, si.base_id, seed))

            elif is_structure_file(fname):
                # Only pick up structure files that are at the root level
                # (not inside AF3 output dirs which are handled as directories)
                if root == input_dir:
                    seed = detect_seed(fname)
                    si = StructureInput(fpath, seed=seed, protein_names=protein_names)
                    groups.setdefault(si.base_id, []).append(si)
                    found_items.append(('file', fname, si.base_id, seed))

        # Check subdirectories at this level for AF3 output folders
        for dname in sorted(dirs):
            dpath = os.path.join(root, dname)
            if _dir_has_model(dpath):
                seed = detect_seed(dname)
                si = StructureInput(dpath, seed=seed, protein_names=protein_names)
                groups.setdefault(si.base_id, []).append(si)
                found_items.append(('dir', dname, si.base_id, seed))

        # Only process dirs at the current walk level (os.walk handles recursion)

    # Sort each group by seed
    for base_id in groups:
        groups[base_id].sort(key=lambda s: s.seed if s.seed is not None else 0)

    n_total = sum(len(v) for v in groups.values())
    logger.info(f"Discovered {n_total} structures "
                f"in {len(groups)} groups from {input_dir}")

    if n_total == 0:
        # Debug: list what we actually saw in the input directory
        try:
            entries = os.listdir(input_dir)
            logger.warning(f"  Input directory contains {len(entries)} entries: "
                           f"{entries[:20]}{'...' if len(entries) > 20 else ''}")
            for e in entries[:10]:
                epath = os.path.join(input_dir, e)
                etype = 'dir' if os.path.isdir(epath) else 'file'
                logger.warning(f"    {etype}: {e}")
                if os.path.isdir(epath):
                    sub_files = [f for f in os.listdir(epath) if os.path.isfile(os.path.join(epath, f))]
                    has_model = _dir_has_model(epath)
                    logger.warning(f"      has_model={has_model}, files: {sub_files[:5]}")
        except Exception as ex:
            logger.warning(f"  Could not list input directory: {ex}")
    else:
        for kind, name, base_id, seed in found_items[:5]:
            logger.debug(f"  [{kind}] {name} → base_id={base_id}, seed={seed}")
        if len(found_items) > 5:
            logger.debug(f"  ... and {len(found_items) - 5} more")

    return groups
