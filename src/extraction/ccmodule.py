"""
ccmodule.py

Module for extracting crystallization conditions from PDB files
(REMARK 280 records). This script parses solvent content,
Matthews coefficient, crystallography types, pH, and organic/salt
precipitates. Results are returned as dictionaries (SQLite-ready).

Relies on a `chemical_list` module for known chemical names.
"""

import os
import re
import difflib
import logging
from chemical_list import *  # contains chemical_list[]

logger = logging.getLogger(__name__)


class PDBCrystallizationConditions:
    """
    Extract crystallization conditions from a PDB file.

    Attributes:
        file (str): Path to the PDB file.
        protein_id (str): 4-letter PDB protein ID.
        filename (str): Base filename of the PDB file.
        conditions (list[str]): Cleaned crystallization condition strings.
        crystallization_conditions_list (list[str]): Raw REMARK 280 condition lines.
        solvent_content (str|None): Extracted solvent percentage.
        matthews_coefficient (str|None): Extracted Matthews coefficient.
        crystallography_type (str|None): Type of crystallography.
        vapor_diffusion_type (str|None): Vapor diffusion type.
        ph (str|None): Extracted pH value.
    """

    # Regex configs for organic vs. salt precipitants
    PRECIPITATE_CONFIG = {
        "organic": {
            "pattern": r'((100((\.|,)[0-9]{1,2})?)|'
                       r'([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?'
                       r'((100((\.|,)[0-9]{1,2})?)|'
                       r'([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[%]',
            "replacements": {",": "", "'": "", ";": "", "SATURATED": ""},
            "slice_extra": 20,
            "exclude_words": None,
            "allow_substring": False,
            "apply_fractional_fix": False,
        },
        "salt": {
            "pattern": r'((100((\.|,)[0-9]{1,2})?)|'
                       r'([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?'
                       r'((100((\.|,)[0-9]{1,2})?)|'
                       r'([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[M][M]?',
            "replacements": {",": "", "'": "", ";": "", "PH ": ""},
            "slice_extra": 15,
            "exclude_words": ['MG/ML', 'DNA', 'PROTEIN', 'ML OF', ' ML'],
            "allow_substring": True,
            "apply_fractional_fix": True,
        },
    }

    # Common regex patterns reused across extractors
    PATTERNS = {
        "number": re.compile(r'[+-]?([0-9]*[.])?[0-9]+'),
        "ph": re.compile(r'PH\s?\b([1-9]|1[0-2])\b[.][0-9]'),
    }

    def __init__(self, pdb_file=None):
        self.file = pdb_file
        self.protein_id = None
        self.filename = None
        self.conditions = []
        self.crystallization_conditions_list = []

        # Extracted values
        self.solvent_content = None
        self.matthews_coefficient = None
        self.crystallography_type = None
        self.vapor_diffusion_type = None
        self.ph = None

    # ---------------------------------------------------------------------
    # Utility helpers
    # ---------------------------------------------------------------------
    def _clean_line(self, line, replacements):
        """Apply text replacements to clean condition lines."""
        for old, new in replacements.items():
            line = line.replace(old, new)
        return line

    def _extract_matches(self, line, pattern, slice_extra):
        """Return regex matches + extra slice of text."""
        return [line[m.start():m.end() + slice_extra] for m in re.finditer(pattern, line)]

    def _create_dicts(self, items, pattern):
        """Convert raw items into list of {match: remainder} dicts."""
        results = []
        for item in items:
            key_match = re.match(pattern, item)
            if not key_match:
                continue
            key = key_match.group()
            value = re.sub(pattern, '', item)
            results.append({key: value})
        return results

    def _remove_short_values(self, dicts, min_len=3):
        """Drop dict entries with very short chemical names."""
        for d in dicts:
            for key in list(d.keys()):
                if len(d[key]) < min_len:
                    del d[key]
        return dicts

    def _remove_if_contains(self, items, words):
        """Remove items containing any excluded words."""
        return [item for item in items if not any(word in item for word in words)]

    def _exact_match(self, chemical, allow_substring=False):
        """Check chemical_list for exact or substring match."""
        for chem in chemical_list:
            if chemical == chem:
                return chemical
            if allow_substring and chem in chemical:
                return chem
        return False

    def _replace_with_most_similar(self, dicts, allow_substring=False):
        """Replace chemical names with closest known match."""
        for d in dicts:
            for key, val in d.items():
                found = self._exact_match(val, allow_substring)
                if not found:
                    close = difflib.get_close_matches(val, chemical_list, 1, 0.2)
                    d[key] = close[0] if close else None
                else:
                    d[key] = found
        return dicts

    def _extract_keyword(self, candidates):
        """Find and return the first candidate keyword in conditions."""
        for cond in self.conditions:
            if cond in candidates:
                self.remove_from_conditions(cond)
                return cond
        return None

    # ---------------------------------------------------------------------
    # File helpers
    # ---------------------------------------------------------------------
    def get_filename_id(self):
        """Extract filename and 4-letter protein ID from path."""
        self.filename = os.path.basename(self.file).replace('.pdb', '')
        self.protein_id = self.filename.replace('pdb', '').upper()

    def remove_from_conditions(self, string_to_remove):
        """Remove a given substring from all conditions."""
        self.conditions = [line.replace(string_to_remove, '') for line in self.conditions]

    # ---------------------------------------------------------------------
    # Extraction routines
    # ---------------------------------------------------------------------
    def extract_remark_280(self):
        """Extract all REMARK 280 lines from the file."""
        with open(self.file, "r", encoding="utf-8", errors="ignore") as f:
            lines = [line for line in f if "REMARK 280" in line]
        self.crystallization_conditions_list = lines

    def extract_solvent_content(self):
        """Extract solvent content percentage."""
        try:
            line = next(item for item in self.crystallization_conditions_list if 'SOLVENT CONTENT' in item)
            line = line.replace("REMARK 280", "")
            match = self.PATTERNS["number"].search(line)
            self.solvent_content = match.group() if match else "No Solvent Content"
        except StopIteration:
            self.solvent_content = "No Solvent Content"

    def extract_matthews_coefficient(self):
        """Extract Matthews coefficient (after colon in REMARK 280 line)."""
        try:
            line = next(item for item in self.crystallization_conditions_list if 'MATTHEWS COEFFICIENT' in item)
            line = line.replace("REMARK 280", "")
            line = line.split(":", 1)[1]  # take text after colon
            match = self.PATTERNS["number"].search(line)
            self.matthews_coefficient = match.group() if match else None
        except StopIteration:
            self.matthews_coefficient = None

    def extract_crystallization_conditions(self):
        """Extract all conditions into a normalized list of strings."""
        self.conditions = []
        for item in self.crystallization_conditions_list:
            if "CRYSTALLIZATION CONDITIONS" in item:
                self.conditions.append([])
            if self.conditions:
                item = item.replace("REMARK 280", "").strip()
                self.conditions[-1].append(item)

        if self.conditions:
            cond_str = "".join(self.conditions[0])
            cond_str = cond_str.replace("CRYSTALLIZATION CONDITIONS: ", "")
            self.conditions = [c.strip() for c in cond_str.split(",")]

    def extract_crystallography_type(self):
        """Extract crystallography method keyword (if present)."""
        types = ['VAPOR DIFFUSION', 'VAPORDIFFUSION', 'MOLECULAR REPLACEMENT',
                 'FIBER SEEDING', 'DIALYSIS', 'BATCH', 'MICRODIALYSIS']
        self.crystallography_type = self._extract_keyword(types)

    def extract_vapor_diffusion_type(self):
        """Extract vapor diffusion type (if present)."""
        types = ['HANGING DROP', 'HANGING-DROP', 'HANGINGDROP',
                 'SITTING DROP', 'SITTING-DROP', 'SITTINGDROP']
        self.vapor_diffusion_type = self._extract_keyword(types)

    def extract_ph(self):
        """Extract pH value from crystallization conditions."""
        self.ph = None
        for line in self.conditions:
            match = self.PATTERNS["ph"].search(line)
            if match:
                self.ph = match.group()
                break

    def extract_precipitates(self, mode="organic"):
        """
        Extract organic or salt precipitants.

        Args:
            mode (str): "organic" or "salt".

        Returns:
            str|None: repr(list of dicts) or None if no matches found.
        """
        cfg = self.PRECIPITATE_CONFIG[mode]
        pattern = cfg["pattern"]

        precip_list = []
        for line in self.conditions:
            if re.search(pattern, line):
                line = self._clean_line(line, cfg["replacements"])
                if cfg["apply_fractional_fix"]:
                    # Normalize fractional M units
                    line = re.sub(
                        r'(\.((100((\.|,)[0-9]{1,2})?)|'
                        r'([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?'
                        r'((100((\.|,)[0-9]{1,2})?)|'
                        r'([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[M])',
                        r'\1M', line)
                precip_list.extend(self._extract_matches(line, pattern, cfg["slice_extra"]))

        if precip_list:
            if cfg["exclude_words"]:
                precip_list = self._remove_if_contains(precip_list, cfg["exclude_words"])
            precip_list = self._create_dicts(precip_list, pattern)
            if mode == "organic":
                precip_list = self._remove_short_values(precip_list)
            precip_list = self._replace_with_most_similar(
                precip_list, allow_substring=cfg["allow_substring"])
        return repr(precip_list) if precip_list else None

    def extract_conditions(self):
        """
        Run full extraction pipeline for this file.

        Returns:
            dict: Extracted crystallization data.
        """
        self.get_filename_id()
        self.extract_remark_280()
        self.extract_solvent_content()
        self.extract_matthews_coefficient()
        self.extract_crystallization_conditions()
        self.extract_crystallography_type()
        self.extract_vapor_diffusion_type()
        self.extract_ph()
        organic = self.extract_precipitates("organic")
        salts = self.extract_precipitates("salt")

        return {
            "Protein_ID": self.protein_id,
            "Solvent_Concentration": self.solvent_content,
            "Matthews_Coefficient": self.matthews_coefficient,
            "Crystallography_Type": self.crystallography_type,
            "Vapor_Diffusion_Type": self.vapor_diffusion_type,
            "pH": self.ph,
            "Organic_Precipitates": organic,
            "Salt_Precipitates": salts,
        }


def process_pdb_file(file):
    """Standalone helper for multiprocessing: extract conditions from one file."""
    extractor = PDBCrystallizationConditions(file)
    return extractor.extract_conditions()
