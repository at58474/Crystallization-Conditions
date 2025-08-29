# Crystallization Conditions Extraction

This project extracts, processes, and stores **crystallization condition data** from PDB files, including solvent content, Matthews coefficient, pH, crystallography type, vapor diffusion type, and organic/salt precipitants. It supports **multiprocessing** for large datasets and optional CSV/database export.

---

## 📂 Project Structure

- **`ccmodule.py`**  
  Core logic for parsing `.pdb` files, extracting crystallization conditions, and managing the main conditions dataframe.  
  Includes detailed methods with docstrings for each extraction step.

- **`controller.py`**  
  Orchestrates the pipeline: runs the extractor on a folder of `.pdb` files, manages multiprocessing, and handles optional CSV export.  

- **`chemical_list.py`**  
  Provides a global list of valid chemical names for fuzzy matching and normalization.  

- **`convert_ent_to_pdb.py`** (optional)  
  Converts `.ent` files into `.pdb` format for compatibility with the extractor.  

- **`extract_sequence.py`** (optional)  
  Extracts FASTA sequences from `.pdb` files into `.fasta` format.  

- **`tests/`**  
  Contains unit tests (`pytest`) verifying the correctness of parsing logic.  

---

## Usage

### 1. Preprocessing (optional)
If your input files are `.ent` or you need FASTA sequences:

```bash
python convert_ent_to_pdb.py --input data/ent_files --output data/pdb_files
python extract_sequence.py --input data/pdb_files --output data/fasta_files
```

### 2. Extract crystallization conditions
Run the main controller:

```bash
python controller.py --pdb-folder data/pdb_files --export-csv true
```

- **`--pdb-folder`**: path containing `.pdb` files.  
- **`--export-csv`**: set `true` to export results as `conditions.csv` in the output folder.  

Results are always stored in **SQLite** for scalability; CSV export is optional.  

---

## Tests

Run all unit tests with:

```bash
pytest
```

- Includes checks for solvent content, Matthews coefficient parsing, crystallography type, vapor diffusion type, pH, and precipitate extraction.  
- All logic is covered by **unit tests**; integration tests with synthetic `.pdb` files are optional.  

---

## Performance

- Multiprocessing is used to process thousands of PDB files efficiently.  
- SQLite database backend ensures scalable storage of results.  
- CSV export is batched for speed.  

For extremely large datasets (e.g., >500k files):
- Adjust **`--processes`** and **`--chunksize`** in `controller.py`.  
- Consider running in HPC or cloud environments.  

---

## Key Extraction Logic

- **Solvent Content**: Extracted from `REMARK 280 SOLVENT CONTENT`.  
- **Matthews Coefficient**: Extracted after `:` in `REMARK 280 MATTHEWS COEFFICIENT`.  
- **Crystallization Conditions**: Extracted from `REMARK 280 CRYSTALLIZATION CONDITIONS`.  
- **Crystallography Type**: Matches against known crystallography methods.  
- **Vapor Diffusion Type**: Identifies hanging/sitting drop methods.  
- **pH**: Regex match for values like `PH 7.2`.  
- **Organic & Salt Precipitates**: Regex search with fuzzy chemical name matching from `chemical_list.py`.  

---

## 🔧 Requirements

- Python 3.8+  
- Dependencies: see `requirements.txt`

Install all requirements:

```bash
pip install -r requirements.txt
```

---

## Example Output

**SQLite table (`conditions`)**:

| Protein_ID | FASTA_Sequence | Solvent_Concentration | Matthews_Coefficient | Crystallography_Type | Vapor_Diffusion_Type | pH  | Organic_Precipitates | Salt_Precipitates |
|------------|----------------|------------------------|----------------------|----------------------|----------------------|-----|----------------------|-------------------|
| 1XYZ       | MTEYKLVVVG...  | 48.5                   | 3.12                 | VAPOR DIFFUSION      | HANGING DROP         | 7.2 | [{"10%": "PEG 4000"}]| [{"0.2M": "NaCl"}] |

---
