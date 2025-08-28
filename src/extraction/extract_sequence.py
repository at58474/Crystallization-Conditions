import os
import re
import time
from pathlib import Path
import multiprocessing


class PDBPreprocessingSequences:
    """
    This class handles preprocessing of PDB files to extract FASTA sequences.

    Notes:
        - The main task is to parse the SEQRES records in .pdb files and extract
          the protein sequences.
        - These sequences are then saved in FASTA format (.fasta files).
        - Multiprocessing is used to speed up processing when handling large numbers
          of files (hundreds of thousands).
        - The number of worker processes and chunksize can be configured to tune
          performance for different hardware setups (e.g. SSD vs HDD).

    Example:
        seq_extractor = PDBPreprocessingSequences("Data/PDB_Files", "Data/Protein_Sequences/SEQRES")
        seq_extractor.handle_multiprocessing()
    """

    def __init__(self, pdb_file_folder, fasta_file_folder, processes=None, chunksize=50):
        """
        Initialize the sequence extractor.

        Args:
            pdb_file_folder (str): Directory containing .pdb files.
            fasta_file_folder (str): Directory where .fasta files will be written.
            processes (int, optional): Number of worker processes to use.
                                       Defaults to (CPU count - 1).
            chunksize (int, optional): Number of files per task given to each worker.
                                       Larger chunks reduce IPC overhead.
        """
        self.pdb_file_folder = pdb_file_folder
        self.fasta_file_folder = fasta_file_folder
        self.processes = processes or max(1, multiprocessing.cpu_count() - 1)
        self.chunksize = chunksize

    def extract_sequences(self, pdb_file):
        """
        Extracts FASTA sequence from a single PDB file.

        - Extracts the 4-character protein ID from the filename.
        - Reads SEQRES lines in the PDB file, which contain amino acid residues.
        - Joins these residues into a continuous FASTA sequence.
        - Writes the sequence to a new .fasta file with the same protein ID.

        Args:
            pdb_file (str): Path to a .pdb file.
        """
        # Extract the protein ID from filename (strip 'pdb' prefix, keep 4-char ID)
        protein_id = os.path.basename(pdb_file).replace(".pdb", "")[3:].upper()
        fasta_path = os.path.join(self.fasta_file_folder, protein_id + ".fasta")

        seq_lines = []
        with open(pdb_file, "r") as f:
            for line in f:
                # SEQRES lines contain the primary sequence information
                if line.startswith("SEQRES"):
                    # Capture residues starting at column 20
                    seq = line[19:].strip().replace(" ", "")
                    seq_lines.append(seq)

        # If sequence lines were found, write them to a FASTA file
        if seq_lines:
            fasta_seq = "".join(seq_lines)
            with open(fasta_path, "w") as f:
                f.write(f">{protein_id}\n{fasta_seq}\n")

    def handle_multiprocessing(self):
        """
        Handles multiprocessing for batch extraction of FASTA sequences.

        - Collects all .pdb files from the pdb_file_folder.
        - Distributes files across worker processes.
        - Each worker calls extract_sequences() on its share of files.
        - Prints a runtime summary at the end.
        """
        tic = time.time()
        pdb_files = list(Path(self.pdb_file_folder).glob("*.pdb"))

        with multiprocessing.Pool(processes=self.processes) as pool:
            pool.map(self.extract_sequences, map(str, pdb_files), chunksize=self.chunksize)

        toc = time.time()
        print(f"✅ Extracted FASTA from {len(pdb_files)} files in {toc - tic:.2f} seconds")
