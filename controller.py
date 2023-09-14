from ccmodule import *

convert_ent_pdb = False
run_extract_sequence = False
run_extract_conditions = True

'''
Preprocessing:
    - (FileZilla) Download all (*.ent.gz) files from ftp.wwpdb.org/pub/pdb/data/structures/all (24hrs)
    - Extract all downloaded files (*.ent.gz -> *.ent) (8-10hrs)
        - (7zip cmd): "C:/Program Files/7-zip/7z.exe" e"*.gz" -o../Extracted/
        - NOTE: all forward slashes need to be converted to backslashes for use in Windows command line...
    - Convert all *.ent files to *.pdb files simply by renaming (4 hours using script below)
        - Do this by setting the convert_ent_pdb variable above to True
        - Also need to set the file path to the directory containing the .ent files, which is below as ent_file_folder
'''
ent_file_folder = 'D:/PDB_extracted'
pdb_file_folder = 'D:/test'
fasta_file_folder = 'D:/test_fasta'

if convert_ent_pdb:
    convert = ConvertToPDB(ent_file_folder)
    if __name__ == '__main__':
        convert.handle_multiprocessing()

if run_extract_sequence:
    '''
    CREATE FASTA SEQUENCE FILES
    
    This section creates a PDBPreprocessingSequences object which is then used to create .fasta files from the SEQRES and
    ATOM sequences, which are extracted using the BioPython library from the .pdb files.
    '''
    atom_folder = Path("D:/Data/Protein_Sequences/ATOM")
    seq_res_folder = Path("D:/Data/Protein_Sequences/SEQRES/")
    extract = PDBPreprocessingSequences(pdb_file_folder, atom_folder, seq_res_folder)
    # This must be included in order to utilize multiprocessing on Windows...
    if __name__ == '__main__':
        extract.extract_sequences()

if run_extract_conditions:
    conditions = PDBCrystallizationConditions(pdb_file_folder)

    # This must be included in order to utilize multiprocessing on Windows...
    if __name__ == '__main__':
        conditions.initialize_dataframe(fasta_file_folder)
    # This must be included in order to utilize multiprocessing on Windows...
    if __name__ == '__main__':
        conditions.handle_multiprocessing()
