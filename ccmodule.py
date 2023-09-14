from Bio import SeqIO
from pathlib import Path
import re
import time
import multiprocessing
import os
import tqdm
import pandas as pd
import numpy as np
import difflib
import csv
from chemical_list import *


class PDBPreprocessingSequences:
    """
    This class is used to extract the protein sequence from a .pdb file. The sequences are extracted and converted into
    FASTA format using the BioPython library. These FASTA sequences are generated from the following:
        - SEQRES: the primary sequence of the protein
        - ATOM: The ATOM coordinated observed in the crystallized protein
            - Often times things like his-tags, chain ends, mobile loops, are cleaved. Thus, the amino acids that make
              up these features will not be present in the ATOM coordinates, but may be observed in the SEQRES sequence.
            - !!NOTE: Many warnings are given "PDBConstructionWarning: WARNING: Chain A is discontinuous at line 14891."
            - !!NOTE: The warnings were produced while extracting the FASTA sequence form the ATOM records. For the time
                      being only the sequences listed in the SEQRES section of the .pdb files will be used. REMARK 465
                      lists each missing residue which was likely removed/truncated. Add a function or class that lists
                      the missing residues in the .fasta files or a separate file.
            - This took roughly 2 hours to extract the SEQRES data, convert to FASTA, and write to new files, the files
              that were not successful are listed in the comments at the end of this module.
    ...
    Attributes
    ----------
    pdb_file_folder : str
        the path to the directory that contains the raw .pdb files
        NOTE: The path should not be in pathlib format! This path is created in the extract_sequences() method
    atom_folder : str
        pathlib path to the atom_folder where the ATOM files will be generated
    seqres_folder : str
        pathlib path to the seqres_folder where the ATOM files will be generated
    file : str
        This is the file that is passed to the create_sequence_file() method from the multiprocessing pool
    filename : str
        The filename generated from the current working file in the get_filename() method

    Methods
    -------
    get_filename(self):
        Grabs the basename of the current file then strips the .pdb file extension and saves into filename attribute
    create_sequence_file(self, file):
        - saves the file passed from the multiprocessing pool into file class attribute
        - calls get_filename method which stores the name of the current file into filename class attribute which is
          used to generate the output filename
        - Opens the file to extract the protein sequences using the BioPython SeqIO.parse method
        - Writes the protein sequence to the corresponding ATOM and SEQRES output file using SeqIO.write method
    extract_sequences(self):
        - Uses the pathlib library to store every file in the PDB directory into the pdb_files method variable
        - Creates a multiprocessing pool which is then used to map each file to the create_sequence_file() method,
          with intent to speed this process up. The reason multiprocessing is utilized is due to the large number
          of files contained in the PDB directory, roughly 300,000 files at 300GB. It is set to use every available
          CPU core. The runtime of this class is then displayed in the program output.

    Outputs
    -------
    *.fasta SEQRES files stored in '/Data/Protein_Sequences/SEQRES'
    *.fasta ATOM files stored in '/Data/Protein_Sequences/ATOM'


    """

    def __init__(self, pdb_file_folder, atom_folder, seqres_folder):
        self.pdb_file_folder = pdb_file_folder
        self.atom_folder = atom_folder
        self.seqres_folder = seqres_folder

        self.file = None
        self.filename = None

    def get_filename(self):
        # Uses os module to get the filename of the file, then strips the file extension
        # This is used to create a new file with the same name, but with .fasta file extension
        self.filename = os.path.basename(self.file.name)
        self.filename = self.filename.replace('.pdb', '')

    def create_sequence_file(self, file):
        self.file = file
        self.get_filename()

        try:
            # Using 'with open' automatically closes the opened file after the block has completed
            # It would be more efficient to only open the file once, but multiple BioPython functions are used to parse
            #   the .pdb file. Leaving this for a future time.
            '''
            !!! This extracts the FASTA sequence from the ATOM section, but with a high rate of failure, so excluding
            with open(file) as handle:
                sequence = next(SeqIO.parse(handle, "pdb-atom"))
            with open(self.atom_folder / f'{self.filename}_atom.fasta', "w") as output_handle:
                SeqIO.write(sequence, output_handle, "fasta")
            '''
            with open(file) as handle:
                sequence = next(SeqIO.parse(handle, "pdb-seqres"))
            with open(self.seqres_folder / f'{self.filename}_seqres.fasta', "w") as output_handle:
                SeqIO.write(sequence, output_handle, "fasta")

            return True

        except Exception as ex:
            file = str(Path(file))
            print(file + ' Not Successful ', ex)
            assert os.path.exists(file)
            return False

    def extract_sequences(self):
        # Using pathlib module to iterate through each file in the PDB directory
        pdb_files = Path(self.pdb_file_folder).glob('*.pdb')
        number_of_files = len([name for name in os.listdir(self.pdb_file_folder) if os.path.isfile(os.path.join(self.pdb_file_folder, name))])

        tic = time.time()
        with multiprocessing.Pool() as pool:
            results = pool.map(self.create_sequence_file, pdb_files)

        toc = time.time()
        print('Done in {:.4f} seconds'.format(toc - tic))

        assert all(results), f'{len(results) - sum(results)} out of {number_of_files} files were not processed'

        pool.close()


class PDBCrystallizationConditions:

    """
    Notes: For the crystallization conditions section, once something has been extracted, might want to consider
           removing that item from the list
                - convert list to string, delete, convert back to list
                - create a function that does this
                - find other common tasks and create functions for those, starting to get complicated

    """

    def __init__(self, pdb_file_folder):
        self.pdb_file_folder = pdb_file_folder

        # Stores the working file
        self.file = None
        # Stores the ID of the working protein
        self.protein_id = None
        # Just the name of the file with the file extension stripped, used for changing file names
        self.filename = None
        # This is the dictionary that contains all the extracted conditions which is sent back to the 'pool'
        self.conditions_dict = {}
        # List where all REMARK 280 lines are stored
        self.crystallization_conditions_list = []
        # List that contains only the data inside the crystallization conditions section of REMARK 280
        self.conditions = []
        # This is the main dataframe where all extracted data is stored
        self.conditions_df = pd.DataFrame(columns=['Protein_ID',
                                                   'FASTA_Sequence',
                                                   'Solvent_Concentration',
                                                   'Matthews_Coefficient',
                                                   'Crystallography_Type',
                                                   'Vapor_Diffusion_Type',
                                                   'pH',
                                                   'Organic_Precipitates',
                                                   'Salt_Precipitates'])
        # Variables that contain the extracted conditions
        self.solvent_content = None
        self.matthews_coefficient = None
        self.crystallography_type = None
        self.vapor_diffusion_type = None
        self.ph = None

    def get_filename_id(self):
        # Uses os module to get the filename of the file, then strips the file extension
        # This is used to create a new file with the same name, different file extension
        self.filename = os.path.basename(self.file.name)
        self.filename = self.filename.replace('.pdb', '')

        # In order to place each extracted condition in the correct row of the main dataframe, the protein ID needs
        #     to be matched. The easiest way is to extract the 4 character ID from the filename, just need to remove
        #     'pdb' and the remaining characters are the ID.
        self.protein_id = self.filename.replace('pdb', '').upper()

    def extract_remark_280(self):
        crystallization_conditions_list = []
        # Opening the given file in read-only mode
        with open(self.file, 'r') as filedata:
            # Traverse each line of the file
            for line in filedata:
                # Checking whether the given string is found in the line data
                if 'REMARK 280' in line:
                    # append the current line as a string to the crystallization_conditions_list
                    crystallization_conditions_list.append(line)
        if len(crystallization_conditions_list) > 0:
            self.crystallization_conditions_list = crystallization_conditions_list
        else:
            print(f'File {self.file} does not contain REMARK 280 content')

    # The most efficient way to append items to a dataframe is to use a list of dictionaries
    # For now:
    #       Step 1: Create a dataframe using the directory containing the extracted FASTA sequences where the unique ID
    #               will be the protein ID, along with the FASTA sequence in the 2nd column. This will create the
    #               base dataframe where the rest of the extracted informaiton will be included.
    #       Step 2: Open the files, extract the information, then append that information where a match is found
    #               based on the unique protein ID / key

    def append_id_fasta_to_dataframe(self, file):
        temp_dict = {}
        with open(file, 'r') as filedata:
            # Read all the lines of the file into a list
            # all_lines = filedata.readlines()
            all_lines = filedata.read().splitlines()

        # Gets the first line from the .fasta file and used a regex to extract only the protein ID from that line
        first_line = all_lines[0]
        protein_id = re.search(r'[A-Z0-9]{4}:[A-Z]', first_line)
        # An example protein ID is 1A0F:A, however only the first 2 characters are needed. Removing the :A
        protein_id = protein_id.group()
        protein_id = protein_id.split(':', 1)[0]

        # Using this conditional to handle a file that does not contain a protein ID.
        if protein_id is not None:
            temp_dict['Protein_ID'] = protein_id
        else:
            # Convert this into a try block and catch the exception
            print(f'File {file} does not have a protein ID')

        # This joins all lines containing the FASTA sequence into a single, unbroken string, and writes to the temp_dict
        temp_dict['FASTA_Sequence'] = ["".join(all_lines[1:])]

        # .append for dataframes is being deprecated so creating a temp dataframe to use .concat
        temp_df = pd.DataFrame(temp_dict)
        self.conditions_df = pd.concat([self.conditions_df, temp_df], axis=0, ignore_index=True)

    def initialize_dataframe(self, fasta_directory):
        # - Open all FASTA files
        # - Extract protein ID and FASTA sequence
        # - Create dataframe
        # Multiprocessing completely breaks down and gives strange and unexpected results here. For the conditions
        #     below, they will need to be saved to a file if I decide to go with multiprocessing, then convert that
        #     or read that file into a dictionary, then add to the dataframe. Otherwise, it will simply just not work.
        fasta_files = Path(fasta_directory).glob('*.fasta')
        for file in fasta_files:
            self.append_id_fasta_to_dataframe(file)

        # Set the index of the dataframe to the Protein_ID column, since each element should be unique, hopefully
        self.conditions_df = self.conditions_df.set_index('Protein_ID')

    def add_to_dataframe(self, column_name, condition):
        # For each condition that has been extracted, add to the dataframe
        self.conditions_df.at[self.protein_id, column_name] = condition

    def remove_from_conditions(self, string_to_remove):
        temp_list = []
        for line in self.conditions:
            line = line.replace(string_to_remove, '')
            temp_list.append(line)
        self.conditions = temp_list

    def extract_solvent_content(self):

        self.solvent_content = next(item for item in self.crystallization_conditions_list if 'SOLVENT CONTENT' in item)

        # remove the string 'REMARK 280' from solvent_content variable
        self.solvent_content = self.solvent_content.replace('REMARK 280', '')
        # extract only the solvent content percentage
        self.solvent_content = re.search(r'[+-]?([0-9]*[.])?[0-9]+', self.solvent_content)
        # check to see if a matching expression was found, if not then flag
        if self.solvent_content is not None:
            # store only the actual match in the solvent_content variable
            self.solvent_content = self.solvent_content.group()
            # print(self.solvent_content)
        else:
            self.solvent_content = 'No Solvent Content'

        # Add extracted solvent content to the dataframe
        self.add_to_dataframe('Solvent_Concentration', self.solvent_content)

    def extract_matthews_coefficient(self):
        self.matthews_coefficient = next(item for item in self.crystallization_conditions_list if 'MATTHEWS COEFFICIENT' in item)
        # remove the string 'REMARK 280' from matthews_coefficient variable
        self.matthews_coefficient = self.matthews_coefficient.replace('REMARK 280', '')
        # since there are multiple numbers will only save the data after the : character
        # Sample line: REMARK 280 MATTHEWS COEFFICIENT, VM (ANGSTROMS**3/DA): 3.07
        self.matthews_coefficient = self.matthews_coefficient.split(":", 1)[1]
        # extract only the matthews coefficient
        self.matthews_coefficient = re.search(r'[+-]?([0-9]*[.])?[0-9]+', self.matthews_coefficient)
        # check to see if a matching expression was found, if not then flag
        if self.matthews_coefficient is not None:
            # store only the actual match in the self.matthews_coefficient variable
            self.matthews_coefficient = self.matthews_coefficient.group()
            # print(self.matthews_coefficient)
        else:
            self.matthews_coefficient = 'No Matthews Coefficient'

        # Add extracted matthew's coefficient to the dataframe
        self.add_to_dataframe('Matthews_Coefficient', self.matthews_coefficient)

    def extract_crystallization_conditions(self):
        # Gets all list items starting with the item that contains the string 'CRYSTALLIZATION CONDITIONS', stores
        #   into list of lists
        self.conditions = []
        for item in self.crystallization_conditions_list:
            if 'CRYSTALLIZATION CONDITIONS' in item:
                item = item.replace('REMARK 280', '')
                # Removes leading or trailing whitespace
                item = item.strip()
                self.conditions.append([])
            if self.conditions:
                item = item.replace('REMARK 280', '')
                # Removes leading or trailing whitespace
                item = item.strip()
                self.conditions[-1].append(item)

        '''
        LIST PROCESSING
        '''
        # Joins all items in the list as a string
        self.conditions = ''.join(self.conditions[0])
        # Removes the string ' CRYSTALLIZATION CONDITIONS: '
        self.conditions = self.conditions.replace('CRYSTALLIZATION CONDITIONS: ', '')
        # Splits the self.conditions string into a list using ',' as the delimiter
        self.conditions = self.conditions.split(",")
        # Strips unwanted whitespace from each element in the list using list comprehension
        self.conditions = [item.strip() for item in self.conditions]

        # print(*self.conditions, sep='\n')

    def extract_crystallography_type(self):
        # The possible types are as follows (VAPOR DIFFUSION, MOLECULAR REPLACEMENT, FIBER SEEDING, DIALYSIS, BATCH,
        #   MICRODIALYSIS)
        # Create a list containing these items and search for these words in the self.condiditons list
        types_of_crystallography_list = ['VAPOR DIFFUSION', 'VAPORDIFFUSION', 'MOLECULAR REPLACEMENT', 'FIBER SEEDING',
                                         'DIALYSIS', 'BATCH', 'MICRODIALYSIS']

        # - If one of the crystallography types were not found then NaN needs to be entered in the dataframe row
        # - If one type was found, great, store it in the variable
        # - If more than one type was found, or multiple of the same type, no matter what store the 1st type
        for cond in self.conditions:
            if cond in types_of_crystallography_list:
                self.crystallography_type = cond
                self.remove_from_conditions(self.crystallography_type)
                break
            else:
                self.crystallography_type = 'No Crystallography Type'

        # Add extracted crystallography type to the dataframe
        self.add_to_dataframe('Crystallography_Type', self.crystallography_type)

    def extract_vapor_diffusion_type(self):
        types_of_vapor_diffusion = ['HANGING DROP', 'HANGING-DROP', 'HANGINGDROP', 'SITTING DROP', 'SITTING-DROP',
                                    'SITTINGDROP']
        for cond in self.conditions:
            if cond in types_of_vapor_diffusion:
                self.vapor_diffusion_type = cond
                self.remove_from_conditions(self.vapor_diffusion_type)
                break
            else:
                self.vapor_diffusion_type = 'No Vapor Diffusion Type'

        # Add extracted vapor diffusion type to the dataframe
        self.add_to_dataframe('Vapor_Diffusion_Type', self.vapor_diffusion_type)

    def extract_ph(self):
        for line in self.conditions:
            self.ph = re.search(r'PH\s?\b([1-9]|1[0-2])\b[.][0-9]', line)
            # check to see if a matching expression was found, if not then flag
            if self.ph is not None:
                # store only the actual match in the solvent_content variable
                self.ph = self.ph.group()
                break
            else:
                self.ph = 'No pH Content'

        # Add extracted solvent content to the dataframe
        self.add_to_dataframe('pH', self.ph)

    def extract_organic_precipitates(self):

        """
        This method executes a for loop that iterates through each line of the self.conditions list:
            - Searches the line to determine if it contains a percentage or percentage range, all of which could be
              either an integer or floating point number. The regular expression created to perform this task is as
              follows:
                - r'((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[%]'
                - this conditional statement uses re.search(expression, string) imported from the re library
                - the re.search() function call is wrapped in a bool function call, so if there is a REGEX match the
                  line will be processed
            - Then the line is split using re.split(expression, string). The expression used is the same as above but
              at the end includes the following:
                - A space or no space
                - Any set of characters which ends at a space
                - A space or no space
                - Any set of characters which ends at a space
              The goal of this is to capture any chemical name that comes after the percentage, since the
              chemical names are highly varied.
            - This regex has numerous capturing groups, a bool filter is applied to the line to remove list items
              such as '', none, NaN, False, etc.
            - Another for loop is used to iterate through each list item in the line list.
                - The same conditional statement is used to determine if the list item contains a percentage/range
                - If it does then that list item is appended to the organic_precipitates_list
        """

        organic_precipitates_list = []
        organic_precipitates_dict = {}

        for line in self.conditions:
            # If the line contains a percentage, then process the line
            if bool(re.search(r'((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[%]', line)):

                # Split the line based on the percentage regex
                # Must wrap the regex in a capturing group... () ... to keep the matched regex in the created list
                line = re.split(r'(((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[%]\s?([^\s]+)?\s?([^\s]+)?\s?([^\s]+)?)', line)
                line = list(filter(bool, line))

                for x in line:
                    if bool(re.search(r'((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[%]', x)):
                        organic_precipitates_list.append(x)

        def remove_bad_words(organic_prec_list):
            # Must be careful none of the words here are found in a chemical name.
            # Could add a space before and after, but sometimes spaces are forgotten in the files and these words are
            #   appended to a chemical name...
            words_to_exclude = ['AND', 'DEGREES', 'CONCENTRATION', 'PRIOR', 'DATA', ' THE', 'AQUEOUS']
            for item in organic_prec_list:
                index = organic_prec_list.index(item)
                for word in words_to_exclude:
                    if item.find(word) != -1:
                        item = item.replace(word, '')
                organic_prec_list[index] = item
            return organic_prec_list

        def create_dict(organic_prec_list):
            temp_dict = {}
            for item in organic_prec_list:
                # Find the percentage and store into percentage variable
                percentage = re.match(r'((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[%]', item)
                # Remove the percentage and store everything else in the item variable
                item = re.sub(r'((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[%]', '', item)
                temp_dict[percentage.group()] = item
            return temp_dict

        if organic_precipitates_list:
            # Remove unwanted words
            organic_precipitates_list = remove_bad_words(organic_precipitates_list)
            # Create list of dictionaries
            organic_precipitates_dict = create_dict(organic_precipitates_list)

            # If the chemical name is 2 characters or less, will assume it is erroneous and will remove from dict
            # Python 3 will not allow iterating through a dictionary and removing key:value pairs, so this looks weird
            for percentage in list(organic_precipitates_dict.keys()):
                if len(organic_precipitates_dict[percentage]) < 3:
                    # print(f'Deleting: {organic_precipitates_dict[percentage]}')
                    del organic_precipitates_dict[percentage]

            # Get the most similar match from the list of chemicals in the PDB and replace to create some sort of
            #   standardization
            for percentage, chemical in organic_precipitates_dict.items():
                find_closest_match = difflib.get_close_matches(chemical, chemical_list, 1, 0.2)
                if find_closest_match:
                    organic_precipitates_dict[percentage] = find_closest_match[0]
                    # print(f'Chemical: {chemical}')
                    # print(f'Closest Match: {find_closest_match[0]}')

            # Add the dictionary to the main conditions dataframe
            # NOTE: To do this the dictionary must be wrapped in a list
            organic_precipitates_list = [organic_precipitates_dict]
            self.add_to_dataframe('Organic_Precipitates', organic_precipitates_list)

            # print(organic_precipitates_dict)
            # print('\n')

    def extract_organic_precipitates_2(self):
        """
        This method executes a for loop that iterates through each line of the self.conditions list:
            - Searches the line to determine if it contains a percentage or percentage range, all of which could be
              either an integer or floating point number. The regular expression created to perform this task is as
              follows:
                - r'((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[%]'
                - this conditional statement uses re.search(expression, string) imported from the re library
                - the re.search() function call is wrapped in a bool function call, so if there is a REGEX match the
                  line will be processed
            - Then the line is split using re.split(expression, string). The expression used is the same as above but
              at the end includes the following:
                - A space or no space
                - Any set of characters which ends at a space
                - A space or no space
                - Any set of characters which ends at a space
              The goal of this is to capture any chemical name that comes after the percentage, since the
              chemical names are highly varied.
            - This regex has numerous capturing groups, a bool filter is applied to the line to remove list items
              such as '', none, NaN, False, etc.
            - Another for loop is used to iterate through each list item in the line list.
                - The same conditional statement is used to determine if the list item contains a percentage/range
                - If it does then that list item is appended to the organic_precipitates_list
        """
        organic_precipitates_list = []
        organic_precipitates_dict = {}

        for line in self.conditions:
            # If the line contains a percentage, then process the line
            if bool(re.search(
                    r'((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[%]',
                    line)):

                line = line.replace(',', '')
                line = line.replace('\'', '')
                line = line.replace(';', '')
                line = line.replace('SATURATED', '')

                # Using finditer() to find each match for molatity and store the match plus the next 15 characters into a list
                count = 0
                for match in re.finditer(
                        r'((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[%]',
                        line):
                    count += 1
                    organic_precipitates_list.append(line[match.start():match.end() + 20])

        def remove_bad_words(organic_prec_list):
            # Must be careful none of the words here are found in a chemical name.
            # Could add a space before and after, but sometimes spaces are forgotten in the files and these words are
            #   appended to a chemical name...
            words_to_exclude = ['AND', 'DEGREES', 'CONCENTRATION', 'PRIOR', 'DATA', ' THE', 'AQUEOUS']
            for word in words_to_exclude:
                for item in organic_prec_list:
                    index = organic_prec_list.index(item)
                    if item.find(word) != -1:
                        del organic_prec_list[index]
            return organic_prec_list

        def create_dict(organic_prec_list):
            temp_list = []
            temp_dict = {}

            for item in organic_prec_list:
                # Find the percentage and store into percentage variable
                percentage = re.match(
                    r'((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[%]',
                    item)

                # Remove the percentage and store everything else in the item variable
                item = re.sub(
                    r'((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[%]',
                    '', item)

                # Create a 1 item dictionary
                temp_dict[percentage.group()] = item

                # Add the dicitonary to the list that will be returned, for some reason have to create a copy first
                copy_dict = dict(temp_dict)
                temp_list.append(copy_dict)

                # Reset the dictionary
                temp_dict.clear()

            return temp_list

        def exact_match(chemical):
            # Iterate through each chemical in the chemical list and if a match is found then replace and send back
            for chem in chemical_list:
                if chemical == chem:
                    return chemical
                else:
                    pass
            return False

        # Need to iterate through each dictionary in the salt_precipitates list, find the closest chemical match for each
        #   value and replace
        def replace_with_most_similar(organic_prec_list):
            # Get the most similar match from the list of chemicals in the PDB and replace to create some sort of
            #   standardization
            for organic_dict in organic_prec_list:
                for molarity, chemical in organic_dict.items():
                    # - See if there is an exact match, if so then go ahead and replace
                    # - If there is not an exach match then use difflib to find the best match and replace with that

                    # find_exact_match will be the replaced chemical name returned by the exact_match function
                    find_exact_match = exact_match(chemical)

                    if not find_exact_match:
                        find_closest_match = difflib.get_close_matches(chemical, chemical_list, 1, 0.2)
                        if find_closest_match:
                            organic_dict[molarity] = find_closest_match[0]
                        else:
                            organic_dict[molarity] = None
                    elif find_exact_match:
                        organic_dict[molarity] = find_exact_match
            return organic_prec_list

        if organic_precipitates_list:
            # Remove unwanted words
            # organic_precipitates_list = remove_bad_words(organic_precipitates_list)

            # Create list of dictionaries
            organic_precipitates_list = create_dict(organic_precipitates_list)

            # If the chemical name is 2 characters or less, will assume it is erroneous and will remove from dict
            # Python 3 will not allow iterating through a dictionary and removing key:value pairs, so this looks weird
            for item in organic_precipitates_list:
                for percentage in list(item.keys()):
                    if len(item[percentage]) < 3:
                        # print(f'Deleting: {organic_precipitates_dict[percentage]}')
                        del item[percentage]

            # Get the most similar match from the list of chemicals in the PDB and replace to create some sort of
            #   standardization
            organic_precipitates_list = replace_with_most_similar(organic_precipitates_list)
            print(f'Organic Precipitates: {organic_precipitates_list} \n')

            # Add the dictionary to the main conditions dataframe
            # NOTE: To do this the dictionary must be wrapped in a list
            # organic_precipitates_list = [organic_precipitates_dict]
            self.add_to_dataframe('Organic_Precipitates', organic_precipitates_list)

            print('\n')

    def extract_salt_precipitates_2(self):
        salt_precipitates_list = []
        salt_precipitates_dict = {}

        for line in self.conditions:
            # If the line contains a molarity unit, then process the line
            if bool(re.search(
                    r'((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[M][M]?',
                    line)):
                temp = line
                # Remove unwanted characters that create confusion
                line = line.replace(',', '')
                line = line.replace('\'', '')
                line = line.replace(';', '')
                line = line.replace('PH ', '')

                # Need to either include fractional M values or convert them to MM first
                # Using sub() to do this
                #   - create a capture group by enclosing the regex in (), then referencing that group with \1, then adding M
                #   - at the end to convert the decimal M to MM, the period will not be included below.
                line = re.sub(
                    r'(\.((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[M])',
                    r'\1M', line)

                # Using finditer() to find each match for molatity and store the match plus the next 15 characters into a list
                count = 0
                for match in re.finditer(
                        r'((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[M][M]?',
                        line):
                    count += 1
                    # print("match", count, match.group(), "start index", match.start(), "End index", match.end())
                    salt_precipitates_list.append(line[match.start():match.end() + 15])

        '''
        !!! NOT USED but keeping to illustrate the error !!!
        def remove_if_contains(salt_prec_list):
            # Instead of making the REGEX even more complex, dealing with the matches where something like 5-7MG/ML was
            #   erroneously extracted
            words_to_exclude = ['MG/ML', 'DNA', 'PROTEIN', 'ML OF', ' ML']
            for item in salt_prec_list:
                index = salt_prec_list.index(item)
                for word in words_to_exclude:
                    if item.find(word) != -1:
                        del salt_prec_list[index]
                        # If item is not changed, if 2 words are found this will delete multiple list items erroneously
                        item = 'None'
            return salt_prec_list
        '''

        # Had to reverse the for loops here, if word is inside item then will either create erroneous deletions unless
        #   item is reset, but then only 1 item can be deleted. This was a better fix.
        def remove_if_contains(salt_prec_list):
            # Instead of making the REGEX even more complex, dealing with the matches where something like 5-7MG/ML was
            #   erroneously extracted
            words_to_exclude = ['MG/ML', 'DNA', 'PROTEIN', 'ML OF', ' ML']
            for word in words_to_exclude:
                for item in salt_prec_list:
                    index = salt_prec_list.index(item)
                    if item.find(word) != -1:
                        del salt_prec_list[index]
            return salt_prec_list

        def create_dict(salt_prec_list):
            temp_list = []
            temp_dict = {}

            for item in salt_prec_list:
                # Find the molarity and store into molarity variable
                molarity = re.match(
                    r'((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[M][M]?',
                    item)

                # Remove the percentage and store everything else in the item variable
                item = re.sub(
                    r'((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))?-?((100((\.|,)[0-9]{1,2})?)|([0-9]{1,2}((\.|,)[0-9]{0,2})?))\s?[M][M]?',
                    '', item)

                # Create a 1 item dictionary
                temp_dict[molarity.group()] = item

                # Add the dicitonary to the list that will be returned, for some reason have to create a copy first
                copy_dict = dict(temp_dict)
                temp_list.append(copy_dict)

                # Reset the dictionary
                temp_dict.clear()

            return temp_list

        def exact_match(chemical):
            # Iterate through each chemical in the chemical list and if a match is found then replace and send back
            for chem in chemical_list:
                if chemical == chem:
                    return chemical
                elif chemical.find(chem) != -1:
                    return chem
                else:
                    pass
            return False

        # Need to iterate through each dictionary in the salt_precipitates list, find the closest chemical match for each
        #   value and replace
        def replace_with_most_similar(salt_prec_list):
            # Get the most similar match from the list of chemicals in the PDB and replace to create some sort of
            #   standardization
            for salt_dict in salt_prec_list:
                for molarity, chemical in salt_dict.items():
                    # - See if there is an exact match, if so then go ahead and replace
                    # - If there is not an exach match then use difflib to find the best match and replace with that

                    # find_exact_match will be the replaced chemical name returned by the exact_match function
                    find_exact_match = exact_match(chemical)
                    if not find_exact_match:
                        find_closest_match = difflib.get_close_matches(chemical, chemical_list, 1, 0.2)
                        if find_closest_match:
                            salt_dict[molarity] = find_closest_match[0]
                        else:
                            salt_dict[molarity] = None
                    elif find_exact_match:
                        salt_dict[molarity] = find_exact_match
            return salt_prec_list

        if salt_precipitates_list:
            salt_precipitates_list = remove_if_contains(salt_precipitates_list)

            salt_precipitates_list = create_dict(salt_precipitates_list)

            salt_precipitates_list = replace_with_most_similar(salt_precipitates_list)

            print(f'Salt Precipitates: {salt_precipitates_list}')
            print('\n')

            self.add_to_dataframe('Salt_Precipitates', salt_precipitates_list)

    def extract_conditions(self, file):
        self.file = file
        self.get_filename_id()

        """
        This section opens the PDB file then stores every line in the REMARK 280 section as a string in a list
        """
        self.extract_remark_280()

        """
        This section finds and stores the solvent content information into the variable solvent_content
        If no solvent content is found in the file, then solvent_content is set to False
        """
        self.extract_solvent_content()

        """
        Extracts matthews coefficient, if available. If this is not found in the file then matthews_coefficient
        is set to False
        """
        self.extract_matthews_coefficient()

        '''
        This section finds and stores the Crystallization Conditions information into new_conditions_list variable 
        that contains the chemical composition of the mother liquor along with pH, temp, and crystallization method.
        
        The purpose is just to simplify the data being worked with and to get some sort of standardized formatting
        before extracting the precipitants, pH, crystallography type, vapor diffusion type, and temperature
        '''
        self.extract_crystallization_conditions()

        self.extract_crystallography_type()

        self.extract_vapor_diffusion_type()

        self.extract_ph()

        self.extract_organic_precipitates_2()

        self.extract_salt_precipitates_2()

        '''
        Things that need to happen:
            1) There are times when the proteins are crystallized in a solution, but then once harvested they are
               then soaked in a different solution. The chemicals that are reused, whether the same concentration
               is used or not, are easy to deal with ->
                    - Just delete the second occurrence since the crystallization solution is listed first in
                      REMARK 280, then the soaking solution is listed second, easy enough.
               The problem is dealing with new chemicals it may be soaked in. There is no evidence without further
               investigation new chemicals are used since the one file pulled lists the soaking solution as being
               the same chemicals but in different concentrations, need to look into this, hopefully new chemicals
               are not used.
               - If they are, then will need to create a function that catches the word 'SOAKED' before any processing
                 occurs and delete everything after that word and hope this causes no new problems.
            2) There are definitely erroneous chemical names, need to find as many as possible and find the best way
               to deal with them.
            3) Need to include temperature in the dataframe, temp is in Kelvin
        '''

    def handle_multiprocessing(self):
        # Using pathlib module to iterate through each file in the PDB directory
        pdb_files = Path(self.pdb_file_folder).glob('*.pdb')
        number_of_files = len([name for name in os.listdir(self.pdb_file_folder) if os.path.isfile(os.path.join(self.pdb_file_folder, name))])
        for file in pdb_files:
            self.extract_conditions(file)

        print(self.conditions_df.to_string())

        # Multiprocessing, more trouble than it's worth on Windows
        '''
        # Using pathlib module to iterate through each file in the PDB directory
        pdb_files = Path(self.pdb_file_folder).glob('*.pdb')
        number_of_files = len([name for name in os.listdir(self.pdb_file_folder) if os.path.isfile(os.path.join(self.pdb_file_folder, name))])

        tic = time.time()
        with multiprocessing.Pool() as pool:
            results = pool.map(self.extract_conditions, pdb_files)
            #print(results)

        # Call the function that appends each extracted condition to its appropriate row in the main dataframe
        # Pass the returned dictionary contained in the results variable
        self.create_dataframe(results)
        #print(self.conditions_df.to_string())

        toc = time.time()
        print('Done in {:.4f} seconds'.format(toc - tic))
        assert all(results), f'{len(results) - sum(results)} out of {number_of_files} files were not processed'

        pool.close()
        '''


class ConvertToPDB:

    def __init__(self, ent_file_folder):
        self.ent_file_folder = ent_file_folder

    @staticmethod
    def change_file_extension(file):
        base = os.path.splitext(file)[0]
        os.rename(file, base + '.pdb')

    def handle_multiprocessing(self):
        tic = time.time()
        ent_files = Path(self.ent_file_folder).glob('*.ent')
        pool = multiprocessing.Pool()
        pool.map(self.change_file_extension, ent_files)
        pool.close()
        toc = time.time()
        print('Done in {:.4f} seconds'.format(toc - tic))



'''
!!! These are the files where the FASTA sequence could not be extracted from SEQRES section
D:\PDB_extracted\pdb1c4s.pdb Not Successful  
D:\PDB_extracted\pdb1c58.pdb Not Successful  
D:\PDB_extracted\pdb1cap.pdb Not Successful  
D:\PDB_extracted\pdb1car.pdb Not Successful  
D:\PDB_extracted\pdb1aga.pdb Not Successful  
D:\PDB_extracted\pdb1ao2.pdb Not Successful  
D:\PDB_extracted\pdb1ao4.pdb Not Successful  
D:\PDB_extracted\pdb1ugt.pdb Not Successful  
D:\PDB_extracted\pdb2hya.pdb Not Successful  
D:\PDB_extracted\pdb2mk1.pdb Not Successful  
D:\PDB_extracted\pdb1hpn.pdb Not Successful  
D:\PDB_extracted\pdb1dey.pdb Not Successful  
D:\PDB_extracted\pdb1hua.pdb Not Successful  
D:\PDB_extracted\pdb1hya.pdb Not Successful  
D:\PDB_extracted\pdb2kqo.pdb Not Successful  
D:\PDB_extracted\pdb2bvk.pdb Not Successful  
D:\PDB_extracted\pdb1kes.pdb Not Successful  
D:\PDB_extracted\pdb2c4s.pdb Not Successful  
D:\PDB_extracted\pdb4hya.pdb Not Successful  
D:\PDB_extracted\pdb3hya.pdb Not Successful  
D:\PDB_extracted\pdb3iri.pdb Not Successful  
D:\PDB_extracted\pdb3irj.pdb Not Successful  
D:\PDB_extracted\pdb3irk.pdb Not Successful  
D:\PDB_extracted\pdb3irl.pdb Not Successful  
D:\PDB_extracted\pdb4hp7.pdb Not Successful  
D:\PDB_extracted\pdb6hrg.pdb Not Successful  'charmap' codec can't decode byte 0x81 in position 1752: character maps to <undefined>
D:\PDB_extracted\pdb7mgy.pdb Not Successful  'charmap' codec can't decode byte 0x8d in position 1916: character maps to <undefined>
D:\PDB_extracted\pdb7ks6.pdb Not Successful
'''