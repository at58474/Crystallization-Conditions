import os
import time
from pathlib import Path


class ConvertToPDB:
    """
    Convert all .ent files in a folder to .pdb by renaming.
    Multiprocessing is not used here because renaming is I/O-bound and
    dominated by filesystem metadata updates, so concurrency does not help.
    """

    def __init__(self, ent_file_folder):
        self.ent_file_folder = ent_file_folder

    @staticmethod
    def change_file_extension(file):
        base = os.path.splitext(file)[0]
        os.rename(file, base + '.pdb')

    def run(self):
        tic = time.time()
        ent_files = list(Path(self.ent_file_folder).glob("*.ent"))

        for f in ent_files:
            self.change_file_extension(str(f))

        toc = time.time()
        print(f"✅ Converted {len(ent_files)} files in {toc - tic:.2f} seconds")
