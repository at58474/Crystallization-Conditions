import os
from pathlib import Path
import time
import multiprocessing

ent_file_folder = 'C:/Users/adamo/OneDrive/Desktop/Projects/Crystallization_Conditions/Data/ENT_Files'


def change_file_extension(file):
    base = os.path.splitext(file)[0]
    os.rename(file, base + '.pdb')


if __name__ == '__main__':
    tic = time.time()
    ent_files = Path(ent_file_folder).glob('*.ent')
    pool = multiprocessing.Pool()
    pool.map(change_file_extension, ent_files)
    pool.close()
    toc = time.time()
    print('Done in {:.4f} seconds'.format(toc - tic))
