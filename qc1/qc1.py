
# In[DESCRIPTION]:
# =============================================================================
# This script: 
# - harmonizes labels
# - removes unnecessary columns
# - calculates efficient sample size
# - creates a coordinator column 
# - creates or updates reference panel
# - (additional) removes ambigues and multiallelic SNPs
# - (additional) removes SNPs with low MAF

# In[PACKAGES]:

import os
import glob
import gc
import pandas as pd
import sys
import argparse
import re

pd.set_option('display.max_columns', 5000)

# In[PARSERS]:

# Create the parser
parser = argparse.ArgumentParser()

# Add the 

parser.add_argument('--sumstats', type=str ,required=True,
                       help='A path to (unzipped) summary statistics; '
                       'only space/tab-separated values; example: "path/to/file/trait.csv"')

parser.add_argument('--cnames', type=str, required=True,
                       help='The path to file with column names pairs for harmonization; '
                       ' IMPORTANT: this file needs to be updated before starting new analysis;' 
                       'only space/tab-separated file; example: "path/to/file/cnames.csv"' )

parser.add_argument('--out', dest='save_dir', type=str, required=True,
                       help='A prefix (with path) for saving corrected sumstats; '
                       'example: "path/for/save" (without output file prefix)')

parser.add_argument('--ncases', dest='n1', nargs='+', type=int, default=0,
                       help='Size of a case sample; necessary if sumstats lacks sample size column(s); ' 
                       '(int or list of int if >1 stratum)')

parser.add_argument('--ncontrols', dest='n2', nargs='+', type=int, default=0,
                       help='Size of a control sample; necessary if sumstats lacks sample size column(s); ' 
                       '(int or list of int if >1 stratum)')

parser.add_argument('--remove-amb', action='store_true',
                    help='Do you want to remove ambigous SNPs?')

parser.add_argument('--remove-misp', action="store_true", 
                    help='Do you want to remove mispelled or multiallelic SNPs?')

parser.add_argument('--reverse-freq', action="store_true",
                    help='Do you want to reverse allelic frequency column values? '
                    '(1-FREQ; apply if only non-effect alleles frequency is present)')

parser.add_argument("--remove-maf", type=float, default=False,
                    help="Do you want to remove SNPs below certain frequency? Give MAF.")

parser.add_argument("--make-reference", type=str, default=None,
                    help="Do you want to create/update reference panel? Give prefix (and path) for saving it")

parser.add_argument("--test", action="store_true", default=False,
                    help="Do you want to load only first 100,000 SNPs for testing?")

# Execute the parse_args() method
args                        = parser.parse_args()

load_sumstats               = args.sumstats
load_cnames                 = args.cnames
save_dir                    = args.save_dir

n1                          = args.n1
n2                          = args.n2

params_remove_ambig         = args.remove_amb
params_remove_misp          = args.remove_misp
params_reverse_freq         = args.reverse_freq
params_remove_maf           = args.remove_maf
params_make_reference       = args.make_reference
params_test                 = args.test


if not os.path.exists(load_sumstats):
    print('ERROR: The specified load path for sumstats does not exist')
    quit()
if not os.path.exists(load_cnames):
    print('ERROR: The specified load path for cnames does not exist')
    quit()
#if not os.path.isdir(save_dir):
    #print('ERROR: The specified save path does not exist')
    #quit()     
#if params_make_reference is not None:
    #if not os.path.isdir(params_make_reference):
        #print('ERROR: The specified save path for reference does not exist')
        #quit()
# In[CUSTOM PACKAGES]

from tools.cnames import missing_cnames, common_coordinate, re_add_col
from tools.useful_tools import convert_dtype32, print_args
from tools.neff import add_neff

# In[LOG]
save_path_log = save_dir + '.qc1.log'

class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open(save_path_log, "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        pass    

# In[FIND COLUMNS NAMES]:

"""
This part finds necessary columns names. If the script cannot find all necessary columns it will stop.
It requires updated cnames file which contains also the description of necessary columns.
"""

sys.stdout = Logger()

print_args(args)

print(f'Quality check of summary statistics (part 1)\n')

cnames_dict = pd.read_csv(load_cnames, header=None, index_col=0, comment='#', squeeze=True, sep="\s+").to_dict()

# TAKE A CHUNK OF DATA TO PICK RIGHT LABELS
         
for trait_chunk in pd.read_csv(load_sumstats, chunksize=10, sep='\s+', error_bad_lines=False, engine='c'):

    print('Columns in the loaded dataset: ', '\n', trait_chunk.head(1) )
    
    # get columns names that interest you 
    correct_cnames = {}

    for col in trait_chunk:
        if col.upper() in cnames_dict:
            correct_cnames[cnames_dict[col.upper()]] = col
    
    keys_cnames = correct_cnames.keys()
    if 'ID' in keys_cnames:
        id_info = bool(re.search(r'(?<!rs)(\d+:\d+)', trait_chunk[correct_cnames['ID']].values[0]))
    else:
        id_info = False
        
    missing = missing_cnames(keys_cnames, id_info, n1, n2)
    missing = ' '.join([str(elem) for elem in missing]) 
    if len(missing) == 1:
        print(f'\nERROR: {missing} is a necessary column and it was not detected in the dataset. \nPlease correct the cnames file and try again.')
        quit()
    elif len(missing) > 1:
        print(f'\nERROR: {missing} are necessary columns and they ware not detected in the dataset. \nPlease correct the cnames file and try again.')
        quit()
    else:
        print(f'\nColumns chosen for further processing: \n {trait_chunk[list(correct_cnames.values())].head(1)}')
        print(f'All necessary columns are present')

        correct_cnames = {value:key for key, value in correct_cnames.items()}

        del trait_chunk                   
   
    break
    
# In[MAIN:
"""
This part harmonize labeles, clean alleles and do the general preparation for further processing.    
"""

chunk_list = []
chunk_num = 0

# LOAD SUMSTATS IN CHUNKS
print(f'\nLoading summary statistics...')
    
for chunk in pd.read_csv(load_sumstats,
                     chunksize=100000,
                     header=0, 
                     sep="\s+", 
                     error_bad_lines=False, 
                     engine='c'):
        
    chunk = chunk[list(correct_cnames.keys())]
    
    # Rename columns using the chosen columns names
    chunk = chunk.rename(columns=correct_cnames)
    
    # Convert all columns with 64 dtypes to 32 to save memory
    chunk = convert_dtype32(chunk, exception='PVAL')
    
    # Convert CHR column to 8bytes to save memory
    if 'CHR' in chunk:
        chunk = chunk[pd.to_numeric(chunk['CHR'], errors='coerce').notnull()]
        chunk = chunk.astype({'CHR': 'int8'})
    
    if 'BP' in chunk:
        chunk = chunk[pd.to_numeric(chunk['BP'], errors='coerce').notnull()]
        chunk = chunk.astype({'BP': 'int32'})
        
    chunk_list.append(chunk)
    
    if params_test:
        break
    
# Add all chunks to dataframe
trait = pd.concat(chunk_list)   

del chunk_list, chunk

print(f'Loaded {trait.shape[0]} SNPs')

if 'RSID' not in keys_cnames and re.search(r'(rs\d+)', trait['ID'][1]):
    re_add_col(trait, "ID", "RSID", r'(rs\d+)')

# Assure that nucleobases are in upper case
trait['A1'] = trait['A1'].str.upper()
trait['A2'] = trait['A2'].str.upper()

# Remove alleles letters different than A, T, C, G
if params_remove_misp:
    trait_len = trait.shape[0]
    
    base_names = ['A', 'C', 'G', 'T']        
    trait = trait.loc[
            trait['A1'].isin(base_names) & trait['A2'].isin(base_names)]
    
    print(f'Removed {trait_len-trait.shape[0]} multi-allelic and misspelled SNPs')
             
# Remove ambiguous SNPs 
if params_remove_ambig:
    trait_len = trait.shape[0]
    
    base_pairs = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    
    trait = trait.loc[
            ((trait['A1'] == 'A') & (trait['A2'] == 'G')) |
            ((trait['A1'] == 'A') & (trait['A2'] == 'C')) |
            ((trait['A1'] == 'T') & (trait['A2'] == 'G')) |
            ((trait['A1'] == 'T') & (trait['A2'] == 'C')) |
            ((trait['A1'] == 'G') & (trait['A2'] == 'A')) |
            ((trait['A1'] == 'G') & (trait['A2'] == 'T')) |
            ((trait['A1'] == 'C') & (trait['A2'] == 'A')) |
            ((trait['A1'] == 'C') & (trait['A2'] == 'T')) 
            ]
        
    print(f'Removed {trait_len-trait.shape[0]} ambiguous SNPs')
    
if params_reverse_freq:
    trait['FREQ'] = 1 - trait['FREQ']
    print(f'Reversed frequency column values')
    
if params_remove_maf:
    trait_len = trait.shape[0]
    trait = trait.loc[
            ~(trait['FREQ'] < params_remove_maf) &
            ~(1 - trait['FREQ'] < params_remove_maf)]

    print(f'Removed {trait_len-trait.shape[0]} SNPs with MAF lower than {params_remove_maf}')

# Add a column with coordinates for each row that consist of chromosome number and bp number
try:
    trait = common_coordinate(trait)
    print(f'Added common coordinates')

except:
    print('WARNING: Cannot create column with coordinates!')
    params_reference = False
    
# Add efficient smaple size or replace original sample size with it
trait = add_neff(trait, n1, n2)

print(f'\nFinal sumstats:', '\n', trait.head(1))

# SAVE SUMSTATS
save_path_trait = save_dir + '.qc1.csv'    

trait.to_csv(save_path_trait, sep=' ', mode='w', index=False)
print(f'Successfully saved {trait.shape[0]} SNPs\n')


# In[REF]:
# Create or update the REFERENCE PANEL   
 
if params_make_reference:
    trait =  trait[['A1', 'A2', 'FREQ', 'COORDINATE']]
        
    load_save_path_ref = params_make_reference + '.ref.csv'

    files_present = glob.glob(load_save_path_ref)

    if not files_present:
        trait_reference = trait
        print(f'Created reference panel')
    
    else: 
        trait_reference = pd.read_csv(load_save_path_ref,sep='\s+', header=0)
        trait_reference = pd.concat([trait_reference, trait], axis=0).drop_duplicates('COORDINATE').reset_index(drop=True)
        os.remove(load_save_path_ref)    
        print(f'Updated reference panel')

    trait_reference.to_csv(load_save_path_ref, sep=' ', mode='a', index=False, columns=['A1', 'A2', 'COORDINATE'])
    print(f'Successfully saved reference panel\n')

    del trait_reference

del trait
gc.collect()

#os.system('shutdown -s')