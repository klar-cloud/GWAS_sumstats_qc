# In[]
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

            
# In[]:
"""
Created on Tue Feb 25 14:27:13 2020

@author: klara

requires: flatten.py

"""

import re
from tools.useful_tools import flatten
def re_add_col(df, col_search, col_add, pattern):
    """ replace the value of one column (or create one if it 
    does not already exist) with regex pattern found in another column
    
    df: dataframe
    col_search: string with name of the column that will be search
    col_replace: string with name of the column that will be replaced
    pattern: regex pattern in the form: r'pattern'
    """
    
    search = []    

    for value in df[col_search]:
        search.append(re.findall(pattern, value))
        
    df[col_add] = flatten(search)
    
    # Add option for NA if its not found
    return df



# In[]
"""
Created on Tue Feb 25 14:27:13 2020

@author: klara
"""
import numpy as np

def common_coordinate(df):
    """ 
    Make a column COORIDNATE which will be a represenattaion of SNP's 
    chromosome number and id number (in format: chrX:XXXX)
    
    """
    
    if 'BP' in df and 'CHR' in df:
        df['COORDINATE'] = df['CHR'].astype(str) + ':' + df['BP'].astype(str)
        
    # elif use regex to find chr and bp 
    else: 
        df = re_add_col(df, 'ID', 'COORDINATE', r'\d+:\d+')
        
        re_add_col(df, 'COORDINATE', 'CHR', r'(\d+):')
        df = df.astype({'CHR': 'int8'})
        
        re_add_col(df, 'COORDINATE', 'BP', r'\d+:(\d+)')
        df = df.astype({'BP': 'int32'})        
        
    df = df.loc[:, df.columns != 'ID'] 
        
    return df

# In[]
"""

@author: klara

"""
def missing_cnames(keys, id_info, n1, n2):
    """ 
    Find if any of the necessary columns is missing. 
    """    
    
    missing = []
    if 'A1' not in keys:
        missing.append('A1')
    if 'A2' not in keys:            
        missing.append('A2')
    if 'FREQ' not in keys:            
        missing.append('FREQ')
    if 'PVAL' not in keys:            
        missing.append('PVAL')
    if 'CHR' not in keys and not id_info:
        missing.append('CHR')
    if 'BP'  not in keys and not id_info:
        missing.append('BP')
    if 'BETA' not in keys:
        if 'ZSCORE' not in keys:
            missing.append('BETA')
    if 'SE' not in keys:
        if 'ZSCORE' not in keys:
            missing.append('SE')
    if 'ZSCORE' not in keys:
        if 'BETA' not in keys and 'SE' not in keys:
            missing.append('ZSCORE')
    if 'NEFF' not in keys and 'N' not in keys:
        if 'NCASES' not in keys or 'NCONTROLS' not in keys:
            if n1 == 0 or n2 ==0 :
                missing.append('N')
                
    return missing   
    