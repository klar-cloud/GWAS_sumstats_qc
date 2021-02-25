#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri May 22 13:11:29 2020

@author: klarg
"""


    
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
    if 'CHR' not in keys:
        if 'ID' not in keys or ('ID' in keys and not id_info):
            missing.append('CHR')
    if 'BP' not in keys:
        if 'ID' not in keys or ('ID' in keys and not id_info):
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
            if n1 is None or n2 is None:
                missing.append('N')
                
    return missing   
    
