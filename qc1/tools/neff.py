#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:23:17 2020

@author: klara
"""

def calculate_neff(n1, n2):
    """ Calculate efficient sample size
    Input: 
    n1 - int or list with number of cases
    n2 - int or list with number of controls 
    if lists n1 and n2 must be of the same length"""
    

    if not isinstance(n1, list):
        n1 = [n1]; n2 = [n2]
    
    neff_total = 0    
    for i in range(len(n1)):
        neff = round(4 * (n1[i] + n2[i]) * (n1[i] / (n1[i] + n2[i])) * (n2[i] / (n1[i] + n2[i])))
        neff_total = neff_total + neff
        
    return(neff_total)
    

def add_neff(df,  n1, n2):
    
    if 'NEFF' in df:
        df = df.rename(columns={"NEFF": "N"})
        print(f'Detected column with efficient sample size')
        
    elif 'NCASES' in df and 'NCONTROLS' in df:
        
        df['N'] = round(4 * (df['NCASES'] + df['NCONTROLS']) * (df['NCASES'] / (df['NCASES'] + df['NCONTROLS'])) * (df['NCONTROLS'] / (df['NCASES'] + df['NCONTROLS'])))
        df = df.drop(['NCASES','NCONTROLS'], 1)
        print(f'Added efficient sample size')
            
    elif 'N' not in df and n1 != 0 and n2 != 0:
        neff = calculate_neff(n1, n2)
        df.insert(4, 'N', neff)
        print(f'Added efficient sample size')
  
    df = df.astype({'N': 'int32'})
    
    return df