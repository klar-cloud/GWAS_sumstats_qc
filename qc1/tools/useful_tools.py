# In[]
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 14:27:13 2020

@author: klara
"""

import pickle as pkl

def save_obj(obj, name):
    with open(name, 'wb') as f:
        pkl.dump(obj, f, pkl.HIGHEST_PROTOCOL)

def load_obj(name):
    with open(name, 'rb') as f:
        return pkl.load(f)

# In[]
"""
Created on Tue Feb 25 14:27:13 2020

@author: from internet
"""


def flatten(lst):
  """ remove all nested list inside a list"""
  new_lst = []
  flatten_helper(lst, new_lst)
  return new_lst
 
def flatten_helper(lst, new_lst):
  for element in lst:
      if isinstance(element, list):
          flatten_helper(element, new_lst)
      else:
          new_lst.append(element)

# In[]
"""
Created on Tue Feb 25 14:27:13 2020

@author: from internet
"""
import numpy as np

def convert_dtype32(df, exception=None):
    
    for col in list(df.columns):
        if col != exception:
            dataTypeObj = df.dtypes[col]

            if dataTypeObj == np.int64:
                df = df.astype({col: np.int32})
            
            if dataTypeObj == np.float64:
                df = df.astype({col: np.float32})
            
    return df

# In[]:
def print_args(args):
    print("\n")
    print('{0:15s} {1:10s}'.format("Argument", "Value"))
    print('{0:15s} {1:10s}'.format("--------", "-----"))
    for arg, value in sorted(vars(args).items()):
        print('{0:15s} {1:10s}'.format(arg, str(value)))
    print("\n")