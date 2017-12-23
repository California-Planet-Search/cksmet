import pandas as pd
import numpy as np

def merge_spec_tables(df1, df2, on=['name','obs'], suffixes=['_sm','_lib']):
    """
    Merge spectroscopic tables
    
    Args:
        df1: first table
        df2: second table
        on: column(s) to perform merge
        suffixes: suffixes to append to merged column name

    Returns:
        pandas.DataFrame: combined table
    """

    comb = pd.merge(df1,df2,on=on,suffixes=suffixes)
    s0 = suffixes[0]
    s1 = suffixes[1]
    csm = [c.split('_')[0] for c in comb.columns if c.count(s0)>0]
    for c in csm:
        comb['%s_diff' % c] = comb["%s%s" % (c,s0)] - comb['%s%s' % (c,s1)]
    comb = comb[np.sort(list(comb.columns))]
    return comb 

def replace_columns(df,s1,s2):
    namemap = {}
    for k in df:
        if k.count(s1)>0:
            namemap[k] = k.replace(s1,s2)
    df = df.rename(columns=namemap)
    return df
