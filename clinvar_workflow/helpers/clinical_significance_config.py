#!/usr/bin/env python
# coding: utf-8

import pandas as pd

# _file_settings = 'clinsig_settings.txt'


## if you modify the clinsig settings file: UPDATE these to match column names
COL_LABEL = 'label'
COL_SORT = 'sort_order'
COL_COLOR = 'color_rgb'
COL_ALIAS = 'alias'


################################################################################
#### Helper functions
################################################################################

def dict_helper(col_key, col_keylow, col_value):
    _d_orig = dict(zip(col_key, col_value))
    _d_lower = dict(zip(col_keylow, col_value))
    return {**_d_orig, **_d_lower}


def clinsig_dict_helper(df, col_key, col_sort=COL_SORT, col_color=COL_COLOR):
    key_lower = df[col_key].copy().str.lower()
    ## create dicts
    sort_dict = dict_helper(df[col_key], key_lower, df[col_sort])
    color_dict = dict_helper(df[col_key], key_lower, df[col_color])
    return sort_dict, color_dict


def alias_helper(df, col_alias=COL_ALIAS):
    ## split & explode alias lists
    df.loc[:, 'alias2'] = df[col_alias].str.split(', ', expand=False)
    df = df.explode('alias2')
    df['alias2'] = df['alias2'].str.strip()
    
    ## create & return alias clinsig dicts
    return clinsig_dict_helper(df, col_key='alias2')



################################################################################
#### Driver function to create Clinical Significance dictionaries
################################################################################

def create_clinsig_dicts(settings_file, col_label=COL_LABEL, col_sort=COL_SORT, col_color=COL_COLOR, col_alias=COL_ALIAS):
    _df = pd.read_csv(settings_file, sep='\t')

    ## modify columns
    for c in [col_label, col_color, col_alias]:
        _df[c] = _df[c].fillna('').astype(str).str.strip()
    _df[col_sort] = _df[col_sort].astype(int)
    
    ## generate clinsig dicts with label column
    dict_sort, dict_color = clinsig_dict_helper(_df, col_key=col_label)
    
    ## check for aliases & add to dicts
    alias_df = _df[_df[col_alias] != ''].copy()
    if alias_df.shape[0] > 0:
        dict_sort_alias, dict_color_alias = alias_helper(alias_df, col_alias)
        dict_sort.update(dict_sort_alias)
        dict_color.update(dict_color_alias)
    
    return dict_sort, dict_color




