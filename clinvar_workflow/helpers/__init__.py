#!/usr/bin/env python
# coding: utf-8
# __init__.py

print(f'Invoking __init__.py for {__name__}')


from pathlib import PurePath
from .clinical_significance_config import create_clinsig_dicts


################################################################################
#### configure Clinical Significance settings
################################################################################

## Clinical Significance settings file
_file_clinsig_settings = 'clinsig_settings.txt'
_path_clinsig_settings = PurePath(__file__).parent.joinpath(PurePath(_file_clinsig_settings))


## generate Clinical Significance setting dictionaries from settings file
clinsig_sort_dict, clinsig_rgb_dict = create_clinsig_dicts(_path_clinsig_settings)


##----------------------------------------------------------------------------##

__all__ = [
		'process_user_inputs',
		'sorting',
		'summary_stats',
		'write_outputs',
		'clinsig_sort_dict',
		'clinsig_rgb_dict'
		]

