#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function

## import ClinVar workflow submodules
from clinvar_workflow.query_clinvar import clinvar_query as cv_query
from clinvar_workflow.helpers.process_user_inputs import process_user_inputs
from clinvar_workflow.helpers.write_outputs import write_output_annotation

## Pandas - setup
import pandas as pd
pd.options.mode.chained_assignment = None
pd.set_option('display.max_columns', None)



################################################################################
#### ClinVar MyVariant & Plotly data viz variables
################################################################################

##----MyVariant ClinVar fields------------------------------------------------##
# COL_CLINSIG = 'clinical_significance'
# COL_COND = 'conditions.name'
# COL_RCV = 'accession'

def run_clinvar_annotation(var_file, out_dir, out_prefix, build, cols_var,
                           cols_input=None, write_output=True, write_excel=True):
	"""
	
	Args:
		var_file:
		out_dir:
		out_prefix:
		build:
		cols_var:
		cols_input:
		write_excel:

	Returns:

	"""
	## Step 1: verify & process user inputs
	print("\nStep 1: verify & process user inputs")
	input_var_df, _out_dir, _col_id, _cols_input = process_user_inputs(var_file=var_file,
																	   output_dir=out_dir,
																	   build=build,
																	   cols_var=cols_var,
																	   cols_input=cols_input)
	
	## Step 2: run MyVariant ClinVar query
	print("\n\nStep 2: run MyVariant ClinVar query")
	cv_df = cv_query.run_clinvar_query(input_var_df, build=build, col_id=_col_id)
	
	if cv_df is None:
		print("\nNo input variants found in ClinVar. Exiting program.")
		return None
	
	## Step 3: process MyVariant ClinVar query results
	print("\n\nStep 3: process MyVariant ClinVar query results")
	result_dict = cv_query.process_clinvar_query(cv_df, input_var_df,
	                                             cols_var=cols_var,
	                                             cols_input=_cols_input,
	                                             col_id=_col_id)
	
	# Step 4: write output files
	if write_output:
		print("\n\nStep 4: write output files")
		write_output_annotation(_out_dir, out_prefix, result_dict, excel=write_excel)
		return result_dict
	
	return {'result_dict':result_dict, '_col_id':_col_id, '_out_dir':_out_dir}






