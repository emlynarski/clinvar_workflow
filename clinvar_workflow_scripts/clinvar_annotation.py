#!/usr/bin/env python
# coding: utf-8

import argparse, os, sys

def run_workflow(pkg_path, var_file, out_dir, out_prefix, build, cols_var, cols_input):
	## import ClinVar exploratory analysis workflow module
	print("\n\t .. importing exploratory analysis module")
	
	#@TODO: remove sys.path.insert
	sys.path.insert(0, os.path.abspath(pkg_path))
	from clinvar_workflow.workflows import annotation_workflow as cv
	
	## run ClinVar exploratory analysis
	print("\n\t .. running exploratory analysis")
	results = cv.run_clinvar_annotation(var_file=var_file,
	                                    out_dir=out_dir,
	                                    out_prefix=out_prefix,
	                                    build=build,
	                                    cols_var=cols_var,
	                                    cols_input=cols_input,
	                                    write_output=True)
	
	#@TODO: test for empty results BEFORE print
	if 'cv_var_summary_df' in results:
		print('\nClinVar annotation:', results['cv_var_summary_df'].head(3))
		
	return results


if __name__ == "__main__":
	print('\n\n\nStarted clinvar_workflow_script.py\n')
	
	parser = argparse.ArgumentParser(description='Driver script for the ClinVar_annotation workflow')
	parser.add_argument('--pkg_path', required=True, default='..',
	                    help='ClinVar workflow Python package absolute or relative path.')
	parser.add_argument('--var_file', required=True,
	                    help='The input variant file absolute or relative path.')
	parser.add_argument('--out_dir', required=True,
	                    help='The output directory absolute or relative path.')
	parser.add_argument('--out_prefix', required=True, default='',
	                    help='The prefix for the output directory & files. Default = \'\'')
	parser.add_argument('--build', required=True, default='hg19', choices=['hg19', 'hg38'],
	                    help='Specify current genome build: hg19 | hg38. Default = hg19.')
	parser.add_argument('--cols_var', required=True, default='CHR,POS,REF,ALT',
	                    help='The 4 Variant columns names: CHR (chromosome), POS (position), REF (reference allele) & ALT (alternative allele). The column names should be comma-separated. Default = \'CHR,POS,REF,ALT\'')
	parser.add_argument('--cols_input', required=False, default='',
	                    help='Optional: string containing a list of input columns to include in the output files. The column names should be comma-separated. Default = \'\'')

	## 1. Parse Args
	print("\n\t .. parsing args")
	pargs = parser.parse_args()
	CV_WKFLW_PKG_PATH = pargs.pkg_path
	VAR_FILE = pargs.var_file
	OUT_DIR = pargs.out_dir
	OUT_PREFIX = pargs.out_prefix
	BUILD = pargs.build

	## split COLS_VAR string --> list([Variant column names])
	COLS_VAR = [c.strip() for c in pargs.cols_var.split(',')]

	## split optional COLS_INPUT string --> list([input file column names])
	if len(pargs.cols_input) > 0:
		COLS_INPUT = [c.strip() for c in pargs.cols_input.split(',')]
	else:
		COLS_INPUT = []
	
	
	## 2. run ClinVar exploratory analysis
	run_workflow(pkg_path=CV_WKFLW_PKG_PATH,
	             var_file=VAR_FILE,
	             out_dir=OUT_DIR,
	             out_prefix=OUT_PREFIX,
	             build=BUILD,
	             cols_var=COLS_VAR,
	             cols_input=COLS_INPUT)
	
	
	## 3. exit
	print('\n\n\nclinvar_annotation.py complete. Goodbye.\n\n')
	exit(0)

