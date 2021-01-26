# process_user_inputs.py

import os

## Pandas - setup
import pandas as pd
pd.options.mode.chained_assignment = None


################################################################################
#### process user input functions
################################################################################

def test_user_input_variant_file(var_file, cols_var):
	"""Test the user supplied input variant file.

	Args:
		var_file:
		cols_var:

	Returns:

	"""
	## check if input file exists
	if os.path.isfile(var_file):
		print("\t.. User specified variant input file exists")
	else:
		print("ERROR: input variant file NOT found!!!")
		raise FileNotFoundError(
			'Input variant file NOT found! Update the "--var_file" value and rerun')
	
	## get input file extension
	ext = var_file.rpartition('.')[2]
	
	########################## @TODO ADD ERROR HANDLING #########################
	## read in input file
	if ext == 'xlsx':
		df = pd.read_excel(var_file)
	elif ext == 'csv':
		df = pd.read_csv(var_file, sep=',')
	else:
		df = pd.read_csv(var_file, sep='\t')
	##############################################################################
	
	## test that input DF contains the user specified variant columns
	cols_var_not_found = [c for c in cols_var if c not in df.columns]
	if len(cols_var_not_found) == 0:
		print("\t.. Input variant file contains the specified variant columns")
	else:
		print("\tcols_var_not_found:", cols_var_not_found)
		raise ValueError(
			'Input variant file does NOT contain the specified variant columns: ' + ', '.join(
				cols_var_not_found))
	return df


def test_user_output_directory(output_dir):
	"""Test the user supplied output directory.
	
	Determine if the user's output directory (a) exists, and (b) is writeable.

	Args:
		output_dir (str): Output directory user supplied path.

	Returns:
		str (path-like object): Output directory absolute path.

	"""
	## check if output directory exists && is writable
	if os.path.isdir(output_dir) & os.access(output_dir, os.W_OK):
		out_path = output_dir
		print("\t.. User specified output directory exists & is writable")
		
	else:
		if (os.path.isdir(output_dir)) & (os.access(output_dir, os.W_OK) == False):
			print("ERROR: cannot write to the specified output directory.")
		else:
			print("ERROR: output directory NOT found!!!")

		## try to write to input variant file directory instead
		if os.access(os.getcwd(), os.W_OK):
			out_path = os.getcwd()
			print("\t.. Output file(s) will be written to the current working directory instead.")
		else:
			raise PermissionError(
				'Write permission denied: specified output directory & current working directory')
	return os.path.abspath(out_path)



################################################################################
#### HGVS ID helper functions
################################################################################

def hgvs_helper_chr_int(variant_row, variant_cols):
	"""Helper function for add_hgvs_id_column(df, cols_var, col_hgvs).
	
	Use when Input variant 'CHROM' column does *NOT* contain 'chr' prefix.
	Apply this to Pandas DataFrame variant columns --> generate HGVS variant ID
	from the CHROM, POSITION, REF & ALT column values.
	
	Args:
		variant_row (Pandas DataFrame Row): Row corresponding to a single variant.
		variant_cols (List[str]): Names of the 4 variant columns.

	Returns:
		str: HGVS formatted ID for the current variant.

	"""
	hgvs_id = 'chr' + variant_row[variant_cols[0]].upper() \
	          + ':g.' + variant_row[variant_cols[1]] \
	          + variant_row[variant_cols[2]] + '>' \
	          + variant_row[variant_cols[3]]
	return hgvs_id


def hgvs_helper_chr_str(variant_row, variant_cols):
	"""Helper function for add_hgvs_id_column(df, cols_var, col_hgvs).
	
	Use when Input variant 'CHROM' column contains 'chr' prefix.
	Apply this to Pandas DataFrame variant columns --> generate HGVS variant ID
	from the CHROM, POSITION, REF & ALT column values.
	
	Args:
		variant_row (Pandas DataFrame Row): Row corresponding to a single variant.
		variant_cols (List[str]): Names of the 4 variant columns.

	Returns:
		str: HGVS formatted ID for the current variant.

	"""
	hgvs_id = variant_row[variant_cols[0]] + ':g.' \
	          + variant_row[variant_cols[1]] \
	          + variant_row[variant_cols[2]] + '>' \
	          + variant_row[variant_cols[3]]
	return hgvs_id


def add_hgvs_id_column(df, cols_var, col_hgvs):
	"""Add column to Input Variant DataFrame containing variant HGVS ID.
	
	Convert the CHROM, POSITION, REF & ALT columns to HGVS-format variant ID.
		Substitution (SNV): 'chr{CHROM}:g.{POSITION}{REF}>{ALT}'
		Deletion/Insertion (INDEL): 'chr{CHROM}:g.{POSITION}del{REF}ins{ALT}' @TODO
	
	Args:
		df (Pandas DataFrame): Input Variant DataFrame.
		cols_var (List[str]): List containing the variant column names.
		col_hgvs (str): Name for the new HGVS ID column.

	Returns:
		Pandas DataFrame: Input Variant DataFrame with new hgvs ID column.

	"""
	## strip any white space from variant columns
	for c in cols_var:
		df[c] = df[c].astype(str).str.strip()
	
	## convert variant columns to HGVS ID
	if 'chr' in str(df.loc[0, cols_var[0]]).lower():
		df[col_hgvs] = df[cols_var].apply(lambda x: hgvs_helper_chr_str(x, cols_var), axis=1)
	else:
		df[col_hgvs] = df[cols_var].apply(lambda x: hgvs_helper_chr_int(x, cols_var), axis=1)
	
	## cast POS column back to int
	df[cols_var[1]] = df[cols_var[1]].astype(int)
	return df



################################################################################
#### Driver function: run process user inputs
################################################################################

def process_user_inputs(var_file, output_dir, build, cols_var, cols_input=None):
	"""

	Args:
		var_file:
		output_dir:
		build:
		cols_var:
		cols_input:

	Returns:
		Pandas DataFrame:
		str: absolute path to output directory
		str: hgvs ID column name
		list: user specified input column names that are in input variant file

	"""
	## test user input variant file
	df = test_user_input_variant_file(var_file, cols_var)
	
	## test user output directory
	output_path = test_user_output_directory(output_dir)
	
	## add HGVS ID column
	col_hgvs = 'hgvs_id.' + build
	df = add_hgvs_id_column(df, cols_var, col_hgvs)
	
	## if user specified [optional] 'cols_input': test if columns in variant DF -> update 'cols_input'
	if cols_input is not None:
		opt_cols_not_found = [c for c in cols_input if c not in df.columns]
		cols_input = [c for c in cols_input if c in df.columns]
		if len(opt_cols_not_found) > 0:
			print("WARNING: user specified \'cols_input\' - column names NOT found: ", ', '.join(opt_cols_not_found))
	
	return df, output_path, col_hgvs, cols_input

