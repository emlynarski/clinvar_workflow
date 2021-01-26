# write_outputs.py

import os
from datetime import datetime


################################################################################
#### Write output files helper functions
################################################################################

def write_output_dir_helper(out_path, out_prefix, workflow_str):
	timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
	results_dir = out_prefix + workflow_str + timestamp
	
	#@TODO: refactor - use out_path2 variable for path
	if not os.path.exists(os.path.join(out_path, results_dir)):
		os.mkdir(os.path.join(out_path, results_dir))
	out_path2 = os.path.join(out_path, results_dir)
	
	return out_path2, timestamp


def write_df_file_helper(df, fname, out_dir, header=True, index=False, excel=True):
	"""
	
	Args:
		df:
		fname:
		out_dir:
		excel:
		header:
		index:

	Returns:

	"""
	## write tab-separated text file
	df.to_csv(os.path.join(out_dir, fname + '.txt'), header=header, index=index, sep='\t')
	
	## write excel file
	if excel:
		df.to_excel(os.path.join(out_dir, fname + '.xlsx'), header=header, index=index, merge_cells=False)



################################################################################
#### Write Annotation workflow Panda Dataframes functions
################################################################################

#@TODO: add param - subdir:False --> ONLY make annotation dir for EXPLORE
def write_annot_df_files(out_path, out_prefix, today, cv_full_df, cv_summ_df, excel=True):
	"""
	
	Args:
		out_path:
		out_prefix:
		today:
		cv_full_df:
		cv_summ_df:
		excel:

	Returns:

	"""
	
	if not os.path.exists(os.path.join(out_path,'annotation')):
		os.mkdir(os.path.join(out_path, 'annotation'))
	out_dir = os.path.join(out_path, 'annotation')
	
	## specify output file names
	fname_cv_full = out_prefix + '_ClinVar_variant_full_' + today
	fname_cv_summ = out_prefix + '_ClinVar_variant_summary_' + today
	
	## write Variant Summary DF output file
	print("\t.. Writing ClinVar Variant summary")
	write_df_file_helper(cv_summ_df, fname_cv_summ, out_dir, excel=excel)
	
	## write full ClinVar DF output file
	print("\t.. Writing ClinVar Variant full detailed DF")
	write_df_file_helper(cv_full_df, fname_cv_full, out_dir, excel=excel)


##----Driver function: write Annotation workflow output files-----------------##
#@TODO: rename fxn
def write_output_annotation(out_path, out_prefix, result_dict, excel=True):
	"""
	
	Args:
		out_path:
		out_prefix:
		result_dict:

	Returns:

	"""
	## create output directory
	output_path, timestamp = write_output_dir_helper(out_path, out_prefix,
	                                                 '_ClinVar_annotation_')
		
	## write DF output files
	write_annot_df_files(output_path, out_prefix=out_prefix, today=timestamp,
	                     cv_full_df=result_dict['cv_full_df'],
	                     cv_summ_df=result_dict['cv_var_summary_df'],
	                     excel=excel)
	


################################################################################
#### Write Exploratory Analysis workflow output files functions
################################################################################

##----write Panda Dataframes output files function----------------------------##
def write_explore_df_files(out_path, out_prefix, today, data_summ_df=None,
                           patho_var_df=None, patho_full_df=None, cs_var_df=None,
                           cs_rcv_df=None, cs_gene_df=None, cs_cond_df=None, excel=True):
	"""
	
	Args:
		out_path:
		out_prefix:
		today:
		data_summ_df:
		patho_var_df:
		patho_full_df:
		cs_var_df:
		cs_rcv_df:
		cs_gene_df:
		cs_cond_df:
		excel:

	Returns:

	"""
	
	if not os.path.exists(os.path.join(out_path,'exploratory_analysis')):
		os.mkdir(os.path.join(out_path, 'exploratory_analysis'))
	out_dir = os.path.join(out_path, 'exploratory_analysis')
	
	## specify output file names
	fname_data_summ = out_prefix + '_ClinVar_data_summary_' + today
	fname_patho_var = out_prefix + '_ClinVar_pathogenic_variants_' + today
	fname_patho_detail = out_prefix + '_ClinVar_pathogenic_variants_detailed_' + today
	fname_cs_var = out_prefix + '_ClinVar_clinsig-Variant_counts_' + today
	fname_cs_rcv = out_prefix + '_ClinVar_clinsig-RCV_counts_' + today
	fname_cs_gene = out_prefix + '_ClinVar_clinsig-Variant_counts_per_gene_' + today
	fname_cs_cond = out_prefix + '_ClinVar_clinsig-Variant_counts_per_condition_' + today
	
	## write DF output files
	_write = [
			[data_summ_df, fname_data_summ, False, True, "Dataset summary"],
			[patho_var_df, fname_patho_var, True, False, "Pathogenic Variant summary"],
			[patho_full_df, fname_patho_detail, True, False, "Pathogenic Variant details"],
			[cs_var_df, fname_cs_var, True, False, "Variant clinsig counts"],
			[cs_rcv_df, fname_cs_rcv, True, False, "RCV clinsig counts"],
			[cs_gene_df, fname_cs_gene, True, True, "Gene clinsig counts"],
			[cs_cond_df, fname_cs_cond, True, True, "Condition clinsig counts"]
	]

	for w in _write:
		if w[0] is not None:
			print("\t.. Writing " + w[4] + " DF files")
			write_df_file_helper(w[0], w[1], out_dir, header=w[2], index=w[3], excel=excel)


##----write Plot output files functions---------------------------------------##
def write_plot_files(out_path, out_prefix='', pie_var=None, pie_rcv=None,
					 bar_cond=None, bar_gene=None, today='',
					 write_plot_fxn=None):
	"""
	
	Args:
		out_path:
		out_prefix:
		pie_var:
		pie_rcv:
		bar_cond:
		bar_gene:
		today:
		write_plot_fxn:

	Returns:

	"""
	## create output directory
	if not os.path.exists(os.path.join(out_path,'images')):
		os.mkdir(os.path.join(out_path, 'images'))
	plot_path = os.path.join(out_path, 'images')
	
	if not os.path.exists(os.path.join(plot_path, 'pdf')):
		os.mkdir(os.path.join(plot_path, 'pdf'))
	
	if not os.path.exists(os.path.join(plot_path, 'png')):
		os.mkdir(os.path.join(plot_path, 'png'))
	
	## specify file names
	fname_pie_var = out_prefix + '_clinsig_Variant_donut_plot_' + today
	fname_pie_rcv = out_prefix + '_clinsig_RCV_donut_plot_' + today
	fname_bar_gene = out_prefix + '_Gene_clinsig_Variant_bar_plot_' + today
	fname_bar_cond = out_prefix + '_Condition_clinsig_Variant_bar_plot_' + today
	
	## write image files
	plot_list = [[pie_var, fname_pie_var, plot_path, "Variant clinsig donut"],
				 [pie_rcv, fname_pie_rcv, plot_path, "RCV clinsig donut"],
				 [bar_gene, fname_bar_gene, plot_path, "Gene Variant clinsig bar"],
				 [bar_cond, fname_bar_cond, plot_path, "Condition Variant clinsig bar"]]
	
	for p in plot_list:
		if p[0] is not None:
			print("\t.. Writing " + p[3] + " plot files")
			write_plot_fxn(p[0], p[1], p[2])


##----Driver function: write Exploratory Analysis workflow output files-------##
#@TODO: rename fxn
def write_output_exploratory_analysis(out_path, out_prefix, result_dict, write_plot_fxn):
	"""
	
	Args:
		out_path:
		out_prefix:
		result_dict:
		write_plot_fxn:

	Returns:

	"""
	## create output directory
	output_path, timestamp = write_output_dir_helper(out_path, out_prefix,
	                                                 '_ClinVar_analysis_')
	
	## write DF output files
	write_annot_df_files(output_path, out_prefix=out_prefix, today=timestamp,
	                     cv_full_df=result_dict['cv_full_df'],
	                     cv_summ_df=result_dict['cv_var_summary_df'])
	
	write_explore_df_files(output_path, out_prefix=out_prefix, today=timestamp,
	                       data_summ_df=result_dict['data_summary_df'],
	                       patho_var_df=result_dict['patho_var_df'],
	                       patho_full_df=result_dict['patho_var_detail_df'],
	                       cs_var_df=result_dict['clinsig_var']['df'],
	                       cs_rcv_df=result_dict['clinsig_rcv']['df'],
	                       cs_gene_df=result_dict['clinsig_var_gene']['df'],
	                       cs_cond_df=result_dict['clinsig_var_cond']['df'])
	
	## write plot image files
	write_plot_files(output_path, out_prefix=out_prefix, today=timestamp,
					 pie_var=result_dict['clinsig_var']['plot'],
					 pie_rcv=result_dict['clinsig_rcv']['plot'],
					 bar_cond=result_dict['clinsig_var_cond']['plot'],
					 bar_gene=result_dict['clinsig_var_gene']['plot'],
					 write_plot_fxn=write_plot_fxn)



