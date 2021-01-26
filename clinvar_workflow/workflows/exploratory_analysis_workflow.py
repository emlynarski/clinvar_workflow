#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function

## import ClinVar workflow submodules
from clinvar_workflow.query_clinvar import clinvar_query as cv_query
from clinvar_workflow.vizualization import viz_static as viz
from clinvar_workflow.helpers import clinsig_sort_dict

# from clinvar_workflow.helpers.process_user_inputs import process_user_inputs
from clinvar_workflow.workflows import annotation_workflow as annot
from clinvar_workflow.helpers.write_outputs import write_output_exploratory_analysis
from clinvar_workflow.helpers import sorting as sort
from clinvar_workflow.helpers import summary_stats as ss


## Pandas - setup
import pandas as pd
pd.options.mode.chained_assignment = None
pd.set_option('display.max_columns', None)



################################################################################
#### ClinVar MyVariant & Plotly data viz variables
################################################################################

##----MyVariant ClinVar fields------------------------------------------------##
COL_CLINSIG = 'clinical_significance'
COL_COND = 'conditions.name'
COL_RCV = 'accession'


################################################################################
#### Exploratory analysis functions
################################################################################

## reformat full variant ClinVar DF for VISUALIZATION
def preprocess_visualization_df(full_df, col_id, col_clinsig):
	"""

	Args:
		full_df:
		col_id:
		col_clinsig:

	Returns:

	"""
	v_df = full_df.copy()
	v_df['condition'] = v_df['conditions.name'].copy()
	v_df['gene'] = v_df['gene.symbol'].copy()
	
	## new variant-RCV ID column used in
	v_df['var_rcv'] = v_df[col_id] + '_' + v_df['accession']
	
	## reformat clinsig string labels
	v_df[col_clinsig] = v_df[col_clinsig].str.strip()\
		.str.replace('Conflicting interpretations of pathogenicity', 'Conflicting')
	v_df[col_clinsig+'.rcv'] = v_df[col_clinsig+'.rcv'].str.strip()\
		.str.replace('Conflicting interpretations of pathogenicity', 'Conflicting')
	
	v_df[col_clinsig] = v_df[col_clinsig].str.strip()\
		.str.replace('Benign/Likely benign', 'Benign / Likely benign')
	v_df[col_clinsig+'.rcv'] = v_df[col_clinsig+'.rcv'].str.strip()\
		.str.replace('Benign/Likely benign', 'Benign / Likely benign')
	
	v_df[col_clinsig] = v_df[col_clinsig].str.strip()\
		.str.replace('Pathogenic/Likely pathogenic', 'Pathogenic / Likely pathogenic')
	v_df[col_clinsig+'.rcv'] = v_df[col_clinsig+'.rcv'].str.strip()\
		.str.replace('Pathogenic/Likely pathogenic', 'Pathogenic / Likely pathogenic')
	
	
	## add new clinsig col with 'Conflicting [some|none] pathogenic'
	v_df[col_clinsig+'2'] = v_df[col_clinsig].copy()
	
	conf_patho = v_df[v_df[col_clinsig] == 'Conflicting']\
					.groupby(col_id)\
					.filter(lambda x: any(x[col_clinsig+'.rcv']\
	                                      .str.contains('Pathogenic')))\
					[col_id].tolist()
	
	v_df.loc[v_df[col_id].isin(conf_patho), col_clinsig+'2'] = 'Conflicting - some pathogenic'
	v_df.loc[(v_df[col_clinsig] == 'Conflicting') &
			 (~v_df[col_id].isin(conf_patho)), col_clinsig + '2'] = 'Conflicting - none pathogenic'
	return v_df
	


def get_dataset_summary(input_df, cv_full_df, col_id, col_cond, col_gene):
	"""
	
	Args:
		input_df:
		cv_full_df:
		col_id:
		col_cond:
		col_gene:

	Returns:

	"""
	## # of variants: input file
	agg_input = input_df[[col_id]].nunique().reset_index(name='nuniq')
	agg_input['row_cat'] = 'Input variants'
	agg_input['row_label'] = 'Total #'

	## # of variants: ClinVar query
	cv_status_label_dict = {'reported':'# in ClinVar database',
							'NOT in ClinVar':'# NOT currently in ClinVar'}
	agg_cv_status = cv_full_df.groupby('clinvar_status')[col_id]\
							.agg([('nuniq', 'nunique')]).reset_index()
	agg_cv_status = agg_cv_status.sort_values(by='clinvar_status', ascending=False)\
								.reset_index(drop=True)
	agg_cv_status['order_label'] = agg_cv_status.index.astype('int') + 1
	agg_cv_status['row_label'] = agg_cv_status['clinvar_status'].map(cv_status_label_dict)
	agg_cv_status['row_cat'] = 'Input variants'

	## # of distinct genes in dataset
	agg_gene = cv_full_df[[col_gene]].nunique().reset_index(name='nuniq')
	agg_gene['row_cat'] = 'Genes'
	agg_gene['row_label'] = '# of distinct symbols'

	## # of distinct conditions in dataset
	agg_cond = cv_full_df[[col_cond]].nunique().reset_index(name='nuniq')
	agg_cond['row_cat'] = 'Conditions'
	agg_cond['row_label'] = '# of distinct ClinVar conditions'

	## concat DFs
	agg_cond_gene = pd.concat([agg_cond, agg_gene], axis=0, sort=False)
	agg_cond_gene = agg_cond_gene.sort_values('index', ascending=True)\
								.reset_index(drop=True)
	agg_cond_gene['order_cat'] = agg_cond_gene.index.astype('int') + 1

	summ_df = pd.concat([agg_input, agg_cv_status.drop(columns=['clinvar_status']),
	                     agg_cond_gene], axis=0, sort=False)
	summ_df['order_label'] = summ_df['order_label'].fillna(0)
	summ_df['order_cat'] = summ_df['order_cat'].fillna(0)
	summ_df = summ_df.sort_values(['order_cat', 'order_label'])\
					.reset_index(drop=True)\
					.drop(columns=['index'])

	## multi index DF for output file
	summ_df_midx = summ_df[['row_cat', 'row_label', 'nuniq']]\
					.set_index(['row_cat', 'row_label'])
	summ_df_midx.rename(columns={'nuniq':'  '}, inplace=True)
	midx = summ_df_midx.index.rename(['', ' '])
	return summ_df_midx.reindex(midx)


def identify_patho_vars(cv_summ_df, cv_full_df, col_id, col_clinsig, cols_var, all_cond=False):
	## determine which pathogenic flag column to use for filtering
	if all_cond:
		col_patho_flag = '_.patho_ALL_cond'
	else:
		col_patho_flag = '_.patho_ANY_cond'

	## filter / select pathogenic variants using the pathogenic flag column
	patho_var_df = cv_summ_df[cv_summ_df[col_patho_flag]==True].copy().reset_index(drop=True)
	patho_detail_df = cv_full_df[cv_full_df[col_patho_flag]==True].copy()

	## specify output columns
	########################## @TODO: make dynamic!
	cols_tmp = [col_id, 'rsid', 'gene.symbol', col_clinsig, col_clinsig+'.rcv.set', 'conditions.name',
				'conditions.name.set', 'conditions.synonyms', 'patho_cond.set',
				'patho_cond.nuniq', 'patho_cond.%', '_.patho_ALL_cond',
				'_.patho_ANY_cond', 'preferred_name', 'hg19.start', 'hg19.end',
				'hg38.start', 'hg38.end']
	return patho_var_df[cols_var+cols_tmp], patho_detail_df


def explore_clinsig(df, clinsig_col, count_col, unreported=False, plot_fxn=viz.plot_donut_annot_legend):
	## generate DF for donut plot
	_cs_totals_df = ss.generate_count_df(df, cnt_col=count_col, grp_cols=[clinsig_col])

	## sort ClinSig rows
	_cs_totals_df = sort.sort_clinsig_df_col(_cs_totals_df, clinsig_col, clinsig_dict=clinsig_sort_dict, reverse=True)
	
	## Capitalize ClinSig rows
	_cs_totals_df[clinsig_col] = _cs_totals_df[clinsig_col].str.capitalize().str.replace('/ likely', '/ Likely')

	## subset count DF for plot +/- 'UNREPORTED' variants
	if unreported:
		donut_plot_df = _cs_totals_df
	else:
		donut_plot_df = _cs_totals_df[_cs_totals_df[clinsig_col].str.lower() != 'unreported']

	## dynamically set up Clinical Significance column name
	if 'rcv' in clinsig_col:
		clinsig_col_out = "RCV-level clinical significance"
	else:
		clinsig_col_out = "Variant clinical significance"

	## dynamically set up count column name
	if 'rcv' in count_col:
		count_col_out = "# of distinct Variant-RCV pairs"
	elif 'rcv' in clinsig_col:
		count_col_out = "# of Variants"
	else:
		count_col_out = "# of distinct Variants"

	## rename output DF columns
	_cs_totals_df.rename(columns={clinsig_col:clinsig_col_out, 'nuniq':count_col_out}, inplace=True)

	## generate clinsig donut plot
	donut_plot_fig = viz.plot_clinsig_donut(donut_plot_df, clinsig_col, plot_fxn)
	return _cs_totals_df, donut_plot_fig


def explore_clinsig_by_gene(df, col_clinsig, col_count, col_gene, axis_order='total'):
	## generate DF for plot
	plot_df = ss.generate_grouped_count_df(df, col_group=col_gene, col_count=col_count, col_clinsig=col_clinsig)

	## sort DF columns
	cols_df = ['Total'] + sort.sort_and_extract_clinsig(plot_df.columns.tolist(), reverse=False)

	## generate stacked bar plot
	bar_plot = viz.plot_clinsig_stacked_bar_by_gene(plot_df, col_clinsig, axis_order=axis_order)
	
	return plot_df[cols_df], bar_plot


def explore_clinsig_by_condition(df, col_clinsig, col_count, col_cond):
	## generate DF for plot
	plot_df = ss.generate_grouped_count_df(df, col_group=col_cond, col_count=col_count, col_clinsig=col_clinsig)

	## sort DF columns
	cols_df = ['Total'] + sort.sort_and_extract_clinsig(plot_df.columns.tolist(), reverse=False)

	## generate stacked bar plot
	bar_plot = viz.plot_clinsig_stacked_bar_by_condition(plot_df, col_clinsig)
	
	return plot_df[cols_df], bar_plot





################################################################################
#### Driver functions
################################################################################

def perform_explore_clinsig(df, col_id, col_clinsig, col_cond='condition', col_gene='gene'):
	## Variant clinsig: # of distinct variants per Variant classification
	print("\t\t.. performing Variant Clinical Significance analysis")
	clinsig_count_df, clinsig_count_plot = explore_clinsig(df, clinsig_col=col_clinsig, count_col=col_id)

	## RCV clinsig: # of NON-DISTINCT variants per RCV classification
	print("\t\t.. performing Variant-Condition (RCV) Clinical Significance analysis")
	rcv_clinsig_count_df, rcv_clinsig_count_plot = explore_clinsig(df, clinsig_col=col_clinsig+'.rcv', count_col=col_id)

	## grouped by GENE: # of distinct variants per clinsig classification, per GENE
	print("\t\t.. performing Gene-based analysis")
	gene_count_df, gene_plot = explore_clinsig_by_gene(df, col_count=col_id, col_clinsig=col_clinsig, col_gene=col_gene)

	## grouped by CONDITION: # of distinct variants per clinsig classification, per CONDITION
	print("\t\t.. performing Condition-based analysis")
	cond_count_df, cond_plot = explore_clinsig_by_condition(df, col_count=col_id, col_clinsig=col_clinsig, col_cond=col_cond)
	
	## assemble resulting DataFrames & Plotly plot figures
	print("\t\t.. assembling Clinical Significance analysis results")
	results_dict = dict(clinsig_var=dict(df=clinsig_count_df,
										 plot=clinsig_count_plot),
						clinsig_rcv=dict(df=rcv_clinsig_count_df,
										 plot=rcv_clinsig_count_plot),
						clinsig_var_gene=dict(df=gene_count_df,
											  plot=gene_plot),
						clinsig_var_cond=dict(df=cond_count_df,
											  plot=cond_plot))
	return results_dict



def perform_exploratory_analysis(result_dict, col_id, col_clinsig=COL_CLINSIG,
                                 cols_var=None, patho_all_cond=False):
	"""
	
	Args:
		result_dict:
		col_id:
		col_clinsig:
		cols_var:
		patho_all_cond:

	Returns:

	"""
	## prep DF for analysis
	print("\t.. converting DF for analyses & visualization")
	df = preprocess_visualization_df(result_dict['cv_full_df'], col_id=col_id,
	                                 col_clinsig=col_clinsig)
	
	## summarize dataset
	print("\t.. generate dataset summary")
	data_summary_df = get_dataset_summary(result_dict['input_df'], df, col_id,
	                                      col_cond='condition', col_gene='gene')
	
	## identify pathogenic variants
	print("\t.. identify pathogenic variants")
	patho_var_df, patho_var_detail_df = identify_patho_vars(result_dict['cv_var_summary_df'], result_dict['cv_full_df'], col_id, col_clinsig, cols_var, all_cond=patho_all_cond)
	
	## add new Data Viz & Pathogenic Variant DataFrames to result_dict
	print("\t.. add data viz & pathogenic variant DataFrames to result dictionary")
	result_dict.update(viz_df=df, data_summary_df=data_summary_df,
	                   patho_var_df=patho_var_df, patho_var_detail_df=patho_var_detail_df)

	## run Clinical Significance analyses
	print("\t.. starting Clinical Significance exploratory analyses")
	clinsig_result_dict = perform_explore_clinsig(df, col_id=col_id, col_clinsig=col_clinsig)
	result_dict.update(clinsig_result_dict)
	return result_dict



def run_clinvar_exploratory_analysis(var_file, out_dir, out_prefix, build, cols_var,
                                     cols_input=None, col_clinsig=COL_CLINSIG,
                                     write_files=True, write_plot_fxn=viz.write_plot_helper):
	"""
	
	Args:
		var_file:
		out_dir:
		out_prefix:
		build:
		cols_var:
		cols_input:
		col_clinsig:
		write_files:
		write_plot_fxn:

	Returns:

	"""
	## Step 1-3: run annotation_workflow
	annot_dict = annot.run_clinvar_annotation(var_file=var_file,
	                                            out_dir=out_dir,
	                                            out_prefix=out_prefix,
	                                            build=build,
	                                            cols_var=cols_var,
	                                            cols_input=cols_input,
	                                            write_output=False)
	## extract annotation workflow outputs
	result_dict = annot_dict['result_dict']
	_col_id = annot_dict['_col_id']
	_out_dir = annot_dict['_out_dir']
	
	## Step 4: run Exploratory analysis
	print("\n\nStep 4: run Exploratory analysis")
	result_dict = perform_exploratory_analysis(result_dict, col_id=_col_id,
	                                           col_clinsig=col_clinsig,
	                                           cols_var=cols_var)
	
	## Step 5: write output files
	if write_files:
		print("\n\nStep 5: write output files")
		write_output_exploratory_analysis(_out_dir, out_prefix, result_dict, write_plot_fxn=write_plot_fxn)
	
	return result_dict






