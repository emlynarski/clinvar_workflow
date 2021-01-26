# summary_stats.py

from clinvar_workflow.helpers import sorting as sort

## Pandas - setup
import pandas as pd
pd.options.mode.chained_assignment = None




################################################################################
#### count DF for exploratory analysis functions
################################################################################

## convert full variant DF --> count DF for plotting
def generate_count_df(df, cnt_col='', grp_cols=None):
	"""

	Args:
		df:
		cnt_col:
		grp_cols:

	Returns:

	"""
	## group DF --> get aggregate counts
	df2 = df.groupby(grp_cols)[cnt_col]\
			.agg([('nuniq', 'nunique')])\
			.reset_index()

	## reshape/pivot DF if plotting per gene | condition
	if len(grp_cols) > 1:
		df2 = df2.pivot_table(index=grp_cols[0], columns=grp_cols[1],
		                      values='nuniq', fill_value=0)
	return df2


def generate_grouped_count_df(df, col_group, col_clinsig, col_count):
	"""

	Args:
		df:
		col_group:
		col_clinsig:
		col_count:

	Returns:

	"""
	## generate initial count DF
	count_df = generate_count_df(df, cnt_col=col_count,
	                             grp_cols=[col_group, col_clinsig])
	
	## add Total column & sort rows
	count_df = sort.sort_grouped_count_df_total(count_df)
	
	## reorder Clinical Significance columns
	cols_clinsig_sort = sort.sort_and_extract_clinsig(count_df.columns.tolist(), reverse=False)
	return count_df[['Total'] + cols_clinsig_sort]





