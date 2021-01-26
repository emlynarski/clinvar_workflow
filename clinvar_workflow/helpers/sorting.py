# sorting.py

from clinvar_workflow.helpers import clinsig_sort_dict

## Pandas - setup
import pandas as pd
pd.options.mode.chained_assignment = None


################################################################################
#### clinsig column / list sorting functions
################################################################################
## custom sort list of clinsig category labels
def sort_clinsig(clinsig_unsorted_list, clinsig_dict=clinsig_sort_dict, reverse=False):
	"""

	Args:
		clinsig_unsorted_list:
		clinsig_dict:
		reverse:

	Returns:

	"""
	return sorted(clinsig_unsorted_list, reverse=reverse,
				  key=lambda x: clinsig_dict[x] if x in clinsig_dict else 24)


def sort_and_extract_clinsig(clinsig_unsorted_list, clinsig_dict=clinsig_sort_dict, reverse=False):
	"""

	Args:
		clinsig_unsorted_list:
		clinsig_dict:
		reverse:

	Returns:

	"""
	clinsig_unsorted_list2 = [c for c in clinsig_unsorted_list if c in clinsig_dict]
	return sorted(clinsig_unsorted_list2, reverse=reverse,
	              key=lambda x: clinsig_dict[x])


################################################################################
#### clinsig Dataframe sorting functions
################################################################################
## sort DF by clinsig column
def sort_clinsig_df_col(df, clinsig_col, clinsig_dict=clinsig_sort_dict, reverse=False):
	"""

	Args:
		df:
		clinsig_col:
		clinsig_dict:
		reverse:

	Returns:

	"""
	return df.loc[df[clinsig_col].map(clinsig_dict).sort_values(
		ascending=reverse).index].reset_index(drop=True).copy()


def sort_grouped_count_df_total(count_df, ascending=True):
	"""

	Args:
		count_df:
		ascending:

	Returns:

	"""
	## add sort columns, if necessary
	if 'Total' not in count_df.columns:
		count_df['Total'] = count_df.sum(axis=1)
	if 'idx' not in count_df.columns:
		count_df['idx'] = count_df.index.copy()
	
	## sort by 1) Total & 2) alphabetically
	count_df = count_df.sort_values(['Total', 'idx'], ascending=[not ascending, ascending])
	return count_df.drop(columns=['idx'])


def sort_grouped_count_df_alpha(count_df, ascending=True, sort_nots_top=False,
								nots_tuple=('not provided', 'not specified')):
	"""
	
	Args:
		count_df:
		ascending:
		sort_nots_top:
		nots_tuple:

	Returns:

	"""
	count_df = count_df.sort_index(axis=0, ascending=ascending)
	if sort_nots_top:
		## check for 'not provided' | 'not specified' --> sort accordingly
		grp_list = [grp for grp in count_df.index if grp not in nots_tuple]
		if len(grp_list) < len(count_df.index.tolist()):
			## sort grp_list --> generate sorting map
			grp_list = sorted(grp_list, reverse=ascending)
			sort_dict = {grp:i for i,grp in enumerate(list(nots_tuple) + grp_list)}
			## sort DF
			count_df['sort_order'] = count_df.index.map(sort_dict)
			count_df = count_df.sort_values(['sort_order'], ascending=[ascending])
			count_df.drop(columns=['sort_order'], inplace=True)
	return count_df



