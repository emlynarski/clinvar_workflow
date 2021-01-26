# clinvar_query.py

## Pandas - setup
import pandas as pd
pd.options.mode.chained_assignment = None
pd.set_option('display.max_columns', None)


## ClinVar query API - myvariant - setup
import myvariant
mv = myvariant.MyVariantInfo()

#TODO: change set() --> set literal --> remove warnings

################################################################################
#### ClinVar MyVariant variables
################################################################################

##----MyVariant ClinVar fields------------------------------------------------##
COL_CLINSIG = 'clinical_significance'
COL_COND = 'conditions.name'
COL_RCV = 'accession'

## integer columns
cv_int_fields = ['number_submitters', 'variant_id', 'hg19.start', 'hg19.end', 'hg38.start', 'hg38.end']

## ClinVar fields to drop:
cv_drop_fields = ['_license', 'allele_id', 'gene.id', 'cytogenic']


def myvariant_run_clinvar_query(input_var_df, col_hgvs, genome_build):
	"""

	Args:
		input_var_df:
		col_hgvs:
		genome_build:

	Returns:

	"""
	## use MyVariant API to query ClinVar
	raw_query_df = mv.getvariants(input_var_df[col_hgvs], fields='clinvar',
								  assembly=genome_build, as_dataframe=1)
		
	## check if all variants were 'notfound'
	if ('notfound' in raw_query_df.columns) & (len(raw_query_df.columns)==1):
		print("\tALERT: NONE of the input variants were found in ClinVar database.")
		return None
	
	## check if some variants were 'notfound'
	if 'notfound' in raw_query_df.columns:
		raw_query_df.drop('notfound', axis=1, inplace=True)
		raw_query_df['_id'] = raw_query_df.index.copy()
	
	## extract reported variants
	if 'clinvar' not in raw_query_df.columns:
		cv_raw_df = raw_query_df[raw_query_df['clinvar.variant_id'].notnull()].copy()
		cv_raw_df.rename(
			columns={c: c.replace('clinvar.', '') for c in raw_query_df.columns},
			inplace=True)
	else:
		## unpack reported variant 'clinvar' fields
		cv_raw_df = raw_query_df[raw_query_df['clinvar'].notnull()].copy()
		cv_raw_df = pd.concat([cv_raw_df.drop(['clinvar'], axis=1),
							   pd.DataFrame(cv_raw_df['clinvar'].values.tolist(),
											index=cv_raw_df.index)], axis=1)
		## unpack other nested columns
		nested_cols = ['gene', 'hg19', 'hg38', 'hgvs']
		for col in nested_cols:
			if col in cv_raw_df.columns:
				cv_raw_df = extract_nested_dict(cv_raw_df, col)
	
	## extract reported variants
	_drop_cv = [c for c in cv_drop_fields if c in cv_raw_df.columns]
	cv_raw_df.drop(_drop_cv, axis=1, inplace=True)
	
	return cv_raw_df


################################################################################
#### ClinVar query postprocessing functions
################################################################################

def extract_nested_dict(df, col_nested, prefix=None):
	"""

	Args:
		df:
		col_nested:
		prefix:

	Returns:

	"""
	## extract nested fields (if they exist)
	mask = df[col_nested].astype('str').str.startswith('{')
	nested = pd.DataFrame(df[mask][col_nested].values.tolist(), index=df[mask].index)
	
	## rename nested extracted columns
	if prefix is None:
		prefix = col_nested+'.'
	nested.rename(columns={c:prefix+c for c in nested.columns}, inplace=True)
	
	## add extracted columns to nested rows
	extract_df = pd.concat([df[mask].drop([col_nested], axis=1), nested],
						   axis=1, sort=False)

	## add unnested rows (if there are any)
	df2 = pd.concat([extract_df, df[~mask].drop(col_nested, axis=1)], sort=False)
	return df2


def extract_nested_list(df, col_nested):
	"""

	Args:
		df:
		col_nested:

	Returns:

	"""
	## extract list items to separate rows
	melt = pd.DataFrame(df[col_nested].copy().values.tolist(), index=df.index)\
			.merge(df.copy(), right_index=True, left_index=True)\
			.melt(id_vars=df.columns.tolist(), value_name=col_nested+'2')\
			.drop(['variable', col_nested], axis=1)\
			.dropna(subset=[col_nested+'2'])\
			.reset_index(drop=True)

	## extract nested fields to separate columns
	extract_df = pd.concat([melt.drop([col_nested + '2'], axis=1),
							pd.DataFrame(melt[col_nested + '2'].values.tolist(),
										 index=melt.index)], axis=1, sort=False)
	return extract_df


def extract_conditions(df):
	"""

	Args:
		df:

	Returns:

	"""
	## extract records with a single conditions entry
	single = df[df['conditions'].astype('str').str.startswith('{')].copy()
	cond_single = pd.concat([single.drop('conditions', axis=1),
							 pd.DataFrame(single['conditions'].values.tolist(),
										  index=single.index)], axis=1, sort=False)

	## extract records with multiple conditions
	multi = df[df['conditions'].astype('str').str.startswith('[')].copy()
	cond_multi = extract_nested_list(multi, 'conditions')
	
	## combine single + multi condition entries
	cond = pd.concat([cond_single, cond_multi], sort=False).reset_index(drop=True)
	
	## rename conditions columns
	cond.rename(columns={c:'conditions.'+c for c in cond.columns if c not in df.columns}, inplace=True)
	
	## extract identifiers (if they exist)
	iden = extract_nested_dict(cond, 'conditions.identifiers', 'id.')
	
	## sort
	iden = iden.sort_values(['_id', 'accession']).reset_index(drop=True)
	return iden


def extract_nested_rcv(cv_raw_df):
	"""

	Args:
		cv_raw_df:

	Returns:

	"""
	## select variants with multiple RCV accessions (have list of RCV dicts in 'clinvar.rcv' column)
	multi_df = cv_raw_df[cv_raw_df['rcv'].astype('str').str.startswith('[')].copy()
	## multi-rcv variants: drop any 'rcv.*' columns
	multi_drop = [c for c in multi_df.columns if 'rcv.' in c]
	multi_df.drop(multi_drop, axis=1, inplace=True)
	
	## select variants with a single NESTED RCV accession (have rcv dict in 'clinvar.rcv' column)
	nested_df = cv_raw_df[cv_raw_df['rcv'].astype('str').str.startswith('{')].copy()
	## nested single-rcv variants: drop any 'conditions.*' columns & rename 'rcv.*' columns
	nested_drop = [c for c in nested_df.columns if 'conditions.' in c]
	nested_df = nested_df.drop(nested_drop, axis=1)\
				.rename(columns={c:c.replace('rcv.','') for c in nested_df.columns})\
				.reset_index(drop=True)
	
	## extract RCV fields --> combine
	rcv_multi = extract_nested_list(multi_df, 'rcv')
	rcv_nested = extract_nested_dict(nested_df, 'rcv', '')
	rcv_df = pd.concat([rcv_multi, rcv_nested], sort=False).reset_index(drop=True)
	return rcv_df


def extract_single_rcv_conditions(cv_raw_df, cols_multi_rcv):
	"""

	Args:
		cv_raw_df:
		cols_multi_rcv:

	Returns:

	"""
	## variants with a single UNNESTED RCV accession ('clinvar.rcv' column is NaN)
	single_df = cv_raw_df[cv_raw_df['rcv'].isnull()].copy()
	
	## single-rcv variants with empty rcv field: drop 'rcv' column
	single_df = single_df.drop(['rcv'], axis=1).reset_index(drop=True)

	## single-rcv variants with empty rcv field: rename 'rcv.*' field columns
	rename_dict = {c: c.replace('rcv.', '').replace('conditions.identifiers.', 'id.') for
				   c in single_df.columns}
	single_df.rename(columns=rename_dict, inplace=True)
	
	## extract single-rcv unnested variants
	if 'conditions' in single_df.columns:
		cond_mask = single_df['conditions'].notnull()
		single_df_cond = extract_conditions(single_df[cond_mask][cols_multi_rcv])

		## combine with unnested 'conditions' variants
		single_df = pd.concat([single_df[~cond_mask].drop('conditions', axis=1),
							   single_df_cond], sort=False)
	return single_df


def add_rcv_condition_name_helper(cv_df, col_id, col_rcv, col_cond):
	"""

	Args:
		cv_df:
		col_id:
		col_rcv:
		col_cond:

	Returns:

	"""
	rcv_cond_df = cv_df[[col_id, col_rcv, col_cond]].copy().dropna()\
					.drop_duplicates()\
					.groupby([col_id, col_rcv])[col_cond]\
					.agg([(col_cond+'.rcv', lambda x: '; '.join(sorted(list(x))))])\
					.reset_index()
	cv_df = cv_df.merge(rcv_cond_df, on=[col_id, col_rcv], how='outer')
	return cv_df


def extract_rcv_clinsig_list(cv_df, col_clinsig):
	tmp = cv_df[cv_df[col_clinsig].fillna('').str.contains(',')]

	if tmp.shape[0] > 0:
		tmp[col_clinsig] = tmp[col_clinsig].apply(lambda x: x.split(', '))
		tmp = tmp.explode(col_clinsig)
		cv_df = pd.concat([cv_df[~cv_df[col_clinsig].str.contains(',')], tmp ])
	return cv_df


def cast_int_str(df, col):
	df.loc[df[col].notnull(), col] = df.loc[df[col].notnull(), col].astype(int).astype(str)
	df.loc[df[col].isnull(), col] = ''
	return df
	

def int_col_cast_helper(df, cols_int):
	int_cols = [c for c in df.columns if c in cols_int]
	for c in int_cols:
		cast_int_str(df, c)
	return df


def myvariant_clinvar_rcv_data_wrangling(cv_raw_df, col_id, col_clinsig, cols_int):
	"""

	Args:
		cv_raw_df:
		col_id:
		col_clinsig:
		cols_int:

	Returns:

	"""
	## extract and combine multi-RCV & single-RCV variants
	rcv_extracted = extract_nested_rcv(cv_raw_df)
	rcv_extracted2 = extract_conditions(rcv_extracted)
	single_df = extract_single_rcv_conditions(cv_raw_df, rcv_extracted.columns)
	cv_df = pd.concat([rcv_extracted2, single_df], sort=False)
	
	## extract RCV clinical significance lists
	cv_df = extract_rcv_clinsig_list(cv_df, col_clinsig)
	
	## cast float cols --> int --> cast int cols to str
	cv_df = int_col_cast_helper(cv_df, cols_int)

	## cast list containing columns to str
	cols_list = ['hgvs.coding', 'hgvs.genomic']
	cv_df[cols_list] = cv_df[cols_list].astype(str).fillna('')

	## add ClinVar status column
	cv_df['clinvar_status'] = 'reported'

	## rename '_id' & 'clinical_significance' column
	cv_df.rename(columns={'_id':col_id, col_clinsig:col_clinsig+'.rcv'},
				 inplace=True)
	
	## add the RCV condition name as 'conditions.name.rcv'
	cv_df = add_rcv_condition_name_helper(cv_df, col_id, 'accession',
										  'conditions.name')
	return cv_df



################################################################################
#### Helper functions
################################################################################

def agg_uniq_sorted_list_dropna(df, grp_col, agg_col):
	"""

	Args:
		df:
		grp_col:
		agg_col:

	Returns:

	"""
	return \
	df[[grp_col, agg_col]].copy().dropna().astype(str).drop_duplicates()\
		.groupby(grp_col, as_index=True)[agg_col]\
		.agg(lambda x: sorted(list(x)) if x.nunique() > 1 else x)


##----convert chromosome column to integer for proper sorting-----------------##
def chrN_df_helper(c, chr_map):
	"""

	Args:
		c:
		chr_map:

	Returns:

	"""
	if c.strip().lower() in chr_map:
		return chr_map[c.strip().lower()]
	return 30


def add_chrN_sort_cols(df, var_cols, chr_map={'x':24, 'y':25, 'mt':26}):
	"""

	Args:
		df:
		var_cols:
		chr_map:

	Returns:

	"""
	## copy original chrom & position column before multiindexing
	df['sort_chr'] = df[var_cols[0]].copy().astype('str')
	df['sort_pos'] = df[var_cols[1]].copy()
	df['sort_date'] = df['last_evaluated'].copy().fillna('')
	df['sort_asc'] = df['accession'].copy().fillna('')

	## map chrom str to int: 'X':24, 'Y':25, 'MT':26
	mask_chr = df['sort_chr'].str.isnumeric()
	df.loc[mask_chr, 'sort_chrN'] = df.loc[mask_chr,'sort_chr'].copy().astype('int')
	df.loc[~mask_chr, 'sort_chrN'] = df['sort_chr'].apply(lambda x: chrN_df_helper(x, chr_map))
	return df



################################################################################
#### Generate variant ClinVar summary from MyVariant RCV-level data functions
################################################################################

def variant_summary_non_clinsig_fields(cv_df, col_id, col_clinsig, cols_non_rcv):
	"""

	Args:
		cv_df:
		col_id:
		col_clinsig:
		cols_non_rcv:

	Returns:

	"""
	######################################################################
	#### #TODO: make this dynamic!!!
	cols_collapse_record = ['last_evaluated', 'accession', 'review_status', 'number_submitters']
	######################################################################
	
	## collapse non-RCV columns
	summ_df = cv_df[[col_id]+cols_non_rcv].copy().drop_duplicates().reset_index(drop=True)
	summ_df.set_index(col_id, drop=False, inplace=True)

	## collapse conditions.names
	summ_df['conditions.name'] = agg_uniq_sorted_list_dropna(cv_df, col_id, 'conditions.name')
	summ_df['conditions.name.set'] = cv_df.groupby(col_id, as_index=True)['conditions.name'].agg(set)

	## collapse conditions.synonyms
	summ_df['conditions.synonyms'] = cv_df[[col_id, 'conditions.synonyms']].copy()\
										.dropna()\
										.astype(str)\
										.drop_duplicates()\
										.groupby(col_id, as_index=True)['conditions.synonyms']\
										.agg(lambda x: ', '.join(sorted(list(x))))

	## dynamically collapse identifier columns
	for col in cv_df.columns:
		if ('id.' in col) & (col!=col_id):
			summ_df[col] = agg_uniq_sorted_list_dropna(cv_df, col_id, col)

	## sort record columns by date --> collapse
	collapse_record_tmp = cv_df[[col_id]+ cols_collapse_record].copy()\
								.drop_duplicates()\
								.reset_index(drop=True)
	collapse_record_tmp = collapse_record_tmp.fillna(value={'last_evaluated':''})\
								.sort_values([col_id, 'last_evaluated'],
											 ascending=[True, False])
	summ_df[cols_collapse_record] = collapse_record_tmp.groupby(col_id, as_index=True)\
														.agg(list)
	
	## collapse any remaining columns
	cols_to_collapse = [c for c in cv_df.columns if
						(c not in [col_clinsig, 'conditions.name.rcv']) & (
									c not in summ_df.columns)]
	for col in cols_to_collapse:
		summ_df[col] = agg_uniq_sorted_list_dropna(cv_df, col_id, col)
	return summ_df


##----clinsig classification rules--------------------------------------------##
def classify_clinsig_rule1_single_rcv(cv_df, col_id, col_clinsig, rule_num=1):
	"""

	Args:
		cv_df:
		col_id:
		col_clinsig:
		rule_num:

	Returns:

	"""
	classified = cv_df.groupby(col_id)\
					.filter(lambda x: x['accession'].nunique()==1)\
					[[col_id, col_clinsig]].drop_duplicates()\
					.reset_index(drop=True)
	classified['rule'] = rule_num
	return classified


def classify_clinsig_rule2_single_distinct(cv_df, col_id, col_clinsig, rule_num=2):
	"""

	Args:
		cv_df:
		col_id:
		col_clinsig:
		rule_num:

	Returns:

	"""
	## rule 2: multi-RCV variants with 1 distinct RCV clinsig -> clinsig = distinct RCV clinsig
	classified = cv_df.groupby(col_id)\
						.filter(lambda x: x[col_clinsig].nunique()==1)\
						[[col_id, col_clinsig]].drop_duplicates()\
						.reset_index(drop=True)
	classified['rule'] = rule_num
	return classified

def classify_clinsig_rule3_not_provided_distinct(cv_df, col_id, col_clinsig, rule_num=3):
	"""

	Args:
		cv_df:
		col_id:
		col_clinsig:
		rule_num:

	Returns:

	"""
	## rule 3: variants with 'not provided' + 1 distinct RCV clinsig -> clinsig = distinct RCV clinsig
	classified = cv_df[(cv_df[col_clinsig]!='not provided')]\
						.groupby(col_id)\
						.filter(lambda x: x[col_clinsig].nunique()==1)\
						[[col_id, col_clinsig]].copy().drop_duplicates()\
						.reset_index(drop=True)
	classified['rule'] = rule_num
	return classified


def classify_clinsig_rule4_expert_review(cv_df, col_id, col_clinsig, rule_num=4):
	## select variants reviewed by expert panel
	cv_df_expert = cv_df[cv_df['review_status'] == 'reviewed by expert panel']
	
	## variants with 1 distinct expert reviewed clinical significance
	expert = cv_df_expert.groupby(col_id)\
					.filter(lambda x: x[col_clinsig].nunique()==1)\
					[[col_id, col_clinsig]].copy().drop_duplicates()
	
	## variants with >1 distinct expert reviewed clinical significance ==> conflicting
	conflict = cv_df_expert[~cv_df_expert[col_id].isin(expert[col_id])]\
					[[col_id]].drop_duplicates()
	conflict[col_clinsig] = 'Conflicting interpretations of pathogenicity'

	classified = pd.concat([expert, conflict], sort=False)\
					.reset_index(drop=True)
	classified['rule'] = rule_num
	return classified


def classify_clinsig_rule5_patho_rcv(cv_df, col_id, col_clinsig, rule_num=5):
	## set up pathogenic variables
	patho_list = ['Pathogenic', 'Pathogenic/Likely pathogenic', 'Likely pathogenic']
	patho_set = set(['pathogenic', 'drug response', 'risk factor', 'affects', 'association', 'other'])
	lp_set = set(['pathogenic', 'likely pathogenic'])
	p_lp_set = lp_set | set(['pathogenic/likely pathogenic', 'drug response'])

	## select all variants with any patho RCV classifications
	p_any = cv_df[cv_df[col_clinsig]!='not provided'].copy()\
				.groupby(col_id)\
				.filter(lambda x: any(x[col_clinsig].isin(patho_list))).copy()\
				.reset_index(drop=True)
	p_any['clinsig'] = p_any[col_clinsig].copy().str.lower().str.strip()

	##--------------------------------------------------------##
	## PART 1: drop rows with 'no assertion' review status --> agg set{clinsig}  --> classify
	agg1 = p_any.loc[~(p_any['review_status'].str.startswith('no assertion'))].copy()\
					.groupby(col_id)['clinsig']\
					.agg([('_set', set)]).reset_index(drop=False)

	## identify & classify "Likely pathogenic" variants
	agg1['_likely'] = agg1['_set'].apply(lambda x: x == set(['likely pathogenic']))
	agg1.loc[agg1['_likely'], col_clinsig] = "Likely pathogenic"
	
	## identify & classify "Pathogenic" variants
	agg1['_patho'] = agg1['_set'].apply(lambda x: x.issubset(patho_set) &
	                                              x.issuperset(set(['pathogenic'])))
	agg1.loc[agg1['_patho'], col_clinsig] = "Pathogenic"
	
	## extract classified variants
	_classified1 = agg1.loc[agg1[col_clinsig].notnull(), [col_id, col_clinsig]]
	
	##--------------------------------------------------------##
	## PART 2: remove classified p_any vars --> agg set{clinsig} --> classify
	agg2 = p_any[~p_any[col_id].isin(_classified1[col_id])].copy()\
					.groupby(col_id)['clinsig']\
					.agg([('_set', set)]).reset_index(drop=False)

	## identify & classify {'Pathogenic', 'Likely pathogenic'} variants
	agg2['_p+Lp'] = agg2['_set'].apply(lambda x: x == lp_set)
	agg2.loc[agg2['_p+Lp'], col_clinsig] = "Likely pathogenic"
	
	## identify & classify 'Pathogenic/Likely pathogenic' variants
	agg2['_p/Lp'] = agg2['_set'].apply(lambda x: x.issubset(p_lp_set) & len(x)>0)
	agg2.loc[(agg2['_p/Lp'] & (~agg2['_p+Lp'])), col_clinsig] \
				= "Pathogenic/Likely pathogenic"
	
	## concat classified pathogenic vars --> classify remaining as Conflicting
	classified = pd.concat([_classified1, agg2[[col_id, col_clinsig]]\
					.fillna('Conflicting interpretations of pathogenicity')])
	classified['rule'] = rule_num
	return classified


def classify_clinsig_rule6_likely_benign(cv_df, col_id, col_clinsig, rule_num=6):
	"""

	Args:
		cv_df:
		col_id:
		col_clinsig:
		rule_num:

	Returns:

	"""
	benign_list = list(['Benign', 'Likely benign', 'Benign/Likely benign'])
	benign_set = set(benign_list)
	
	## select all variants with any benign RCV classifications --> agg set{clinsig}
	tmp = cv_df[cv_df[col_clinsig]!='not provided'].copy()\
			.groupby(col_id)\
			.filter(lambda x: any(x[col_clinsig].isin(benign_list)))\
			.groupby(col_id)\
			[col_clinsig].agg([('_set', set)])\
			.reset_index(drop=False)
	
	## identify possible 'Benign/Likely benign' variants
	tmp['_likely'] = tmp['_set'].apply(lambda x: x.issubset(benign_set) & len(x)>0)
	benign_likely_vars = tmp[tmp['_likely']][col_id].copy()
	
	## TEST for 'Benign/Likely benign' EDGE CASES
	def likely_edge_case_grp(g):
		return g.groupby(col_clinsig).apply(
			lambda x: any(x['conditions.name'] == 'not provided') & (
				x['accession'].nunique() == 1))

	edge_case_vars = cv_df[cv_df[col_id].isin(benign_likely_vars)].copy()\
						.groupby(col_id)\
						.filter(lambda x: any(likely_edge_case_grp(x)))\
						[col_id].unique().tolist()
	
	## extract & classify'Benign/Likely benign' variants
	vars_benign = [v for v in benign_likely_vars if v not in edge_case_vars]
	tmp.loc[tmp[col_id].isin(vars_benign), col_clinsig] = "Benign/Likely benign"
	classified = tmp[[col_id, col_clinsig]]\
					.fillna('Conflicting interpretations of pathogenicity')
	classified['rule'] = rule_num
	return classified


def classify_clinsig_ruleX_default_conflicting(cv_df, col_id, col_clinsig, rule_num):
	"""

	Args:
		cv_df:
		col_id:
		col_clinsig:
		rule_num:

	Returns:

	"""
	## rule 6) classify any remaining variants as 'Conflicting'
	classified = cv_df[[col_id]].drop_duplicates()\
					.reset_index(drop=True)
	classified[col_clinsig] = 'Conflicting interpretations of pathogenicity'
	classified['rule'] = rule_num
	return classified
##----------------------------------------------------------------------------##

def variant_summary_classify_rcv_clinsig(cv_df, col_id, col_clinsig):
	## clinsig classification rules helper functions
	rule_fxn_list = [(classify_clinsig_rule1_single_rcv, 1),
					 (classify_clinsig_rule2_single_distinct, 2),
					 (classify_clinsig_rule3_not_provided_distinct, 3),
					 (classify_clinsig_rule4_expert_review, 4),
					 (classify_clinsig_rule5_patho_rcv, 5),
					 (classify_clinsig_rule6_likely_benign, 6),
					 (classify_clinsig_ruleX_default_conflicting, 7)]
	
	var_list, rule_df_list = [], []
	cols_classify = [col_id, col_clinsig, 'conditions.name', 'accession',
	                 'review_status', 'rsid', 'variant_id']
	
	cv_df['review_status'] = cv_df['review_status'].str.strip().str.lower()
	rcv_tmp = cv_df.loc[cv_df['review_status']!='no assertion provided',
	                    cols_classify].copy()\
					.drop_duplicates()\
					.reset_index(drop=True)
	
	## apply clinical significance classification rules
	for rule_fxn, rule_num in rule_fxn_list:
		rule_df = rule_fxn(rcv_tmp[~rcv_tmp[col_id].isin(var_list)],
						   col_id, col_clinsig, rule_num)
		rule_df_list.append(rule_df)
		var_list = var_list + rule_df[col_id].tolist()
	classified = pd.concat(rule_df_list, sort=False).reset_index(drop=True)

	## rename clinsig column
	col_clinsig2 = col_clinsig.rsplit('.rcv', maxsplit=1)[0]
	classified.rename(columns={col_clinsig: col_clinsig2}, inplace=True)

	## aggregate RCV clinsig as set
	rcv_clinsig_set = cv_df[[col_id, col_clinsig]].copy().drop_duplicates()\
							.sort_values([col_id, col_clinsig])\
							.groupby(col_id)[col_clinsig]\
							.agg([(col_clinsig2 + '.rcv.set', lambda x: set(x))])

	clinsig_classified = classified.merge(rcv_clinsig_set, on=col_id, how='outer')
	clinsig_classified = clinsig_classified.sort_values(col_id)\
											.set_index(col_id, drop=True)
	return clinsig_classified


def generate_clinvar_variant_summary_df(cv_df, col_id, col_clinsig, cols_non_rcv):
	"""

	Args:
		cv_df:
		col_id:
		col_clinsig:
		cols_non_rcv:

	Returns:

	"""
	## classify variant's clinical significance
	var_clinsig = variant_summary_classify_rcv_clinsig(cv_df, col_id, col_clinsig)
	
	## collapse/summarize remaining columns
	var_summ_tmp = variant_summary_non_clinsig_fields(cv_df, col_id, col_clinsig,
													  cols_non_rcv=cols_non_rcv)
	
	## merge variant summary
	var_summ_df = var_clinsig.drop(['rule'], axis=1).merge(var_summ_tmp,
	                                                       left_index=True,
	                                                       right_index=True,
	                                                       how='outer')
	
	###################################
	#### #TODO: reorder columns
	###################################
	
	return var_summ_df






################################################################################
#### Process ClinVar data functions
################################################################################


def flag_condition_clinsig_conflicts_helper(cv_df, col_id, col_sig, col_cond, col_rcv):
	"""

	Args:
		cv_df:
		col_id:
		col_sig:
		col_cond:
		col_rcv:

	Returns:

	"""
	filt_df = cv_df.groupby([col_id, col_cond])\
					.filter(lambda x: x[col_sig].nunique() > 1)\
					.reset_index(drop=True)
	if filt_df.shape[0] > 0:
		conflict_df = filt_df.groupby([col_id, col_cond, col_sig])[col_rcv]\
						.agg(lambda x: ', '.join(sorted(list(x.unique()))))\
							.reset_index(name="rcv_str")\
						.groupby([col_id, col_cond])\
							.apply(lambda x: dict(zip(x[col_sig],x['rcv_str'])))\
							.reset_index(name='dict_clinsig')\
						.groupby([col_id])\
							.apply(lambda x: dict(zip(x[col_cond], x['dict_clinsig'])))\
							.reset_index(name='FLAG.condition_conflict.dict')
	
		conflict_df['FLAG.condition_conflict'] = True
		conflict_df = conflict_df[[col_id, 'FLAG.condition_conflict', 'FLAG.condition_conflict.dict']]
	else:
		conflict_df = pd.DataFrame(columns=[col_id, 'FLAG.condition_conflict', 'FLAG.condition_conflict.dict'])
	return conflict_df


def flag_condition_duplicated_helper(cv_df, col_id, col_cond, col_rcv, col_rcv_cond):
	"""

	Args:
		cv_df:
		col_id:
		col_cond:
		col_rcv:
		col_rcv_cond:

	Returns:

	"""
	## extract subset of ClinVar DF
	tmp_df = cv_df[[col_id, col_rcv, col_cond, col_rcv_cond]].copy().dropna().drop_duplicates()
	
	## grpby id, rcv -> agg dict(rcv:rcv cond) => grpby id, cond -> agg & merge {rcv:rcv cond} dicts => merge
	agg_df = tmp_df.merge(
					tmp_df.groupby([col_id, col_rcv])\
						.apply(lambda x: dict(zip(x[col_rcv],x[col_rcv_cond])))\
						.reset_index(name='rcv_cond_dict'),
					on=[col_id, col_rcv], how='left'
					).groupby([col_id, col_cond])\
						.apply(lambda x: {k: v for d in list(x['rcv_cond_dict']) for k, v in d.items()})\
						.reset_index(name='cond_dict')

	## identify duplicated conditions -> grpby cond -> agg dict(cond: {rcv:rcv cond})
	dup_cond_tmp = agg_df[agg_df['cond_dict'].apply(lambda x: len(x) > 1)]
	
	if dup_cond_tmp.shape[0] > 0:
		flag_df = dup_cond_tmp.groupby(col_id)\
							.apply(lambda x: dict(zip(x[col_cond],x['cond_dict'])))\
							.reset_index(name='FLAG.condition_duplicated.dict')
		
		## add Bool FLAG col
		flag_df['FLAG.condition_duplicated'] = True
		flag_df = flag_df[[col_id, 'FLAG.condition_duplicated', 'FLAG.condition_duplicated.dict']]
	else:
		flag_df = pd.DataFrame(columns=[col_id, 'FLAG.condition_duplicated', 'FLAG.condition_duplicated.dict'])
	return flag_df


def add_agg_stats_summary_df(cv_df, cv_var_summary_df, col_id, col_clinsig, col_cond, col_rcv):
	"""

	Args:
		cv_df:
		cv_var_summary_df:
		col_id:
		col_clinsig:
		col_cond:
		col_rcv:

	Returns:

	"""
	## copy subset of ClinVar full DF for aggregations
	cv_tmp_df = cv_df[[col_id, col_clinsig, col_rcv, col_cond]].copy().drop_duplicates()

	cv_var_agg = cv_tmp_df.groupby(col_id).agg({col_cond:[('cond.nuniq', 'nunique')],
												col_rcv:[('rcv.nuniq', 'nunique')],
												col_clinsig:[('clinsig.nuniq', 'nunique')]})
	cv_var_agg.columns = cv_var_agg.columns.droplevel()
	cv_var_agg = cv_var_agg[sorted(cv_var_agg.columns.tolist())]
	
	
	## identify 'Pathogenic' variants --> agg count and Boolean indicator columns
	tmp_patho_agg = cv_tmp_df[cv_tmp_df[col_clinsig]=='Pathogenic'].copy()\
							.drop_duplicates()\
							.groupby(col_id)[col_cond]\
							.agg({('patho_cond.nuniq', 'nunique'), ('patho_cond.set', set)})\
							.reset_index().set_index(col_id, drop=True)
	cv_var_agg = cv_var_agg.merge(tmp_patho_agg, left_index=True, right_index=True, how='left')
	cv_var_agg[['patho_cond.nuniq']] = cv_var_agg[['patho_cond.nuniq']].fillna(0).astype(int)
	cv_var_agg['patho_cond.%'] = cv_var_agg['patho_cond.nuniq'] / cv_var_agg['cond.nuniq']
	
	## add Boolean indicator columns
	cv_var_agg['_.multi_rcv'] = cv_var_agg['rcv.nuniq'] > 1
	cv_var_agg['_.multi_cond'] = cv_var_agg['cond.nuniq'] > 1
	cv_var_agg['_.multi_clinsig'] = cv_var_agg['clinsig.nuniq'] > 1
	
	patho_vars = cv_var_summary_df[cv_var_summary_df[
									   col_clinsig.rsplit('.rcv', maxsplit=1)[
										   0]].str.strip().str.lower() == 'pathogenic'][
		col_id].unique().tolist()
	cv_var_agg['_.patho_ALL_cond'] = (cv_var_agg['patho_cond.nuniq'] == cv_var_agg[
		'cond.nuniq']) | (cv_var_agg.index.isin(patho_vars))
	cv_var_agg['_.patho_ANY_cond'] = (cv_var_agg['patho_cond.nuniq'] > 0)
	
	## merge with summary DF
	cv_var_summary_df = cv_var_summary_df.merge(cv_var_agg, on=col_id, left_index=False,
												right_index=True, how='outer')
	
	## add FLAGS to variant summary
	flag_cond_dup = flag_condition_duplicated_helper(cv_df, col_id, col_cond, col_rcv,
													 col_cond + '.rcv')
	flag_conflict = flag_condition_clinsig_conflicts_helper(cv_df, col_id, col_clinsig,
															col_cond, col_rcv)
	cv_var_summary_df = cv_var_summary_df.merge(
		flag_cond_dup.merge(flag_conflict, on=col_id, how='outer'), on=col_id, how='left')
	
	## FLAG & Boolean indicator columns - fillna:
	cols_flag = [c for c in
				 flag_cond_dup.columns.tolist() + flag_conflict.columns.tolist() if
				 ('FLAG.' in c) & ('.dict' not in c) & ('.list' not in c)]
	cols_bool = [c for c in cv_var_agg.columns if c.startswith('_.')]
	cv_var_summary_df[cols_flag + cols_bool] = cv_var_summary_df[
		cols_flag + cols_bool].fillna(False)
	return cv_var_summary_df


## add variant summary DF columns to full ClinVar DF
def add_clinvar_variant_summary_columns(cv_summary_df, cv_df, col_id, col_clinsig,
										cols_cv, cols_var, cols_input):
	"""

	Args:
		cv_summary_df:
		cv_df:
		col_id:
		col_clinsig:
		cols_cv:
		cols_var:
		cols_input:

	Returns:

	"""
	## extract summary DF columns to add to new full reported DF
	###################### @TODO: make dynamic!
	cols_clinsig_agg = ['clinical_significance.rcv.set', 'clinsig.nuniq']
	cols_misc = ['rcv.nuniq', 'cond.nuniq', 'conditions.name.set', 'patho_cond.nuniq',
				 'patho_cond.%', 'patho_cond.set']
	##################################################################
	
	cols_bool = [c for c in cv_summary_df.columns if c.startswith('_.')]
	cols_flag = [c for c in cv_summary_df.columns if 'FLAG.' in c]
	cols_to_add = [col_clinsig] + cols_var + cols_input + cols_clinsig_agg + cols_misc + cols_bool + cols_flag + ['clinvar_status']
	
	summ_to_add_df = cv_summary_df[[col_id]+cols_to_add].copy()
	
	## merge summary DF columns with new full DF
	cols_drop = [c for c in cv_df.columns if c in cols_to_add]
	cv_full_df = cv_df.copy().drop(columns=cols_drop)\
					.merge(summ_to_add_df, on=col_id, how='outer')

	## reorder full DF columns
	cols_reorder_agg = [c for c in cv_df.columns if c not in cols_cv]
	cols_reorder_tmp = [col_id, col_clinsig, col_clinsig+'.rcv'] + cols_drop
	cols_cv2 = [c for c in cv_df.columns if (c in cols_cv) and (c not in cols_reorder_tmp)]
	reorder_full_cols = cols_var + cols_reorder_tmp + cols_cv2 + cols_clinsig_agg + cols_misc + cols_reorder_agg + cols_bool + cols_flag + cols_input
	return cv_full_df[reorder_full_cols]




def clinvar_summary_add_input_cols(cv_summary_df, input_var_df, col_id, col_clinsig, cols_var, cols_input):
	"""

	Args:
		cv_summary_df:
		input_var_df:
		col_id:
		col_clinsig:
		cols_var:
		cols_input:

	Returns:

	"""
	cv_var_df = input_var_df[cols_var + cols_input + [col_id]].merge(
		cv_summary_df.copy().reset_index(drop=True), on=col_id, how='outer')

	cv_var_df['clinvar_status'].fillna('NOT in ClinVar', inplace=True)
	cv_var_df[[col_clinsig]] = cv_var_df[[col_clinsig]].fillna('UNREPORTED', axis=1)
	
	cv_var_df.loc[cv_var_df[col_clinsig] == 'UNREPORTED', col_clinsig + '.rcv.set'] = \
	cv_var_df[col_clinsig].apply(lambda x: set([x]))

	## add columns for sorting full variant ClinVar DataFrame
	if 'sort_chrN' not in cv_var_df.columns:
		cv_var_df = add_chrN_sort_cols(cv_var_df, cols_var)
	
	## sort full variant ClinVar DataFrame
	sort_cols = [c for c in cv_var_df.columns if 'sort_' in c]
	cv_var_df = cv_var_df.sort_values(['sort_chrN', 'sort_pos'],
									  ascending=[True, True]).reset_index(drop=True)
	cv_var_df.drop(sort_cols, axis=1, inplace=True)
	return cv_var_df



################################################################################
#### Driver functions
################################################################################

def run_clinvar_query(input_var_df, build, col_id, col_clinsig=COL_CLINSIG,
                      cols_int=cv_int_fields):
	"""

	Args:
		input_var_df:
		build:
		col_id:
		col_clinsig:
		cols_int:

	Returns:

	"""
	## run ClinVar query
	print("\t.. run ClinVar query")
	cv_raw_df = myvariant_run_clinvar_query(input_var_df, col_hgvs=col_id, genome_build=build)
	
	if cv_raw_df is None:
		return None
	
	## ClinVar query data wrangling
	print("\t.. ClinVar query data wrangling")
	cv_df = myvariant_clinvar_rcv_data_wrangling(cv_raw_df, col_id, col_clinsig, cols_int)
	return cv_df



def process_clinvar_query(cv_df, input_var_df, cols_var, cols_input, col_id,
						  col_clinsig=COL_CLINSIG, col_rcvclinsig=COL_CLINSIG+'.rcv',
						  col_cond=COL_COND, col_rcv=COL_RCV, col_gene='gene.symbol'):
	"""

	Args:
		cv_df:
		input_var_df:
		cols_var:
		cols_input:
		col_id:
		col_clinsig:
		col_rcvclinsig:
		col_cond:
		col_rcv:
		col_gene:

	Returns:

	"""
	## single value per variant columns
	###################### @TODO: MAKE DYNAMIC!!!
	cols_collapse_non_rcv = ['clinvar_status', 'preferred_name', 'variant_id', 'type', 'rsid', 'gene.symbol', 'chrom', 'ref', 'alt', 'hg19.start', 'hg19.end', 'hg38.start', 'hg38.end', 'hgvs.coding', 'hgvs.genomic']
	############################################
	cols_cv_extracted = cv_df.columns.tolist()
	
	## generate variant ClinVar summary (1 row per variant)
	print("\t.. generating variant ClinVar summary (1 row per variant)")
	cv_var_summary_df = generate_clinvar_variant_summary_df(cv_df, col_id=col_id, col_clinsig=col_rcvclinsig, cols_non_rcv=cols_collapse_non_rcv)

	## add input columns to variant summary DF --> update UNREPORTED variants
	print("\t.. adding input columns to variant summary DF --> update UNREPORTED variants")
	cv_var_summary_df = clinvar_summary_add_input_cols(cv_var_summary_df, input_var_df, col_id=col_id, col_clinsig=col_clinsig, cols_var=cols_var, cols_input=cols_input)

	## add aggregation stats, Boolean indicator columns & FLAG columns to CV variant summary DF
	print("\t.. adding aggregation stats, Boolean indicator columns & FLAG columns to CV variant summary DF")
	cv_var_summary_df = add_agg_stats_summary_df(cv_df, cv_var_summary_df, col_id=col_id, col_clinsig=col_rcvclinsig, col_cond=col_cond, col_rcv=col_rcv)

	## add variant summary DF columns to full ClinVar DF
	print("\t.. adding variant summary DF columns to full ClinVar DF")
	cv_full_df = add_clinvar_variant_summary_columns(cv_var_summary_df, cv_df, col_id=col_id, col_clinsig=col_clinsig, cols_cv=cols_cv_extracted, cols_var=cols_var, cols_input=cols_input)
	
	return dict(cv_var_summary_df=cv_var_summary_df, cv_full_df=cv_full_df, input_df=input_var_df)

	
	


