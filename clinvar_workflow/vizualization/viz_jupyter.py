#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function
import os

## Pandas - setup
import pandas as pd
pd.options.mode.chained_assignment = None
pd.set_option('display.max_columns', None)

## Plotly data viz
import plotly.graph_objs as go
from plotly.offline import init_notebook_mode
init_notebook_mode(connected=True)
## ipywidgets
from ipywidgets import HBox, VBox, AppLayout, Accordion, Layout




################################################################################
#### Data Viz: Plotly Table functions
################################################################################

def get_clinsig_table(count_df, clinsig_type):
	"""

	Args:
		count_df:
		clinsig_type:

	Returns:

	"""
	headerColor = 'grey'
	rowEvenColor = 'lightgrey'
	rowOddColor = 'white'
	col_clinsig = clinsig_type + ' classification'

	count_df = count_df.copy().rename(columns={count_df.columns[0]:col_clinsig})
	count_df = count_df[count_df[col_clinsig] != 'UNREPORTED']
	head_vals = ['<b>'+str(c)+'</b>' for c in count_df.columns]
	cell_vals = [count_df[c] for c in count_df.columns]

	fill_color=[[rowOddColor,rowEvenColor] * (int(count_df.shape[0]/2)+1)
	            for _ in cell_vals]
	height_h = 28
	height_c = 24
	height = height_h + (height_c * len(cell_vals[0])) + 15
	
	table = go.Table(
		header=dict(values=head_vals,
					fill=dict(color=headerColor),
					line=dict(width=0),
					align=['left','center'],
					font=dict(color='white', size=12),
					height=height_h),
		cells=dict(values=cell_vals,
				   fill=dict(color=fill_color),
				   line=dict(width=0),
				   align=['left', 'center'],
				   font=dict(size = 12),
				   height=height_c)
	)
	go_fig = go.Figure(data=[table])
	go_fig.layout.update({'height':height, 'margin':{'t':5, 'b':0, 'l':0, 'autoexpand':False}})
	return go_fig


# def get_clinsig_table_color(count_df, clinsig_type, clinsig_color_dict=clinsig_rgb_dict):
def get_clinsig_table_color(count_df, clinsig_type, clinsig_color_dict):
	## rename DF columns for Table & drop UNREPORTED
	col_clinsig = clinsig_type + ' classification'
	count_df = count_df.copy().rename(columns={count_df.columns[0]:col_clinsig})
	count_df = count_df[count_df[col_clinsig].str.upper() != 'UNREPORTED']

	nrow = count_df.shape[0]

	## set up Table colors
	clinsig = count_df[count_df.columns[0]].tolist()
	clinsig_colors = [clinsig_color_dict[c] for c in clinsig]
	white = ['#FFFFFF'] * nrow

	## specify Cell values
	cell_vals = [[''] * nrow] + [clinsig, count_df[count_df.columns[1]].tolist()]

	## specify heights
	height_c = 20
	margin = 20
	height = (height_c *(nrow+1)) + 2*margin

	## create Plotly Table trace
	table = go.Table(
		columnwidth=[1, 10, 2],
		header=dict(height=0, fill=dict(color=None), line=dict(width=0)),
		cells=dict(values=cell_vals,
				   fill=dict(color=[clinsig_colors, white, white]),
				   line=dict(width=[8, 0, 0], color=[white, white, white]),
				   align='left',
				   font=dict(size=12, color='black'),
				   height=height_c)
	)
	## create Plotly Plot Figure with Table trace
	go_fig = go.Figure(data=[table])
	go_fig.layout.update({'height':height,
						  'width':350,
						  'autosize':True,
						  'margin':{'t':0, 'b':margin, 'l':0, 'autoexpand':True}})
	return go_fig



def get_grouped_clinsig_table(count_df):
	## select non-zero columns & reset index
	cols_keep = [c for c in count_df.columns if count_df[c].sum() > 0]
	count_df = count_df[cols_keep].copy().reset_index(drop=False)

	## rename DF columns for Table
	r_dict = {'gene':'Gene symbol<br>', 'condition':'Condition name<br>',
			  'Total':'Total<br>', 'Pathogenic':'Pathogenic<br>',
			  'Pathogenic / Likely pathogenic':'Pathogenic<br>/ Likely<br>pathogenic',
			  'Likely pathogenic':'Likely<br>pathogenic',
			  'drug response':'Drug<br>response', 'risk factor':'Risk<br>factor',
			  'Affects':'Affects<br>', 'association':'Association<br>',
			  'other':'Other<br>', 'protective':'Protective<br>',
			  'Conflicting':'Conflicting<br>',
			  'Benign':'Benign<br>', 'Likely benign':'Likely<br>benign',
			  'Benign / Likely benign':'Benign<br>/ Likely<br>benign',
			  'Uncertain significance':'Uncertain<br>significance'}
	count_df.rename(columns=r_dict, inplace=True)

	## specify Cell & Header data
	cell_vals = [count_df[c] for c in count_df.columns]
	head_vals = ['<b>'+str(c)+'</b>' for c in count_df.columns]

	## set up Table colors
	headerColor = 'grey'
	rowEvenColor = 'lightgrey'
	rowOddColor = 'white'
	fill_color=[[rowOddColor,rowEvenColor] * (int(count_df.shape[0]/2)+1)
	            for _ in cell_vals]

	## set column widths
	if 'Gene' in count_df.columns[0]:
		width_group = 100
	else:
		width_group = 230
	width_data = [min((10+max(len(s) for s in i.split('<br>'))*6), 80)
	              for i in count_df.columns[1:]]
	col_widths = [width_group] + width_data
	width = sum(col_widths) + 30

	## create Plotly Table trace
	table = go.Table(
		columnwidth=col_widths,
		header=dict(values=head_vals,
					fill=dict(color=headerColor),
					line=dict(width=0),
					align=['left','center'],
					font=dict(color='white', size=10)),
		cells=dict(values=cell_vals,
				   fill=dict(color=fill_color),
				   line=dict(width=0),
				   align=['left', 'center'],
				   font=dict(size = 11))
	)

	## create Plotly Plot Figure with Table trace
	go_fig = go.Figure(data=[table])
	go_fig.layout.update({'height':230,
						  'width':width,
						  'margin':{'t':10, 'b':0, 'l':0, 'r':25, 'autoexpand':True},
						  'autosize':True})
	return go_fig



################################################################################
#### Data Viz: Jupyter + ipywidgets + Plotly figure widgets
################################################################################

def get_donut_plot_figure_widget(plot_fig, hide_textinfo=True):
	"""

	Args:
		plot_fig:
		hide_textinfo:

	Returns:

	"""
	w = plot_fig.layout.width
	r = plot_fig.layout.margin['r']

	fig_widget = go.FigureWidget(plot_fig)
	fig_widget.layout.update({'title':{'text':''},
							  'margin': {'autoexpand': False, 'b': 0,
							             'l': 0, 't': 0, 'r':r+20},
						   'autosize':True,
						   'width':w+20})
	if hide_textinfo:
		fig_widget.data[0].update({'textinfo':'none'})
	return fig_widget


def get_clinsig_plot_annot_legend_figure(cs_plot):
	## reformat plot layout
	h = cs_plot.layout.height
	
	cs_plot2 = go.Figure(cs_plot)
	cs_plot2.layout.update(dict(margin=dict(autoexpand=False, t=10, b=10),
								autosize=False,
								height=h - 40,
								title=dict(text='', pad=dict(b=0, t=0))))
	
	## convert into FigureWidgets
	cs_plot_widget = go.FigureWidget(go.Figure(cs_plot2))
	
	## generate ipywidget
	plot_annot_fig = HBox([cs_plot_widget])
	return plot_annot_fig


def get_bar_plot_figure_widget(plot_fig):
	## reformat plot layout
	h = plot_fig.layout.height
	margin_b = plot_fig.layout.margin['b']
	margin_r = plot_fig.layout.margin['r']
	
	plot_fig2 = go.Figure(plot_fig)
	plot_fig2.layout.update(dict(margin=dict(autoexpand=True, t=0,
											 b=margin_b-20, r=margin_r+30),
								 autosize=True,
								 height=h-(40 + 10),
								 title=dict(text='', pad=dict(b=0, t=0))))
	
	## convert Plotly Figure to FigureWidget
	fig_widget = go.FigureWidget(plot_fig2)
	return fig_widget



################################################################################
#### Plotly combined Table & Plot figure fxns
################################################################################

def get_clinsig_plot_table_figure(count_df, cs_plot, clinsig_type='Variant',
								  table_fxn=get_clinsig_table):
	## generate clinsig table
	cs_table = table_fxn(count_df.copy(), clinsig_type)

	## reformat copy of plot
	r = cs_plot.layout.margin['r']
	w = cs_plot.layout.width
	h = cs_plot.layout.height
	
	fig = go.Figure(cs_plot)
	fig.data[0].update(dict(showlegend=False))
	fig.layout.update(dict(margin=dict(autoexpand=False, t=10, b=10, r=10),
						   autosize=False,
						   width=w - r + 10,
						   height=h - 40,
						   title=dict(text='', pad=dict(b=0, t=0))))
	
	## convert into FigureWidgets
	cs_plot_widget = go.FigureWidget(fig)
	cs_table_widget = go.FigureWidget(cs_table)
	
	## setup Widget Layout
	height = max(cs_table.layout.height, cs_plot_widget.layout.height)
	width = cs_table.layout.width + cs_plot_widget.layout.width + 50
	_layout = Layout(
		display='flex',
		flex_flow='row',
		justify_content='space-around',
		height=str(height)+'px',
		width=str(width)+'px')
	
	## generate ipywidget
	plot_table_fig = HBox([cs_plot_widget, cs_table_widget], layout=_layout)
	
	##############################
	#### @TODO: REMOVE ME AFTER DEV
	tmp_return = dict(table=cs_table, table_widget=cs_table_widget, plot_widget=cs_plot_widget)
	##############################
	# return plot_table_fig
	return plot_table_fig, tmp_return


def get_grouped_clinsig_plot_table_figure(count_df, count_plot):

	## convert Plotly stacked bar Figure to FigureWidget
	cs_plot_widget = get_bar_plot_figure_widget(go.Figure(count_plot))

	## generate Plotly Table Figure & convert to FigureWidget
	cs_table = get_grouped_clinsig_table(count_df.copy())
	cs_table_widget = go.FigureWidget(cs_table)
	cs_table_widget.layout.update(dict(autosize=True, margin=dict(autoexpand=True)))

	curr_width = cs_table_widget.layout.width
	
	plot_table_fig = AppLayout(header=None,
							  left_sidebar=None,
							  center=HBox([cs_plot_widget],
										  layout=Layout(min_width=str(cs_plot_widget.layout.width)+'px',
														overflow='scroll hidden')),
							  right_sidebar=None,
							  footer=HBox([cs_table_widget],
										  layout=Layout(overflow='scroll',
														min_width=str(curr_width+20)+'px')),
							  pane_widths=["0px", 10, "0px"],
							  pane_heights=["0px", 3, "250px"])

	##############################
	#### @TODO: REMOVE ME AFTER DEV
	tmp_return = dict(table=cs_table, table_widget=cs_table_widget, plot_widget=cs_plot_widget)
	##############################
	# return plot_table_fig
	return plot_table_fig, tmp_return


def display_css(widget, css_class):
	widget.add_class(css_class)
	return widget


def display_clinsig_exploratory_analysis(result_dict):
	## generate Variant Clinical Significance Plotly Table
	print("\t\t.. generating Variant Clinical Significance Plotly Table")
	clinsig_var = get_clinsig_plot_annot_legend_figure(result_dict['clinsig_var']['plot'])

	## generate Variant-Condition (RCV) Clinical Significance Plotly Table
	print("\t\t.. generating Variant-Condition (RCV) Clinical Significance Plotly Table")
	clinsig_rcv = get_clinsig_plot_annot_legend_figure(result_dict['clinsig_rcv']['plot'])

	## generate Gene Variant Clinical Significance Plotly Table
	print("\t\t.. generating Gene Variant Clinical Significance Plotly Table")
	gene = get_grouped_clinsig_plot_table_figure(result_dict['clinsig_var_gene']['df'],
												 result_dict['clinsig_var_gene']['plot'])

	## generate Condition Variant Clinical Significance Plotly Table
	print("\t\t.. generating Condition Variant Clinical Significance Plotly Table")
	cond = get_grouped_clinsig_plot_table_figure(result_dict['clinsig_var_cond']['df'],
												 result_dict['clinsig_var_cond']['plot'])
	
	## generate containers to display widgets
	print("\t\t.. generating containers to display widgets")
	clinsig_fig_list = [clinsig_var, clinsig_rcv, gene[0], cond[0]]
	acc_titles = ['Variant Clinical Significance',
				  'Variant-Condition (RCV) Clinical Significance',
				  'Variant Clinical Significance, grouped by Gene',
				  'Variant Clinical Significance, grouped by Condition']

	## generate list of Accordion containers
	acc_list = [Accordion(children=[plot_table_widget]) for plot_table_widget in clinsig_fig_list]
	for i in range(len(acc_list)):
		acc_list[i].set_title(0, acc_titles[i])
	
	print("\t\t.. assembling results to display")
	## assemble side-by-side donut plots
	donut_acc = HBox([acc_list[0], acc_list[1]])
	
	## assemble results to display
	display_acc = VBox([donut_acc, acc_list[2], acc_list[3]])
	
	##############################
	#### @TODO: REMOVE ME AFTER DEV
	cs_keys = ['clinsig_var', 'clinsig_rcv', 'clinsig_var_gene', 'clinsig_var_cond']
	viz_list = [clinsig_var, clinsig_rcv, gene[0], cond[0]]
	for k,v in zip(cs_keys, viz_list):
		result_dict[k]['plot_table_widget'] = v
	##############################
	
	return result_dict, display_acc




