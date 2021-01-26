#!/usr/bin/env python
# coding: utf-8
from __future__ import print_function
import os

## Pandas - setup
import pandas as pd
pd.options.mode.chained_assignment = None
pd.set_option('display.max_columns', None)

## Plotly data viz
from plotly import __version__
import plotly.io as pio
import plotly.graph_objs as go

## ClinVar workflow imports
# from clinvar_workflow.helpers import _COLOR_LINE, _COLOR_BG, clinsig_rgb_dict
from clinvar_workflow.helpers import clinsig_rgb_dict
from clinvar_workflow.helpers import sorting as sort
from clinvar_workflow.vizualization import _COLOR_LINE, _COLOR_BG



################################################################################
#### REQUIRED write plot file helper function
################################################################################

def write_plot_helper(fig, fig_file, out_dir, height=None, width=None):
	"""

	Args:
		fig:
		fig_file:
		out_dir:
		height:
		width:
	"""
	## specify image height
	if ('height' in fig.layout) & (fig.layout['height'] is not None):
		h = fig.layout['height']
	elif height is not None:
		h = height
	else:
		h = 450
	
	## specify image width
	if ('width' in fig.layout) & (fig.layout['width'] is not None):
		w = fig.layout['width']
	elif width is not None:
		w = width
	else:
		if fig.data[0]['type'] == 'bar':
			w = 650
		elif fig.data[0]['type'] == 'pie':
			w = 600
		else:
			w=600
	
	## write image files - pdf & png
	pio.write_image(fig, os.path.join(out_dir, 'png', fig_file+'.png'), height=h, width=w)
	pio.write_image(fig, os.path.join(out_dir, 'pdf', fig_file+'.pdf'), height=h, width=w)




################################################################################
#### clinical significance (clinsig) category color helper functions
################################################################################

## helper fxn to convert rgb string to rgba
def rgba_helper(rgb, alpha=1.0):
	"""

	Args:
		rgb:
		alpha:

	Returns:

	"""
	alpha_str = "," + str(alpha) + ")"
	rgba = rgb.replace('rgb', 'rgba').replace(')', alpha_str)
	return rgba


## bar plot: update clinsig trace color settings
def clinsig_color_bar(fig, clinsig_color_dict=clinsig_rgb_dict):
	"""

	Args:
		fig:
		clinsig_color_dict:
	"""
	for trace in fig.data:
		cs = trace['name']
		cs_marker = {'color':rgba_helper(clinsig_color_dict[cs], 0.999),
					 'line':{'width':0,
							 'color':rgba_helper(clinsig_color_dict[cs], 1.0)}}
		trace['marker'] = cs_marker


def clinsig_color_donut(fig, clinsig_color_dict=clinsig_rgb_dict):
	## remove count annot from clinsig labels (if present)
	labels = [i.split(' (')[0].strip() for i in fig['data'][0]['labels']]
	
	## update donut plot marker colors
	colors = [rgba_helper(clinsig_color_dict[l]) for l in labels]
	update_dict = dict(marker=dict(colors=colors,
								   line=dict(color='#FFFFFF', width=1)))
	fig.data[0].update(update_dict)



################################################################################
#### Data Viz: Plotly Clinical Significance by GROUP stacked bar plot functions
################################################################################

def plot_stacked_bar_horizontal(df, title, title_group, title_count, axis_order='total',
								bg_color=_COLOR_BG, line_color=_COLOR_LINE):
	#### generate Plotly traces
	y_values=df.index.tolist()
	traces = [go.Bar(name=c, y=y_values, x=df[c].tolist(), orientation='h')
					for c in df.columns if c.lower()!='total']

	#### setup plot layout
	## specify plot & formatting variables
	axis_title_font = dict(size=14)
	grp_max = pd.Series(df.index).apply(lambda x: len(x)).max()
	grp_offset = (grp_max*6 + 35)
	margin_l = grp_offset + 15
	margin_r = 30
	
	## replace y-axis title (cannot currently modify axis title position)
	_annot = [dict(showarrow=False, font=axis_title_font, borderpad=10,
				  xref='paper', yref='paper',
				  text='<b>' + title_group + '</b>', textangle=-90,
				  x=0, xanchor='center', xshift=-grp_offset+5,
				  y=0.5, yanchor='middle', yshift=20)]
	_axis_y = dict(title=dict(text=''), type='category',
				  showgrid=False, showline=True, linecolor=line_color,
				  showticklabels=True, dtick=1,
				  autorange=True, automargin=False)
	
	## format x-axis
	_axis_x = dict(title=dict(text='<b>' + title_count + '</b>',
							  font=axis_title_font),
				   showgrid=True, showline=True, linecolor=line_color,
				   tickcolor=line_color, ticks='outside',
				   automargin=False, mirror=False)
	
	## format Title
	_title = dict(text=title, pad=dict(b=45, t=15),
	              yanchor='top', y=0.99, yref='container')
	
	## format legend
	_legend = dict(traceorder='reversed', bgcolor=bg_color, itemsizing='constant',
	               x=1.02, xanchor='right',
	               y=0.5, yanchor='middle',
	               font=dict(size=10, color='black'))
	if axis_order!='total':
		_legend.update(dict(y=0.95, yanchor='auto'))
	
	## generate plot Layout dict
	layout_dict = dict(annotations=_annot,
					   yaxis=_axis_y,
					   xaxis=_axis_x,
					   title=_title,
					   bargap=0.2,
					   barmode='stack',
					   paper_bgcolor=bg_color,
					   plot_bgcolor=bg_color,
					   margin=dict(autoexpand=False, l=margin_l,
								   r=margin_r, b=85, t=70),
					   height=450+(len(df.index.tolist()) * 4),
					   width=margin_l+520+margin_r,
					   autosize=False,
					   legend=_legend)
	
	## generate Ploty Figure
	bar_plot = go.Figure(data=traces, layout=layout_dict)

	## update plot colors
	clinsig_color_bar(bar_plot)
	return bar_plot



def plot_stacked_bar_vertical(df, title, title_group, title_count, axis_order='total',
							  bg_color=_COLOR_BG, line_color=_COLOR_LINE):
	#### generate Plotly traces
	x_values = df.index.tolist()
	traces = [go.Bar(name=c, x=x_values, y=df[c].tolist(), orientation='v')
					for c in df.columns if c.lower()!='total']

	#### setup plot layout
	## specify plot & formatting variables
	axis_title_font = dict(size=14)
	grp_tickangle = -45
	grp_max = pd.Series(df.index).apply(lambda x: len(x)).max()
	grp_offset = grp_max * 9.5 * (-grp_tickangle/100) + 35
	margin_b = grp_offset + 20
	margin_r = 70
	
	## replace x-axis title (cannot currently modify axis title position)
	_annot = [dict(showarrow=False, font=axis_title_font, borderpad=10,
				  xref='paper', yref='paper',
				  text='<b>' + title_group + '</b>', textangle=0,
				  x=0.5, xanchor='center', xshift=0,
				  y=0, yanchor='middle', yshift=-grp_offset+15)]
	_axis_x = dict(title=dict(text=''), type='category',
				  showgrid=False, showline=True, linecolor=line_color,
				  showticklabels=True, dtick=1, tickangle=grp_tickangle,
				  autorange=True, automargin=False)
	
	## reformat y-axis
	_axis_y = dict(title=dict(text='<b>' + title_count + '</b>', font=axis_title_font),
				  showgrid=True, showline=True, linecolor=line_color,
				  tickcolor=line_color, ticks='outside', automargin=False)
	
	## format title
	_title = dict(text=title, yanchor='top', y=0.99, yref='container',
				  pad=dict(b=45, t=15))
	
	## format legend
	_legend = dict(traceorder='reversed', bgcolor=bg_color, itemsizing='constant',
				  x=1.02, xanchor='right', font=dict(size=10, color='black'))
	if axis_order!='total':
		_legend.update(dict(y=0.95, yanchor='auto'))
	
	## update plot formatting
	layout_dict = dict(annotations=_annot,
					   yaxis=_axis_y,
					   xaxis=_axis_x,
					   title=_title,
					   bargap=0.2,
					   barmode='stack',
					   paper_bgcolor=bg_color,
					   plot_bgcolor=bg_color,
					   margin=dict(autoexpand=False, l=70, t=70,
								   r=margin_r, b=margin_b),
					   height=470+margin_b-70,
					   width=900,
					   autosize=False,
					   legend=_legend)
	
	## generate Ploty Figure
	bar_plot = go.Figure(data=traces, layout=layout_dict)

	## update plot colors
	clinsig_color_bar(bar_plot)
	return bar_plot


def plot_clinsig_stacked_bar_by_gene(plot_df, col_clinsig, axis_order=None):
	if axis_order is None:
		axis_order='total'
	
	## set up plot & axis titles
	if 'rcv' in col_clinsig:
		count_title = '# of variants'
		clinsig_type = 'Variant-Condition (RCV)'
	else:
		count_title = '# of distinct variants'
		clinsig_type = 'Variant'
	plot_title = 'ClinVar ' + clinsig_type + ' Clinical Significance classifications<br>' + count_title + ', grouped by <i>Gene</i>'
	
	## filter top N genes for plot
	num_gene = 25
	if plot_df.shape[0] > num_gene:
		tmp_df = sort.sort_grouped_count_df_total(plot_df)
		plot_df = tmp_df.head(num_gene)
		
	## capitalize Clinsig columns
	plot_df.rename(columns={c:c.capitalize().replace('/ like', '/ Like')
	                        for c in plot_df.columns}, inplace=True)
	
	## extract non-zero clinsig columns & sort
	cols_keep = [c for c in plot_df.columns if plot_df[c].sum() > 0]
	cols_plot = sort.sort_and_extract_clinsig(cols_keep, reverse=True)
	
	## sort DF rows
	if axis_order == 'total':
		plot_df = sort.sort_grouped_count_df_total(plot_df)
	else:
		plot_df = sort.sort_grouped_count_df_alpha(plot_df)

	## call stacked bar function
	bar_plot = plot_stacked_bar_vertical(plot_df[cols_plot],
	                                     title_group='Gene',
										 title=plot_title,
										 title_count=count_title,
										 axis_order=axis_order)
	return bar_plot


def plot_clinsig_stacked_bar_by_condition(plot_df, col_clinsig):
	## set up plot & axis titles
	if 'rcv' in col_clinsig:
		clinsig_type = 'Variant-Condition (RCV)'
	else:
		clinsig_type = 'Variant'
	count_title = '# of variants'
	plot_title = 'ClinVar ' + clinsig_type + ' Clinical Significance classifications<br>' + count_title + ', grouped by <i>Condition</i>'
	
	## exclude 'not provided', 'not specified' rows
	plot_df = plot_df[~plot_df.index.isin(['not provided', 'not specified'])].copy()
	
	## filter top N cond for plot
	num_cond = 25
	if plot_df.shape[0] > num_cond:
		tmp_df = sort.sort_grouped_count_df_total(plot_df)
		plot_df = tmp_df.head(num_cond)
	
	## capitalize Clinsig columns
	plot_df.rename(columns={c:c.capitalize().replace('/ like', '/ Like')
	                        for c in plot_df.columns}, inplace=True)
	
	## extract non-zero clinsig columns & sort
	cols_keep = [c for c in plot_df.columns if plot_df[c].sum() > 0]
	cols_plot = sort.sort_and_extract_clinsig(cols_keep, reverse=True)
	
	## sort DF rows
	plot_df = sort.sort_grouped_count_df_total(plot_df, ascending=False)

	## call stacked bar function
	bar_plot = plot_stacked_bar_horizontal(plot_df[cols_plot],
	                                       title_group='Condition',
										   title=plot_title,
										   title_count=count_title,
										   axis_order='total')
	return bar_plot



################################################################################
#### Data Viz: Plotly Clinical Significance donut plot functions
################################################################################

def plot_donut(df, col_label, col_value, title, bg_color=_COLOR_BG):
	## generate Pie trace
	trace= go.Pie(labels=df[col_label],
				  values=df[col_value],
				  hole=.45,
				  textinfo='none',
				  hoverinfo='value+percent+text',
				  hovertext=df[col_label],
				  sort=False)
	
	## generate Plotly Figure with Pie trace
	fig_donut = go.Figure(data=[trace])

	## add center text to plot
	annot = [dict(text=str(df[col_value].sum()),
				  font=dict(size=36), showarrow=False,
				  x=0.5, y=0.5, xanchor='center')]
	
	## specify plot layout
	margin_r = 10 + (7 * df[df.columns[0]].apply(lambda x: len(x)).max())
	margin_l = 10
	margin_b = 20
	margin_t = 50
	layout = dict(annotations=annot,
				  paper_bgcolor=bg_color,
				  title=dict(text=title, xref='container',
							 x=0.5, xanchor='center'),
				  legend=dict(bgcolor=bg_color, itemsizing='constant',
							  x=1.15, tracegroupgap=5),
				  margin=dict(b=margin_b, l=margin_l, r=margin_r,
							  t=margin_t, autoexpand=False),
				  height=360,
				  width=margin_l + 330 + margin_r,
				  autosize=False)

	## update donut plot figure layout
	fig_donut.layout.update(layout)
	
	return fig_donut



def plot_donut_annot_legend(df, col_label, col_value, title, bg_color=_COLOR_BG):
	## add value to clinsig label str
	df = df.copy()
	df['label'] = df[col_label] + ' (' + df[col_value].astype(str) + ')'
	
	## generate Pie trace
	trace= go.Pie(labels=df['label'],
				  values=df[col_value],
				  hole=.45,
				  textinfo='none',
				  hoverinfo='value+percent+text',
				  hovertext=df[col_label],
				  sort=False)
	
	## generate Plotly Figure with Pie trace
	fig_donut = go.Figure(data=[trace])

	## add center text to plot
	annot = [dict(text=str(df[col_value].sum()),
				  font=dict(size=36), showarrow=False,
				  x=0.5, y=0.5, xanchor='center')]
	
	## specify plot layout
	margin_r = 35 + (7 * df[df.columns[0]].apply(lambda x: len(x)).max())
	margin_l = 10
	margin_b = 20
	margin_t = 50
	layout = dict(annotations=annot,
				  paper_bgcolor=bg_color,
				  title=dict(text=title, xref='container',
				             x=0.5, xanchor='center'),
				  legend=dict(bgcolor=bg_color, itemsizing='constant',
				              x=1.15, tracegroupgap=5),
				  margin=dict(b=margin_b, l=margin_l, r=margin_r, t=margin_t,
				              autoexpand=False),
				  height=360,
				  width=margin_l + 330 + margin_r,
				  autosize=False)

	## update donut plot figure layout
	fig_donut.layout.update(layout)
	
	return fig_donut

	

def plot_clinsig_donut(count_df, clinsig_col, plot_fxn=plot_donut_annot_legend):
	## dynamically set up plot title
	if 'rcv' in clinsig_col:
		donut_title = "RCV-level Clinical Significance classifications"
	else:
		donut_title = "Variant Clinical Significance classifications"

	## generate donut plot figure
	donut_plot_fig = plot_fxn(count_df, col_label=clinsig_col,
							  col_value=count_df.columns[1], title=donut_title)

	## update to clinsig colors
	clinsig_color_donut(donut_plot_fig)
	
	return donut_plot_fig


