
# from ..EconTools import xlsx_to_latex
import sys
import os
import pandas as pd
import numpy as np
import re

def drop_lines(text, linenos):
	lines = text.splitlines()
	linenos.sort(reverse=True)

	for lineno in linenos:
		del lines[lineno]
	return '\n'.join(lines)

def replace_line(text, lineno, newline):
	lines = text.splitlines()
	lines[lineno] = newline
	return '\n'.join(lines)

def nlines(text):
	lines = text.splitlines()
	return len(lines)

def apply_float_formatting(df):
	for col in df.columns:
		if col != 'decimals':
			new_col = df.apply(lambda x: float2string(x, col), axis=1)
			df[col] = new_col

	del df['decimals']

	return df

def float2string(data, colname):
	precision = int(data['decimals'])
	try:
		if np.isnan(data[colname]):
			return ''
		else:
			return f'{data[colname]:.{int(data["decimals"])}f}'
	except:
		return data[colname]

def header_panel(filepath):
	df = pd.read_excel(filepath, index_col=0, header=0)
	df = apply_float_formatting(df)
	colfmt = 'l'
	for ic in range(len(df.columns)):
		colfmt += 'c'
	tex = df.to_latex(float_format="%.1f", escape=False, na_rep='', column_format=colfmt)

	n = nlines(tex)
	lines_to_drop = [3, n-2, n-1]
	tex = drop_lines(tex, lines_to_drop)
	tex = re.sub(r'__v\d+__', '', tex)

	return tex

def other_panel(dirpath, table, panel, panelname):
	filename = f'table{table}_panel{panel}.xlsx'
	filepath = os.path.join(dirpath, filename)

	df = pd.read_excel(filepath, index_col=0, header=0)
	df = apply_float_formatting(df)
	colfmt = 'l'
	for ic in range(len(df.columns)):
		colfmt += 'c'
	tex = df.to_latex(float_format="%.1f", escape=False, na_rep='', column_format=colfmt)

	cols = len(df.columns)
	newline = ''
	for i in range(cols):
		newline += ' & '

	newline += ' \\\\'

	n = nlines(tex)
	tex = replace_line(tex, 0, newline)
	tex = drop_lines(tex, [3, n-2, n-1])

	headername = f'Panel {panel}: {panelname}'
	newline = f'\\multicolumn{{{cols+1}}}{{c}}{{\\textbf{{{headername}}}}}\\\\'
	tex = replace_line(tex, 2, newline)
	tex = re.sub(r'__v\d+__', '', tex)

	return tex

def save_tex_table1_panels(dirpath):
	table1_header = header_panel(os.path.join(dirpath, 'table1_header.xlsx'))
	table1_panelA = other_panel(dirpath, 1, 'A', 'Income Statistics')
	table1_panelB = other_panel(dirpath, 1, 'B', 'Wealth Statistics')
	table1_panelC = other_panel(dirpath, 1, 'C', 'MPC Size Effects')
	table1_panelD = other_panel(dirpath, 1, 'D', 'MPC Sign Effects')

	tex = '\n'.join([
		table1_header,
		table1_panelA,
		table1_panelB,
		table1_panelC,
		table1_panelD,
		])
	tex += '\n\\bottomrule';
	tex += '\n\\end{tabular}'
	
	texfilepath = os.path.join(dirpath, 'table1.tex')
	fobj = open(texfilepath, 'w')
	fobj.write(tex)
	fobj.close()

def save_tex_table2_panels(dirpath):
	table2_header = header_panel(os.path.join(dirpath, 'table2_header.xlsx'))
	table2_panelA = other_panel(dirpath, 2, 'A', 'Decomposition of Mean MPC, around 0')
	table2_panelB = other_panel(dirpath, 2, 'B', 'Decomposition of Mean MPC, around 0.01')
	table2_panelC = other_panel(dirpath, 2, 'C', 'Decomposition of Mean MPC, around 0.05')
	table2_panelD = other_panel(dirpath, 2, 'D', 'Decomposition of Mean MPC - MPC$_{RA}$')

	tex = '\n'.join([
		table2_header,
		table2_panelA,
		table2_panelB,
		table2_panelC,
		table2_panelD,
		])
	tex += '\n\\bottomrule';
	tex += '\n\\end{tabular}'

	texfilepath = os.path.join(dirpath, 'table2.tex')
	fobj = open(texfilepath, 'w')
	fobj.write(tex)
	fobj.close()

def save_tex_experiment_table(dirpath, tableno):
	table_header = header_panel(os.path.join(dirpath, f'table{tableno}_header.xlsx'))
	table_panelA = other_panel(dirpath, tableno, 'A', 'Quarterly MPC Decomp w.r.t. Baseline')
	table_panelA2 = other_panel(dirpath, tableno, 'A2',
		'Quarterly MPC Decomp as \\% of $E[m_1]-E[m_0]$')
	table_panelB = other_panel(dirpath, tableno, 'B', 'Wealth Statistics')
	table_panelC = other_panel(dirpath, tableno, 'C', 'MPC Size Effects')
	table_panelD = other_panel(dirpath, tableno, 'D', 'MPC Sign Effects')
	tex = '\n'.join([
		table_header,
		table_panelA,
		table_panelA2,
		table_panelB,
		table_panelC,
		table_panelD,
		])
	tex += '\n\\bottomrule';
	tex += '\n\\end{tabular}'

	texfilepath = os.path.join(dirpath, f'table{tableno}.tex')
	fobj = open(texfilepath, 'w')
	fobj.write(tex)
	fobj.close()

if __name__ == '__main__':
	save_tex_table1_panels(sys.argv[1])
	save_tex_table2_panels(sys.argv[1])
	save_tex_experiment_table(sys.argv[1], 3)
	save_tex_experiment_table(sys.argv[1], 4)
	save_tex_experiment_table(sys.argv[1], 5)
	save_tex_experiment_table(sys.argv[1], 6)
	save_tex_experiment_table(sys.argv[1], 7)
	save_tex_experiment_table(sys.argv[1], 8)
	save_tex_experiment_table(sys.argv[1], 9)