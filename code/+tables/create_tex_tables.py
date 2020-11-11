
# from ..EconTools import xlsx_to_latex
import sys
import os
import pandas as pd

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

def header_panel(filepath):
	df = pd.read_excel(filepath, index_col=0, header=0)
	tex = df.to_latex(float_format="%.1f")

	n = nlines(tex)
	lines_to_drop = [3, n-2, n-1]
	tex = drop_lines(tex, lines_to_drop)
	# tex = drop_lines(tex, 3)
	# tex = drop_line(tex, 0)
	# tex = drop_line(tex, -1)
	# tex = drop_line(tex, -1)

	return tex

def other_panel(dirpath, table, panel, panelname):
	filename = f'table{table}_panel{panel}.xlsx'
	filepath = os.path.join(dirpath, filename)

	df = pd.read_excel(filepath, index_col=0, header=0)
	tex = df.to_latex(float_format="%.1f")

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
	# tex = drop_line(tex, 3)
	# tex = drop_line(tex, -1)
	# tex = drop_line(tex, -1)
	# tex = drop_line(tex, -1)

	return tex

def save_tex_table_panels(dirpath):
	table1_header = header_panel(os.path.join(dirpath, 'table1_header.xlsx'))
	table1_panelA = other_panel(dirpath, 1, 'A', 'Income Statistics')
	table1_panelB = other_panel(dirpath, 1, 'B', 'Wealth Statistics')

	tex = '\n'.join([table1_header, table1_panelA, table1_panelB])
	tex += '\n\\end{tabular}'
	
	texfilepath = os.path.join(dirpath, 'table1.tex')
	fobj = open(texfilepath, 'w')
	fobj.write(tex)
	fobj.close()

if __name__ == '__main__':
	save_tex_table_panels(sys.argv[1])