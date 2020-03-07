function cell_out = cellround(cell_in, ndigits)
	cell_out = num2cell(cellfun(@(x) round(x, ndigits), cell_in));
end