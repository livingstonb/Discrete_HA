function tables = create_final_tables(params, results)
	freqencies = [1, 4]

	all_tables = cell(2, 1);
	for ifreq = 1:2
		freq = frequencies(ifreq);

		all_tables(ifreq,1) = tables.TableFinal_Main(...
			params, results, freq);
	end
end