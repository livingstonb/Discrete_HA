function out = wealth_panel(values, panel_name)
	if nargin < 2
		panel_name = 'Panel B: Wealth tables';
	end
	out = tables.TableGen.new_table_with_header(panel_name);

	% Mean assets and saving
	new_labels = {	'Mean assets, a'
					'Median a'
		};
	new_entries = {	values.direct.mean_a
					values.direct.median_a
		};
	out = tables.TableGen.append_to_table(out, new_entries, new_labels);

	% Fraction with saving = 0
	new_labels = {'s == 0'};
	new_entries = values.direct.s0;
	out = tables.TableGen.append_to_table(out, new_entries, new_labels);

	% Fraction with assets or cash<= some value
	new_labels = {	'a <= 0'
		            'a <= 0.5% mean ann inc'
		            'a <= 1% mean ann inc'
		            'a <= 2% mean ann inc'
		            'a <= 5% mean ann inc'
		            'a <= 10% mean ann inc'
		            'a <= 15% mean ann inc'
		            'a <= 1/6 own quarterly income'
		            'a <= 1/12 own quarterly income'
		};
	tmp = num2cell(values.direct.constrained(:));
	new_entries = {	values.direct.a_lt_sixth
	                values.direct.a_lt_twelfth
		};
	new_entries = [tmp; new_entries];
	out = tables.TableGen.append_to_table(out, new_entries, new_labels);

	% Percentiles
	new_labels = {	'10th percentile'
		            '25th percentile'
		            '50th percentile'
		            '75th percentile'
		            '90th percentile'
		            '95th percentile'
		            '99th percentile'
		            '99.9th percentile'
		};
	new_entries = num2cell(values.direct.wpercentiles(:));
	out = tables.TableGen.append_to_table(out, new_entries, new_labels);

	% Other stats
	new_labels = {	'Top 10% share'
					'Top 1% share'
					'Gini coefficient'
		};
	new_entries = {	values.direct.top10share
					values.direct.top1share
					values.direct.wealthgini
		};
	out = tables.TableGen.append_to_table(out, new_entries, new_labels);
end