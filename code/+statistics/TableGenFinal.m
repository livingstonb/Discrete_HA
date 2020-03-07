classdef TableGenFinal < statistics.TableGen
	properties
		default_fname = 'Final_Output.csv';
		included_names = {
			'Quarterly'
			'Annual'
		};
	end

	methods
		function obj = TableGenFinal(params, results, freq, use_all)
			obj.freq = freq;

			this_freq = find([params.freq] == freq);
			if nargin < 4
				use_all = false;
			end

			if use_all
				obj.selected_cases = this_freq;
			else
				all_names = [params.name];
				inames = [];
				for ii = 1:numel(obj.included_names)
					name = obj.included_names(ii);
					inames = [inames, find(ismember(all_names, name))];
				end

				obj.selected_cases = intersect(this_freq, inames);
			end

			obj.check_options(params, results);
			obj.outdir = params(1).outdir;
		end

		function output_table = create(obj, params, results, freq)
			output_table = table();
			if isempty(obj.selected_cases)
			    return;
			end

			shocks_labels = params(1).shocks_labels;

			for ip = obj.selected_cases
				p = params(ip);
				result_structure = results(ip);

				new_column = intro_panel(result_structure, p);

				temp = panel_A_Income(result_structure, p);
				new_column = [new_column; temp];

				temp = panel_B_Wealth(result_structure);
				new_column = [new_column; temp];

				temp = panel_C_MPCSize(result_structure, shocks_labels);
				new_column = [new_column; temp];

				temp = panel_D_MPCSign(result_structure, shocks_labels);
				new_column = [new_column; temp];

				column_label = sprintf('Specification%d', p.index);
				new_column.Properties.VariableNames = {column_label};
				output_table = [output_table, new_column];
			end

			obj.output = output_table;
		end
	end
end

function out = intro_panel(values, p, shocks_labels)
	out = table({p.name},...
		'VariableNames', {'results'},...
		'RowNames', {'Model'});

	new_labels = {	'Quarterly MPC (%)'
		            'Annual MPC (%)'
		            'Beta (Annualized)'
		};
	new_entries = {	values.direct.mpcs(5).avg_quarterly * 100
                    values.direct.mpcs(5).avg_annual * 100
                    values.direct.beta_annualized    
		};
	out = statistics.TableGen.append_to_table(out,...
		new_entries, new_labels);
end

function out = panel_A_Income(values, p)
	out = statistics.TableGen.new_table_with_header(...
		'Panel A: Income Statistics');

	new_labels = {	'Mean gross annual income'
		            'Stdev log annual gross income'
		            'Stdev log annual net income'
		};
	new_entries = {	values.direct.mean_grossy_A
                    values.direct.stdev_loggrossy_A
                    values.direct.stdev_lognety_A   
		};
	out = statistics.TableGen.append_to_table(out,...
		new_entries, new_labels);
end

function out = panel_B_Wealth(values)
	out = statistics.TableGen.new_table_with_header(...
		'Panel B: Wealth Statistics');

	% Mean assets and saving
	new_labels = {	'Mean assets, a'
					'Median a'
		};
	new_entries = {	values.direct.mean_a
					values.direct.median_a
		};
	out = statistics.TableGen.append_to_table(out, new_entries, new_labels);

	% Fraction with saving = 0
	new_labels = {'s == 0'};
	new_entries = values.direct.s0;
	out = statistics.TableGen.append_to_table(out, new_entries, new_labels);

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
	out = statistics.TableGen.append_to_table(out, new_entries, new_labels);

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
	out = statistics.TableGen.append_to_table(out, new_entries, new_labels);

	% Other stats
	new_labels = {	'Top 10% share'
					'Top 1% share'
					'Gini coefficient'
		};
	new_entries = {	values.direct.top10share
					values.direct.top1share
					values.direct.wealthgini
		};
	out = statistics.TableGen.append_to_table(out, new_entries, new_labels);
end

function out = panel_C_MPCSize(values, shocks_labels)
	out = statistics.TableGen.new_table_with_header(...
		'Panel C: MPC Size Effects');
	new_labels = {	sprintf('Annual MPC (%%), %s', shocks_labels{4})
		            sprintf('Annual MPC (%%), %s', shocks_labels{6})
		            sprintf('Quarterly MPC (%%), %s', shocks_labels{4})
		            sprintf('Quarterly MPC (%%), %s', shocks_labels{6})
		};
	new_entries = {	values.direct.mpcs(4).avg_annual * 100
					values.direct.mpcs(6).avg_annual * 100
					values.direct.mpcs(4).avg_quarterly * 100
					values.direct.mpcs(6).avg_quarterly * 100
		};
	out = statistics.TableGen.append_to_table(out, new_entries, new_labels);
end

function out = panel_D_MPCSign(values, shocks_labels)
	out = statistics.TableGen.new_table_with_header(...
		'Panel D: MPC Sign Effects');
	new_labels = {	sprintf('Annual MPC (%%), %s', shocks_labels{1})
					sprintf('Annual MPC (%%), %s', shocks_labels{2})
		            sprintf('Annual MPC (%%), %s', shocks_labels{3})
		            sprintf('Quarterly MPC (%%), %s', shocks_labels{1})
		            sprintf('Quarterly MPC (%%), %s', shocks_labels{2})
		            sprintf('Quarterly MPC (%%), %s', shocks_labels{3})
		};
	new_entries = {	values.direct.mpcs(1).avg_annual * 100
					values.direct.mpcs(2).avg_annual * 100
					values.direct.mpcs(3).avg_annual * 100
					values.direct.mpcs(1).avg_quarterly * 100
					values.direct.mpcs(2).avg_quarterly * 100
					values.direct.mpcs(3).avg_quarterly * 100
		};
	out = statistics.TableGen.append_to_table(out, new_entries, new_labels);
end