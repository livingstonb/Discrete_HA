classdef TableFinal_BaselineDecomps < tables.TableGen
	properties
		default_fname = '';
		included_names = {
			'Quarterly'
			'Annual'
		};
	end

	methods
		function obj = TableFinal_BaselineDecomps(...
			params, results, table_num, use_all)
            if nargin < 4
                use_all = false;
            end

            frequencies = [1, 4];
			obj = obj@tables.TableGen(...
				params, results, frequencies, use_all);

			obj.default_fname = sprintf(...
            	'Table%d_baseline_decompositions.csv', table_num);
		end

		function output_table = create(obj, params, results)
			output_table = table();
			if isempty(obj.selected_cases)
			    return;
			end

			shocks_labels = params(1).shocks_labels;

			for ip = obj.selected_cases
				p = params(ip);
				result_structure = results(ip).direct;

				new_column = tables.DecompPanels.decomp_intro(...
					result_structure, p);

				decomp_structure = results(ip).decomp_norisk;
				prefixes = {'A', 'B', 'C'};
				for ithresh = 1:numel(params(1).abars)
					panel_prefix = strcat('Panel', " ", prefixes{ithresh});
					temp = tables.DecompPanels.mean_mpc_decomp(...
						decomp_structure, p, ithresh, panel_prefix);
					new_column = [new_column; temp];
				end

				panel_prefix = 'Panel D';
				decomp_structure = results(ip).decomp_RA;
				temp = tables.DecompPanels.ra_mpc_decomp(...
						decomp_structure, panel_prefix);
				new_column = [new_column; temp];

				column_label = sprintf('Specification%d', p.index);
				new_column.Properties.VariableNames = {column_label};
				output_table = [output_table, new_column];
			end

			obj.output = output_table;
		end
	end
end