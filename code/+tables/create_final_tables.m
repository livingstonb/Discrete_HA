function [tables_out, table_gens] = create_final_tables(...
	params, results, decomps_baseline, save_tables)
	ntable = 1;
	table_gens{ntable} = tables.TableFinal_Baselines(...
		params, results, ntable);
	tables_out{ntable} = table_gens{ntable}.create(params, results);

	ntable = 2;
	table_gens{ntable} = tables.TableFinal_BaselineDecomps(...
		params, results, ntable);
	tables_out{ntable} = table_gens{ntable}.create(params, results);

	% Quarterly models
	ntable = 3;
	included_names = {	'Quarterly'
						'A/Y = 1'
						'A/Y = 0.5'
						'A/Y = 0.25'
						'r = 0% p.a.'
						'r = 5% p.a.'
						'No Death'
						'No Bequests'
						'Annuities'
						};
	table_gens{ntable} = tables.TableFinal_Experiments(...
		params, results, 3, included_names);
	tables_out{ntable} = table_gens{ntable}.create(params, results,...
		decomps_baseline);

	for itable = 1:numel(table_gens)
		if save_tables
			table_gens{itable}.save_table();
		end
	end
end