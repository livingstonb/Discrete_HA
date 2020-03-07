function [tables_out, table_gens] = create_final_tables(...
	params, results, save_tables)
	ntable = 1;
	table_gens{ntable} = tables.TableFinal_Main(...
		params, results, ntable);

	ntable = 2;
	table_gens{ntable} = tables.TableFinal_BaselineDecomps(...
		params, results, ntable);

	for itable = 1:numel(table_gens)
		tables_out{itable} = table_gens{itable}.create(params, results);

		if save_tables
			table_gens{itable}.save_table();
		end
	end
end