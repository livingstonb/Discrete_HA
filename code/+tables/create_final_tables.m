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
	group = 'Q1';
	include_annual_mpcs = false;
	table_gens{ntable} = tables.TableFinal_Experiments(...
		params, results, ntable, group);
	tables_out{ntable} = table_gens{ntable}.create(params, results,...
		decomps_baseline, include_annual_mpcs);

	% Quarterly, beta heterogeneity
	ntable = 4;
	group = 'Q2';
	include_annual_mpcs = false;
	table_gens{ntable} = tables.TableFinal_Experiments(...
		params, results, ntable, group);
	tables_out{ntable} = table_gens{ntable}.create(params, results,...
		decomps_baseline, include_annual_mpcs);
	table_gens{ntable}.default_fname = ...
		sprintf('Table%d_beta_heterogeneity.csv', ntable);

	% Quarterly, CRRA experiments
	ntable = 5;
	group = 'Q3';
	include_annual_mpcs = false;
	table_gens{ntable} = tables.TableFinal_Experiments(...
		params, results, ntable, group);
	tables_out{ntable} = table_gens{ntable}.create(params, results,...
		decomps_baseline, include_annual_mpcs);
	table_gens{ntable}.default_fname = ...
		sprintf('Table%d_crra_experiments.csv', ntable);

	% Quarterly, EZ experiments
	ntable = 6;
	group = 'Q4';
	include_annual_mpcs = false;
	table_gens{ntable} = tables.TableFinal_Experiments(...
		params, results, ntable, group);
	tables_out{ntable} = table_gens{ntable}.create(params, results,...
		decomps_baseline, include_annual_mpcs);
	table_gens{ntable}.default_fname = ...
		sprintf('Table%d_ez_experiments.csv', ntable);

	% Quarterly, temptation experiments
	ntable = 7;
	group = 'Q5';
	include_annual_mpcs = false;
	table_gens{ntable} = tables.TableFinal_Experiments(...
		params, results, ntable, group);
	tables_out{ntable} = table_gens{ntable}.create(params, results,...
		decomps_baseline, include_annual_mpcs);
	table_gens{ntable}.default_fname = ...
		sprintf('Table%d_temptation_experiments.csv', ntable);

	% Quarterly, returns experiments
	ntable = 8;
	group = 'Q6';
	include_annual_mpcs = false;
	table_gens{ntable} = tables.TableFinal_Experiments(...
		params, results, ntable, group);
	tables_out{ntable} = table_gens{ntable}.create(params, results,...
		decomps_baseline, include_annual_mpcs);
	table_gens{ntable}.default_fname = ...
		sprintf('Table%d_returns_experiments.csv', ntable);

	% Mixed experiments
	ntable = 9;
	group = 'Q7';
	include_annual_mpcs = true;
	table_gens{ntable} = tables.TableFinal_Experiments(...
		params, results, ntable, group);
	tables_out{ntable} = table_gens{ntable}.create(params, results,...
		decomps_baseline, include_annual_mpcs);
	table_gens{ntable}.default_fname = ...
		sprintf('Table%d_mixed_experiments.csv', ntable);

	for itable = 1:numel(table_gens)
		if save_tables
			table_gens{itable}.save_table();
		end
	end
end