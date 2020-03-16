classdef StatsTable < tables.BaseTable
	methods
		function output_table = create(obj, params, stats)
			obj.output = table();
			for ip = 1:obj.n_cols
				p_ip = params(ip);
				stats_ip = stats{ip};

				obj.current_column = table();
				obj.intro_stats_table(p_ip, stats_ip);

				% obj.income_stats_table();
				% obj.wealth_stats_table(stats_ip);
				% obj.mpc_size_table(stats_ip);
				% obj.mpc_sign_table(stats_ip);

				% obj.decomp_norisk_table(p_ip, stats_ip);

				% shock = stats_ip.mpcs(5).shock.value;
				% obj.mpc_comparison(stats_ip, shock);
				% obj.mpc_comparison_pct(stats_ip, shock);

				% if obj.one_asset_only
				% 	assets = {'b'};
				% else
				% 	assets = {'w', 'b', 'a'};
				% end
				% obj.percentiles_tables(stats_ip, assets);
    %             obj.illiquid_mpcs_table(stats_ip);
    %             obj.adj_costs_table(stats_ip);
    %             obj.other_params_table(stats_ip);
    %             obj.other_stats_table(stats_ip);

                obj.add_column(ip)
			end

			output_table = obj.output;
		end

		function intro_stats_table(obj, p, stats)
			out = table({p.label},...
				'VariableNames', {'results'},...
				'RowNames', {'Model'});

			new_entries = {
				stats.mpcs(5).quarterly
				stats.mpcs(5).annual
				stats.beta_A
				stats.std_log_gross_y_annual
				stats.std_log_net_y_annual
				stats.mean_gross_y_annual
			};

			obj.update_current_column(out, new_entries);
		end
	end
end